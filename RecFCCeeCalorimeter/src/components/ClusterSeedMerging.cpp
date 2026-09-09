#include "ClusterSeedMerging.h"

#include "edm4hep/Cluster.h"
#include "edm4hep/MutableCluster.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// ============================================================
//  ClusterSeedMerging
// ============================================================

ClusterSeedMerging::ClusterSeedMerging(const std::string& name, ISvcLocator* svcLoc)
    : ClusterSeedingBase(name, svcLoc,
                         {KeyValues("CaloDrivenSeeds", {}), // allow multiple input collections to be merged together
                          KeyValues("TrackDrivenSeeds", {})},
                         {KeyValue("OutputMergedSeeds", "MergedSeeds")}) {}

// ------------------------------------------------------------

StatusCode ClusterSeedMerging::initialize() {
  info() << "ClusterSeedMerging: MergeDistance = " << m_mergeDistance.value() << " mm" << endmsg;

  return StatusCode::SUCCESS;
} // initialize

// ------------------------------------------------------------

std::tuple<edm4hep::ClusterCollection>
ClusterSeedMerging::operator()(const std::vector<const edm4hep::ClusterCollection*>& caloSeeds,
                               const std::vector<const edm4hep::ClusterCollection*>& trackSeeds) const {

  edm4hep::ClusterCollection mergedOut;

  // ------------------------------------------------------------------
  // Step 1: Build a flat list of nodes: (type, original cluster index,
  //         theta, phi).
  // ------------------------------------------------------------------
  struct Node {
    int type;               // seed types: 1 = calo-driven A, 2 = calo-driven B, 4 = track-driven C
    int collIdx;            // index of the source collection in the input vector
    int srcIdx;             // index in the source collection
    float x, y, z;          // position of the seed (mm)
    uint64_t trackSeedCell; // Type-C only: the crystal the track extrapolates onto, else 0

    // Scratch for the Step 3 split: the opening angle to the track seed that owns this
    // node, and which one (index into that component's trackNodes).  Each node belongs
    // to exactly one component, so this is written at most once.
    float reachDist{std::numeric_limits<float>::max()};
    size_t reachOwner{0};
  };

  // One merged group of seeds.  The fields are filled in stages: nodes and
  // hasTrackSeed when the group is formed in Step 3, cellIDs and absorbedBy by the
  // absorption in Step 4, and types/hits/state by Step 5.
  struct Component {
    std::vector<int> nodes;                 // indices into nodes[]
    bool hasTrackSeed{false};               // at most one Type-C seed per group
    uint64_t anchorCell{0};                 // its track seed's crystal: never redistributed
    std::unordered_set<uint64_t> cellIDs{}; // union of its seeds' cell IDs
    int absorbedBy{-1};                     // >= 0: swallowed by that component
    int types{0};                           // OR of its seeds' Cluster::type
    ClusterSeedingBase::Hitmap hits{};      // deduplicated hits, refined by Step 5.5
    ClusterState state{};                   // frozen barycenter used for redistribution
  };

  std::vector<Node> nodes;
  size_t total = 0;
  for (const auto* c : caloSeeds)
    total += c->size();
  for (const auto* c : trackSeeds)
    total += c->size();
  nodes.reserve(total);

  auto addNodes = [&nodes](const edm4hep::ClusterCollection* coll, int collIdx) {
    for (int i = 0; i < static_cast<int>(coll->size()); ++i) {
      const auto& cl = (*coll)[i];
      const auto& p = cl.getPosition();
      const int type = cl.getType();
      // TrackDrivenClusterSeeding attaches the crystal the track points at ahead of the
      // neighbourhood, so for a Type-C seed it is the first hit.
      const uint64_t seedCell =
          (ClusterSeeding::hasSeed(type, ClusterSeeding::SeedType::TrackDrivenC) && cl.hits_size() > 0)
              ? (*cl.getHits().begin()).getCellID()
              : 0;
      nodes.push_back({type, collIdx, i, p.x, p.y, p.z, seedCell});
    }
  }; // lambda addNodes

  for (size_t i = 0; i < caloSeeds.size(); ++i)
    addNodes(caloSeeds[i], i);

  for (size_t i = 0; i < trackSeeds.size(); ++i)
    addNodes(trackSeeds[i], i + caloSeeds.size()); // track seeds come after calo seeds

  const int n = static_cast<int>(nodes.size());
  if (n == 0) {
    debug() << "ClusterSeedMerging: no input seeds - returning empty collections." << endmsg;
    return std::make_tuple(std::move(mergedOut));
  }

  // ------------------------------------------------------------------
  // Step 2: Adjacency predicate.
  //   Two seeds i and j are adjacent when the opening angle between their
  //   position vectors is less than arctan(MergeDistance / |pos_i|):
  //
  //     alpha(i,j) < arctan(d / min(r_i, r_j))
  //
  //   The C+C rule is enforced later in BFS.
  // ------------------------------------------------------------------
  const float mergeDist = m_mergeDistance.value();

  auto adjacent = [&nodes, mergeDist](int i, int j) -> bool {
    const Node& ni = nodes[i];
    const Node& nj = nodes[j];
    const float r2i = ni.x * ni.x + ni.y * ni.y + ni.z * ni.z;
    const float r2j = nj.x * nj.x + nj.y * nj.y + nj.z * nj.z;
    if (r2i <= 0.f || r2j <= 0.f)
      return false;

    const float dotij = ni.x * nj.x + ni.y * nj.y + ni.z * nj.z;
    if (dotij <= 0.f)
      return false; // opposite hemisphere - never adjacent

    // opening angle between the two position vectors
    const float cosAlpha = dotij / std::sqrt(r2i * r2j);
    const float alpha = std::acos(std::max(-1.f, std::min(1.f, cosAlpha)));

    // Union of the two cones: adjacent when either seed lies inside the other's.
    // The smaller radius gives the wider cone.
    return alpha < std::atan2(mergeDist, std::sqrt(std::min(r2i, r2j)));
  }; // lambda adjacent

  // ------------------------------------------------------------------
  // Step 3: BFS to find connected components, then enforce the invariant that
  //   each component holds AT MOST ONE Type-C node.  A component holding several
  //   is split into one component per Type-C node, and every other seed joins the
  //   closest in opening angle among the Type-C nodes that can reach it along the
  //   chain of adjacent seeds.
  // ------------------------------------------------------------------
  auto isTrackSeed = [&nodes](int idx) {
    return ClusterSeeding::hasSeed(nodes[idx].type, ClusterSeeding::SeedType::TrackDrivenC);
  };

  std::vector<Component> components;
  std::vector<bool> visited(n, false);

  for (int i = 0; i < n; ++i) {
    if (visited[i])
      continue;

    std::vector<int> comp;
    std::vector<int> trackNodes;
    std::queue<int> q;
    q.push(i);

    while (!q.empty()) {
      const int cur = q.front();
      q.pop();
      if (visited[cur])
        continue;

      visited[cur] = true;
      comp.push_back(cur);
      if (isTrackSeed(cur))
        trackNodes.push_back(cur);

      // add all adjacent nodes to the queue
      for (int j = 0; j < n; ++j) {
        if (!visited[j] && adjacent(cur, j))
          q.push(j);
      }
    } // while queue not empty

    if (trackNodes.size() <= 1) {
      const uint64_t anchor = trackNodes.empty() ? 0 : nodes[trackNodes.front()].trackSeedCell;
      components.push_back({std::move(comp), !trackNodes.empty(), anchor});
      continue;
    }

    // Several track seeds: one component each.  A track seed can claim only the seeds it
    // reaches through the adjacency chain, and among the candidates the nearest in angle wins.
    const size_t firstOfSplit = components.size();
    for (const int t : trackNodes)
      components.push_back({{t}, true, nodes[t].trackSeedCell});

    for (size_t k = 0; k < trackNodes.size(); ++k) {
      const Node& track = nodes[trackNodes[k]];
      std::unordered_set<int> reached;
      std::queue<int> frontier;
      frontier.push(trackNodes[k]);

      while (!frontier.empty()) {
        const int cur = frontier.front();
        frontier.pop();

        for (const int nb : comp) {
          if (nb == cur || isTrackSeed(nb) || !adjacent(cur, nb))
            continue; // not on the chain, or another track seed
          if (!reached.insert(nb).second)
            continue; // this track has already walked through it

          const float d = openingAngleDist({track.x, track.y, track.z, 0.f}, nodes[nb].x, nodes[nb].y, nodes[nb].z);
          if (d < nodes[nb].reachDist) {
            nodes[nb].reachDist = d;
            nodes[nb].reachOwner = k;
          }

          frontier.push(nb);
        } // loop over the other seeds of this component
      } // while the chain still grows
    } // loop over track seeds

    for (const int idx : comp) {
      if (!isTrackSeed(idx))
        components[firstOfSplit + nodes[idx].reachOwner].nodes.push_back(idx);
    }
  } // loop over nodes

  // ------------------------------------------------------------------
  // Step 4: Absorb completely-overlapping components.
  //   If every hit of component i is already present in component j (i ⊆ j),
  //   then i is absorbed into j: its node type bits are merged into j and it
  //   produces no independent output cluster.  We iterate until stable to
  //   handle transitive chains of subset relations.
  // ------------------------------------------------------------------
  const int nComp = static_cast<int>(components.size());

  // Helper: resolve source collection from collIdx
  auto retrieveSrcColl = [&](int collIdx) -> const edm4hep::ClusterCollection* {
    if (collIdx < static_cast<int>(caloSeeds.size()))
      return caloSeeds[collIdx];

    return trackSeeds[collIdx - static_cast<int>(caloSeeds.size())];
  };

  // Build per-component cellID sets for subset testing
  for (int ci = 0; ci < nComp; ++ci) {
    for (const int idx : components[ci].nodes) {
      const Node& nd = nodes[idx];
      for (const auto& hit : (*retrieveSrcColl(nd.collIdx))[nd.srcIdx].getHits())
        components[ci].cellIDs.insert(hit.getCellID());
    }
  }

  bool anyAbsorbed = true;
  while (anyAbsorbed) {
    anyAbsorbed = false;
    for (int i = 0; i < nComp; ++i) {
      if (components[i].absorbedBy >= 0 || components[i].cellIDs.empty())
        continue;

      for (int j = 0; j < nComp; ++j) {
        if (i == j || components[j].absorbedBy >= 0)
          continue;
        if (components[i].hasTrackSeed && components[j].hasTrackSeed)
          continue; // absorbing would put two track seeds in one group
        if (components[j].cellIDs.size() < components[i].cellIDs.size())
          continue; // j must be at least as large as i

        // Test i ⊆ j
        bool isSubset = true;
        for (const uint64_t cid : components[i].cellIDs) {
          if (!components[j].cellIDs.count(cid)) {
            isSubset = false;
            break;
          }
        }

        if (isSubset) {
          components[i].absorbedBy = j;
          // j inherits the track seed and its anchor; the guard above means j has none of its own
          if (components[i].hasTrackSeed) {
            components[j].hasTrackSeed = true;
            components[j].anchorCell = components[i].anchorCell;
          }
          // Merge i's nodes into j so the type bitmask is collected in Step 6
          for (const int idx : components[i].nodes)
            components[j].nodes.push_back(idx);

          anyAbsorbed = true;
          break;
        }
      } // loop over j
    } // loop over i
  } // while anyAbsorbed

  const int nSurviving = static_cast<int>(
      std::count_if(components.begin(), components.end(), [](const Component& c) { return c.absorbedBy < 0; }));
  debug() << "ClusterSeedMerging: " << nSurviving << " output components after absorbing " << (nComp - nSurviving)
          << " fully-overlapping components." << endmsg;

  // ------------------------------------------------------------------
  // Step 5: Partial-overlap hit redistribution via opening-angle distance.
  //   For every cellID that appears in two or more non-absorbed components,
  //   assign it exclusively to the component with the smallest opening-angle
  //   distance:  d = 1 - cos(angle)
  //   E_clus and the cluster direction come from the cluster state
  //   *before* redistribution (frozen snapshot).
  // ------------------------------------------------------------------

  // Build per-component deduplicated hit maps, frozen cluster states, and
  // the cellOwners index (reused in Step 7).
  std::unordered_map<uint64_t, std::vector<int>> cellOwners;

  for (int ci = 0; ci < nComp; ++ci) {
    if (components[ci].absorbedBy >= 0)
      continue;

    for (const int idx : components[ci].nodes) {
      const Node& nd = nodes[idx];
      components[ci].types |= nd.type;

      for (const auto& hit : (*retrieveSrcColl(nd.collIdx))[nd.srcIdx].getHits())
        components[ci].hits.try_emplace(hit.getCellID(), hit);
    } // loop over nodes in component to build hit map

    if (!components[ci].hits.empty()) {
      auto bary = calcBarycenter(components[ci].hits);
      components[ci].state = {bary.x, bary.y, bary.z, bary.energy};

      for (const auto& [cellID, hit] : components[ci].hits)
        cellOwners[cellID].push_back(ci);
    } // if component has any hits
  } // loop over components

  // Resolve each contested hit
  int nRedistributed = 0;
  for (auto& [cellID, owners] : cellOwners) {
    if (owners.size() < 2)
      continue; // not contested

    // Retrieve the hit from the first owner (all copies are identical)
    const edm4hep::CalorimeterHit& hit = components[owners[0]].hits.at(cellID);
    const edm4hep::Vector3f& hpos = hit.getPosition();

    int winner = owners[0];
    float bestDist = std::numeric_limits<float>::max();

    for (const int ci : owners) {
      if (components[ci].hasTrackSeed && components[ci].anchorCell == cellID) {
        winner = ci;
        break; // a track never loses the crystal it points at
      }

      const ClusterState& cs = components[ci].state;
      const float d = openingAngleDist(cs, hpos.x, hpos.y, hpos.z);

      if (d < bestDist) {
        bestDist = d;
        winner = ci;
      }
    } // loop over owners to find winner

    // Remove from all losers
    for (const int ci : owners) {
      if (ci != winner) {
        components[ci].hits.erase(cellID);
        ++nRedistributed;
      }
    } // loop over owners to remove losers
  } // loop over contested hits

  debug() << "ClusterSeedMerging: redistributed " << nRedistributed << " contested hits via opening-angle distance."
          << endmsg;

  // ------------------------------------------------------------------
  // Step 6: Build output collections.
  //   Non-absorbed components -> one merged cluster with barycenter position.
  //   Absorbed components     -> skipped (their hits+types are in the absorbing component).
  //   Uses the per-component hits built/refined by Step 5.5.
  // ------------------------------------------------------------------
  for (int ci = 0; ci < nComp; ++ci) {
    if (components[ci].absorbedBy >= 0)
      continue; // swallowed by another component

    auto& hitMap = components[ci].hits;
    if (hitMap.empty())
      continue; // no hits left after redistribution - skip

    auto merged = mergedOut.create();

    auto bary = calcBarycenter(hitMap);
    merged.setPosition({bary.x, bary.y, bary.z});
    merged.setEnergy(bary.energy);
    merged.setType(components[ci].types);

    for (const auto& [cellID, hit] : hitMap)
      merged.addToHits(hit);
  } // loop over components

  debug() << "ClusterSeedMerging: " << mergedOut.size() << " merged groups." << endmsg;

  return std::make_tuple(std::move(mergedOut));
}

DECLARE_COMPONENT(ClusterSeedMerging)
