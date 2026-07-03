#include "ClusterSeedMerging.h"

#include "edm4hep/Cluster.h"
#include "edm4hep/MutableCluster.h"

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
    int type;      // seed types: 1 = calo-driven A, 2 = calo-driven B, 4 = track-driven C
    int collIdx;   // index of the source collection in the input vector
    int srcIdx;    // index in the source collection
    float x, y, z; // position of the seed (mm)
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
      nodes.push_back({type, collIdx, i, p.x, p.y, p.z});
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
  //     alpha(i,j) < arctan(d / r_i)
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

    return alpha < std::atan2(mergeDist, std::sqrt(r2i));
  }; // lambda adjacent

  // ------------------------------------------------------------------
  // Step 3: BFS to find connected components.
  //   Invariant: each component holds AT MOST ONE Type-C node.
  //   If a second C node is encountered while growing a component, it is
  //   skipped (left unvisited) so it will seed its own component later.
  // ------------------------------------------------------------------
  std::vector<std::vector<int>> components;
  std::vector<bool> visited(n, false);

  for (int i = 0; i < n; ++i) {
    if (visited[i])
      continue;

    std::vector<int> comp;
    bool compHasC = false;
    std::queue<int> q;
    q.push(i);

    while (!q.empty()) {
      const int cur = q.front();
      q.pop();
      if (visited[cur])
        continue;

      // Enforce at-most-one-C rule
      if (nodes[cur].type == static_cast<int>(ClusterSeeding::SeedType::TrackDrivenC)) {
        if (compHasC)
          continue; // defer this C node to its own component

        compHasC = true;
      }

      visited[cur] = true;
      comp.push_back(cur);

      // add all adjacent nodes to the queue
      for (int j = 0; j < n; ++j) {
        if (!visited[j] && adjacent(cur, j))
          q.push(j);
      }
    } // while queue not empty

    components.push_back(std::move(comp));
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
  std::vector<std::unordered_set<uint64_t>> compCellIDs(nComp);
  for (int ci = 0; ci < nComp; ++ci) {
    for (const int idx : components[ci]) {
      const Node& nd = nodes[idx];
      for (const auto& hit : (*retrieveSrcColl(nd.collIdx))[nd.srcIdx].getHits())
        compCellIDs[ci].insert(hit.getCellID());
    }
  }

  // absorbed[i] = j means component i is swallowed by component j
  std::vector<int> absorbed(nComp, -1);

  bool anyAbsorbed = true;
  while (anyAbsorbed) {
    anyAbsorbed = false;
    for (int i = 0; i < nComp; ++i) {
      if (absorbed[i] >= 0 || compCellIDs[i].empty())
        continue;

      for (int j = 0; j < nComp; ++j) {
        if (i == j || absorbed[j] >= 0)
          continue;
        if (compCellIDs[j].size() < compCellIDs[i].size())
          continue; // j must be at least as large as i

        // Test i ⊆ j
        bool isSubset = true;
        for (const uint64_t cid : compCellIDs[i]) {
          if (!compCellIDs[j].count(cid)) {
            isSubset = false;
            break;
          }
        }

        if (isSubset) {
          absorbed[i] = j;
          // Merge i's nodes into j so the type bitmask is collected in Step 6
          for (const int idx : components[i])
            components[j].push_back(idx);

          anyAbsorbed = true;
          break;
        }
      } // loop over j
    } // loop over i
  } // while anyAbsorbed

  debug() << "ClusterSeedMerging: " << std::count(absorbed.begin(), absorbed.end(), -1)
          << " output components after absorbing " << (nComp - std::count(absorbed.begin(), absorbed.end(), -1))
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
  std::vector<ClusterState> frozenState(nComp, {0.f, 0.f, 0.f, 0.f});
  std::vector<ClusterSeedingBase::Hitmap> compHitMaps(nComp);
  std::vector<int> compTypes(nComp, 0);
  std::unordered_map<uint64_t, std::vector<int>> cellOwners;

  for (int ci = 0; ci < nComp; ++ci) {
    if (absorbed[ci] >= 0)
      continue;

    for (const int idx : components[ci]) {
      const Node& nd = nodes[idx];
      compTypes[ci] |= nd.type;

      for (const auto& hit : (*retrieveSrcColl(nd.collIdx))[nd.srcIdx].getHits())
        compHitMaps[ci].try_emplace(hit.getCellID(), hit);
    } // loop over nodes in component to build hit map

    if (!compHitMaps[ci].empty()) {
      auto bary = calcBarycenter(compHitMaps[ci]);
      frozenState[ci] = {bary.x, bary.y, bary.z, bary.energy};

      for (const auto& [cellID, hit] : compHitMaps[ci])
        cellOwners[cellID].push_back(ci);
    } // if component has any hits
  } // loop over components

  // Resolve each contested hit
  int nRedistributed = 0;
  for (auto& [cellID, owners] : cellOwners) {
    if (owners.size() < 2)
      continue; // not contested

    // Retrieve the hit from the first owner (all copies are identical)
    const edm4hep::CalorimeterHit& hit = compHitMaps[owners[0]].at(cellID);
    const edm4hep::Vector3f& hpos = hit.getPosition();

    int winner = owners[0];
    float bestDist = std::numeric_limits<float>::max();

    for (const int ci : owners) {
      const ClusterState& cs = frozenState[ci];
      const float d = openingAngleDist(cs, hpos.x, hpos.y, hpos.z);

      if (d < bestDist) {
        bestDist = d;
        winner = ci;
      }
    } // loop over owners to find winner

    // Remove from all losers
    for (const int ci : owners) {
      if (ci != winner) {
        compHitMaps[ci].erase(cellID);
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
  //   Uses compHitMaps built/refined by Step 5.5.
  // ------------------------------------------------------------------
  for (int ci = 0; ci < nComp; ++ci) {
    if (absorbed[ci] >= 0)
      continue; // swallowed by another component

    auto& hitMap = compHitMaps[ci];
    if (hitMap.empty())
      continue; // no hits left after redistribution - skip

    auto merged = mergedOut.create();

    auto bary = calcBarycenter(hitMap);
    merged.setPosition({bary.x, bary.y, bary.z});
    merged.setEnergy(bary.energy);
    merged.setType(compTypes[ci]);

    for (const auto& [cellID, hit] : hitMap)
      merged.addToHits(hit);
  } // loop over components

  debug() << "ClusterSeedMerging: " << mergedOut.size() << " merged groups." << endmsg;

  return std::make_tuple(std::move(mergedOut));
}

DECLARE_COMPONENT(ClusterSeedMerging)
