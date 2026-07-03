#include "ClusterSeedGrower.h"
#include "ClusterSeedMerging.h" // for enum

#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/MutableCluster.h"

#include "DD4hep/Detector.h"
#include "DDSegmentation/Segmentation.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <queue>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// ============================================================
//  ClusterSeedGrower
// ============================================================

ClusterSeedGrower::ClusterSeedGrower(const std::string& name, ISvcLocator* svcLoc)
    : ClusterSeedingBase(name, svcLoc, {KeyValues("InputSeeds", {}), KeyValues("InputHits", {})},
                         {KeyValue("OutputClusters", "TopoGrownClusters")}) {}

// ------------------------------------------------------------

StatusCode ClusterSeedGrower::initialize() {
  if (ClusterSeedingBase::initialize().isFailure())
    return StatusCode::FAILURE;

  if (m_fieldStringsToInclude.size() != m_fieldValuesToInclude.size()) {
    error() << "ClusterSeedGrower: FieldStringsToInclude and FieldValuesToInclude must have the same number of entries "
            << "(got " << m_fieldStringsToInclude.size() << " vs " << m_fieldValuesToInclude.size() << ")." << endmsg;
    return StatusCode::FAILURE;
  }

  try {
    for (const auto& field : m_fieldStringsToInclude.value())
      m_decoder->index(field); // check field existence
  } catch (const std::exception& ex) {
    error() << "ClusterSeedGrower: Error in cell ID filtering configuration: " << ex.what() << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_hardThreshold.value() > m_growThreshold.value()) {
    warning() << "ClusterSeedGrower: HardThreshold (" << m_hardThreshold.value() << " GeV) > GrowThreshold ("
              << m_growThreshold.value() << " GeV). No isolated hits will ever be assigned in the post-BFS phase."
              << endmsg;
  }

  return StatusCode::SUCCESS;
}

// ------------------------------------------------------------

auto ClusterSeedGrower::connectedSeeds(const Hitmap& pool, unsigned int vnDist, float threshold,
                                       unsigned int minHits) const -> std::vector<Hitmap> {
  // visited tracks which cells in *pool* have already been consumed
  // into a component.  Only cells present in *pool* are ever visited.
  std::unordered_map<uint64_t, bool> visited;
  visited.reserve(pool.size());
  for (const auto& [cid, hit] : pool)
    visited[cid] = false;

  std::vector<ClusterSeedingBase::Hitmap> result;

  for (const auto& [startCID, startHit] : pool) {
    if (visited[startCID])
      continue; // already consumed into a previous component
    if (startHit.getEnergy() < threshold)
      continue; // below seeding threshold

    // BFS: grow one connected component from startCID.
    // Only cells inside *pool* with E >= threshold are traversed.
    ClusterSeedingBase::Hitmap component;
    std::queue<uint64_t> bfsQueue;
    bfsQueue.push(startCID);

    while (!bfsQueue.empty()) {
      const uint64_t cur = bfsQueue.front();
      bfsQueue.pop();

      if (visited[cur])
        continue;
      visited[cur] = true;

      const auto& curHit = pool.at(cur);
      if (curHit.getEnergy() < threshold)
        continue; // below threshold - do not include

      component.emplace(cur, curHit);

      // Expand via VN distance vnDist.  vonNeumannNeighbors already
      // applies passSelection, so only geometrically valid neighbours
      // are returned.
      const std::set<uint64_t> vnn = vonNeumannNeighbors(cur, vnDist);
      for (const uint64_t nb : vnn) {
        if (nb == cur)
          continue;
        const auto vit = visited.find(nb);
        if (vit == visited.end() || vit->second)
          continue; // not in pool, or already visited
        if (pool.at(nb).getEnergy() >= threshold)
          bfsQueue.push(nb);
      }
    } // BFS

    if (component.size() >= minHits)
      result.push_back(std::move(component));
  } // loop over pool cells

  return result;
} // connectedSeeds

// ------------------------------------------------------------

auto ClusterSeedGrower::grow(ContestStrategy strategy, const Hitmap& pool, const std::vector<Hitmap>& seeds,
                             unsigned int vnDist, float growThreshold) const -> std::vector<Hitmap> {
  const int n = static_cast<int>(seeds.size());
  if (n == 0)
    return {};

  // Working copy: starts identical to seeds and accumulates new hits layer by layer.
  // Returned to the caller when done (possibly with fewer entries for MergeClusters).
  std::vector<ClusterSeedingBase::Hitmap> clusters(seeds);

  // --- O(1) "is this cell already assigned?" bookkeeping ---
  std::unordered_set<uint64_t> assigned;
  assigned.reserve(pool.size());
  for (const auto& hm : clusters)
    for (const auto& [cid, h] : hm)
      assigned.insert(cid);

  // --- Per-cluster BFS frontier: the cells added in the last layer. ---
  // In each new layer we expand outward from these cells.
  std::vector<std::set<uint64_t>> frontier(n);
  for (int i = 0; i < n; ++i)
    for (const auto& [cid, h] : clusters[i])
      frontier[i].insert(cid);

  // --- ResolveByDist: barycenter state per cluster, kept up-to-date each layer ---
  std::vector<ClusterState> states;
  if (strategy == ContestStrategy::ResolveByDist) {
    states.resize(n);
    for (int i = 0; i < n; ++i)
      states[i] = calcBarycenter(clusters[i]);
  }

  // --- Union-find for MergeClusters strategy ---
  // rep[i] is the "representative" of cluster i's group.
  // Initially every cluster is its own representative: rep[i] = i.
  //
  //   findRep(i): walk up the rep[] chain to the root (representative) of i's
  //               group. Uses path-halving (rep[x] = rep[rep[x]] each step) to
  //               keep chains short without extra storage.
  //
  //   mergeReps(a,b): make the two groups one by pointing b's representative
  //                   to a's representative.
  //
  // Two cluster indices that share the same representative after some merges
  // belong to the same logical cluster and will be folded together at the end.
  std::vector<int> rep(n);
  std::iota(rep.begin(), rep.end(), 0);
  auto findRep = [&](int x) {
    while (rep[x] != x) {
      rep[x] = rep[rep[x]];
      x = rep[x];
    }

    return x;
  };
  auto mergeReps = [&](int a, int b) {
    a = findRep(a);
    b = findRep(b);
    if (a != b)
      rep[b] = a;
  };

  // --- Helper shared by both strategy branches ---
  // Returns all cells in pool that neighbour fid (within vnDist), are not yet
  // assigned to any cluster, and are above growThreshold.
  auto eligibleNeighbors = [&](uint64_t fid) {
    std::vector<uint64_t> result;
    for (const uint64_t nb : vonNeumannNeighbors(fid, vnDist)) {
      if (nb == fid || assigned.count(nb))
        continue;

      const auto it = pool.find(nb);
      if (it == pool.end() || it->second.getEnergy() < growThreshold)
        continue;

      result.push_back(nb);
    }

    return result;
  };

  // --- Layer-by-layer BFS ---
  bool anyAdded = true;
  while (anyAdded) {
    anyAdded = false;

    if (strategy == ContestStrategy::ResolveByDist) {

      // Collect candidates: cellID -> list of claiming cluster indices.
      // A cluster may appear multiple times (via different frontier cells),
      // so we deduplicate before resolving.
      std::unordered_map<uint64_t, std::vector<int>> candidates;
      for (int i = 0; i < n; ++i)
        for (const uint64_t fid : frontier[i])
          for (const uint64_t nb : eligibleNeighbors(fid))
            candidates[nb].push_back(i);

      for (auto& [cid, claimers] : candidates) {
        std::sort(claimers.begin(), claimers.end());
        claimers.erase(std::unique(claimers.begin(), claimers.end()), claimers.end());
      }

      // Assign: winner is the cluster with the smallest distance.
      // All other claimers get nothing this layer.
      std::vector<std::set<uint64_t>> newFrontier(n);
      std::vector<bool> gained(n, false);
      for (const auto& [cid, claimers] : candidates) {
        const auto& hit = pool.at(cid);
        const auto& pos = hit.getPosition();
        int winner = claimers[0];
        if (claimers.size() > 1) {
          float best = std::numeric_limits<float>::max();
          for (const int i : claimers) {
            const float d = openingAngleDist(states[i], pos.x, pos.y, pos.z);
            if (d < best) {
              best = d;
              winner = i;
            }
          } // loop over claimers to find winner
        } // if multiple claimers

        assigned.insert(cid);
        clusters[winner].emplace(cid, hit);
        newFrontier[winner].insert(cid);
        gained[winner] = true;
        anyAdded = true;
      } // loop over candidates

      frontier = std::move(newFrontier);
      // Update barycenters for clusters that received new hits this layer.
      for (int i = 0; i < n; ++i) {
        if (gained[i])
          states[i] = calcBarycenter(clusters[i]);
      } // loop over clusters to update states
    } else { // MergeClusters
      // Several cluster indices may now share the same representative (after
      // earlier merges). Collapse their individual frontiers into one combined
      // frontier per representative so we grow from the full merged boundary.
      std::unordered_map<int, std::set<uint64_t>> repFrontier;
      for (int i = 0; i < n; ++i) {
        const int r = findRep(i);
        for (const uint64_t cid : frontier[i])
          repFrontier[r].insert(cid);
      }

      // Collect candidates: cellID -> set of representative indices that can reach it.
      std::unordered_map<uint64_t, std::set<int>> candidates;
      for (const auto& [r, fr] : repFrontier)
        for (const uint64_t fid : fr)
          for (const uint64_t nb : eligibleNeighbors(fid))
            candidates[nb].insert(r);

      // Assign the hit to the (merged) representative.
      // If multiple representatives can reach the same cell, merge them all
      // into one group — the cell then belongs to that merged cluster.
      std::unordered_map<int, std::set<uint64_t>> newFrontier;
      for (const auto& [cid, reps] : candidates) {
        const auto& hit = pool.at(cid);
        auto it = reps.begin();
        int winner = findRep(*it++);
        // If multiple representatives claim this cell, merge their clusters.
        for (; it != reps.end(); ++it)
          mergeReps(winner, findRep(*it));

        const int r = findRep(winner);
        assigned.insert(cid);
        clusters[r].emplace(cid, hit);
        newFrontier[r].insert(cid);
        anyAdded = true;
      } // loop over candidates

      // Propagate the new frontier back to per-index structures for next layer.
      for (int i = 0; i < n; ++i) {
        const int r = findRep(i);
        const auto it = newFrontier.find(r);
        frontier[i] = (it != newFrontier.end()) ? it->second : std::set<uint64_t>{};
      } // loop over clusters to update frontiers
    } // if MergeClusters strategy
  } // while anyAdded

  // For MergeClusters: fold every non-representative cluster into its
  // representative, then erase the now-empty entries.  The caller receives
  // exactly one Hitmap per surviving (representative) cluster.
  if (strategy == ContestStrategy::MergeClusters) {
    for (int i = 0; i < n; ++i) {
      const int r = findRep(i);
      if (r == i)
        continue;

      // representative r gets all of i's hits
      // i becomes empty and will be erased.
      for (auto& [cid, h] : clusters[i])
        clusters[r].emplace(cid, h);

      clusters[i].clear();
    }

    // Erase empty absorbed clusters.
    clusters.erase(std::remove_if(clusters.begin(), clusters.end(),
                                  [](const ClusterSeedingBase::Hitmap& hm) { return hm.empty(); }),
                   clusters.end());
  } // if MergeClusters

  return clusters;
} // grow

// ------------------------------------------------------------

void ClusterSeedGrower::attachIsolatedHits(const ClusterSeedingBase::Hitmap& pool,
                                           std::vector<ClusterSeedingBase::Hitmap>& seededClusters,
                                           std::vector<ClusterSeedingBase::Hitmap>& unseededClusters) const {
  // --- Snapshot barycenters for all clusters (seeded + unseeded) ---
  const int nSeeded = static_cast<int>(seededClusters.size());
  const int nUnseeded = static_cast<int>(unseededClusters.size());
  const int nTotal = nSeeded + nUnseeded;

  // for i < nSeeded   -> seededClusters[i]
  // for i >= nSeeded  -> unseededClusters[i - nSeeded]
  std::vector<ClusterState> statesBeforeAttach;
  statesBeforeAttach.reserve(nTotal);
  for (const auto& hm : seededClusters)
    statesBeforeAttach.push_back(calcBarycenter(hm));
  for (const auto& hm : unseededClusters)
    statesBeforeAttach.push_back(calcBarycenter(hm));

  const float cosGate = std::cos(m_maxOpeningAngle.value());

  // Build the set of already-owned hitPool cells.
  std::unordered_set<uint64_t> ownedCells;
  ownedCells.reserve(pool.size());
  for (const auto& hm : seededClusters)
    for (const auto& [cid, h] : hm)
      ownedCells.insert(cid);
  for (const auto& hm : unseededClusters)
    for (const auto& [cid, h] : hm)
      ownedCells.insert(cid);

  int nAttached = 0;
  for (const auto& [cid, hit] : pool) {
    if (ownedCells.count(cid))
      continue; // already in a cluster

    const auto& pos = hit.getPosition();
    const float rh = std::sqrt(pos.x * pos.x + pos.y * pos.y + pos.z * pos.z);
    if (rh < std::numeric_limits<float>::epsilon())
      continue; // degenerate hit position

    int bestIdx = -1;
    float bestDist = std::numeric_limits<float>::max();

    for (int i = 0; i < nTotal; ++i) {
      const ClusterState& cs = statesBeforeAttach[i];
      const float rc = std::sqrt(cs.x * cs.x + cs.y * cs.y + cs.z * cs.z);
      if (rc < std::numeric_limits<float>::epsilon())
        continue;

      // Opening-angle preselection: cosine of angle between cluster
      // barycenter direction and hit direction.
      const float cosAngle = (cs.x * pos.x + cs.y * pos.y + cs.z * pos.z) / (rc * rh);
      if (cosAngle < cosGate)
        continue; // opening angle > MaxOpeningAngle

      const float dist = openingAngleDist(cs, pos.x, pos.y, pos.z);
      if (dist < bestDist) {
        bestDist = dist;
        bestIdx = i;
      }
    } // loop over all clusters

    if (bestIdx < 0)
      continue; // no cluster within angle gate

    // Attach hit to the winning cluster.
    if (bestIdx < nSeeded)
      seededClusters[bestIdx].emplace(cid, hit);
    else
      unseededClusters[bestIdx - nSeeded].emplace(cid, hit);

    ++nAttached;
  } // loop over hitPool cells

  debug() << "ClusterSeedGrower: Step 5 attached " << nAttached << " previously unowned hit(s) to existing clusters "
          << "(MaxOpeningAngle = " << m_maxOpeningAngle.value() << " rad)." << endmsg;
} // attachIsolatedHits

// ------------------------------------------------------------

std::tuple<edm4hep::ClusterCollection>
ClusterSeedGrower::operator()(const std::vector<const edm4hep::ClusterCollection*>& seedColls,
                              const std::vector<const edm4hep::CalorimeterHitCollection*>& hitColls) const {

  edm4hep::ClusterCollection output;

  // ------------------------------------------------------------------
  // Step 1: Build unified hit pool: cellID -> CalorimeterHit.
  //   Contains only hits that pass field selection AND E >= HardThreshold.
  //   (Hits with HardThreshold <= E < GrowThreshold block BFS propagation
  //    but are eligible for isolated-hit assignment in Step 4.)
  // ------------------------------------------------------------------
  ClusterSeedingBase::Hitmap hitPool;
  ClusterSeedingBase::Hitmap leftoverHits;
  for (const auto* coll : hitColls) {
    if (!coll)
      continue;

    for (const auto& hit : *coll) {
      if (passSelection(hit.getCellID()) && hit.getEnergy() >= m_hardThreshold.value())
        hitPool.emplace(hit.getCellID(), hit);
      else if (hit.getEnergy() >= m_hardThreshold.value())
        leftoverHits.emplace(hit.getCellID(), hit);
    }
  } // loop over input hit collections

  // ------------------------------------------------------------------
  // Step 2: Count total clusters; initialise per-cluster data.
  // ------------------------------------------------------------------
  int nClusters = 0;
  for (const auto* coll : seedColls) {
    if (coll)
      nClusters += static_cast<int>(coll->size());
  }

  if (nClusters == 0) {
    warning() << "ClusterSeedGrower: no input seeds — returning empty collection." << endmsg;
    return std::make_tuple(std::move(output));
  }

  // Per-cluster hit maps (cellID → hit); initial frontier is the seed hits.
  // seedTypes records each input seed's type so Step 6 can reproduce it.
  // A temporary flat set handles shared seed hits: first seed wins.
  std::vector<ClusterSeedingBase::Hitmap> initClusterHits(nClusters);
  std::vector<int32_t> seedTypes(nClusters);
  {
    std::unordered_set<uint64_t> initAssigned;
    initAssigned.reserve(hitPool.size());
    int ci = 0;
    for (const auto* coll : seedColls) {
      if (!coll)
        continue;

      for (const auto& seed : *coll) {
        seedTypes[ci] = seed.getType();
        for (const auto& hit : seed.getHits()) {
          const uint64_t cid = hit.getCellID();
          if (!initAssigned.insert(cid).second)
            continue; // shared hit: first seed wins
                      // should have been resolved by ClusterSeedMerging upstream (if used)
          initClusterHits[ci].emplace(cid, hit);
        }
        ++ci;
      } // loop over seeds in this collection
    } // loop over seed collections
  } // wrapper scope for initAssigned and ci

  // ------------------------------------------------------------------
  // Step 3: BFS growth of seeded clusters (opening-angle distance contest resolution).
  // ------------------------------------------------------------------
  auto seededClusters =
      grow(ContestStrategy::ResolveByDist, hitPool, initClusterHits, m_vnDistSeeded.value(), m_growThreshold.value());

  // ------------------------------------------------------------------
  // Step 4: Unseeded cluster formation from unassigned hits.
  //
  //   a) Collect all hitPool cells not yet assigned to any seeded cluster.
  //   b) connectedSeeds finds groups of >= MinUnseededHits connected cells
  //      all above UnseededThreshold within VN distance VNDistUnseeded.
  //   c) grow() grows them (E >= GrowThreshold,
  //      VN distance VNDistUnseeded); a contested hit fully merges the
  //      two competing clusters (union-find).
  // ------------------------------------------------------------------
  const unsigned int vnUnseeded = m_vnDistUnseeded.value();
  const float unseededThr = m_unseededThreshold.value();
  const unsigned int minUHits = m_minUnseededHits.value();

  // Build unowned pool: hitPool cells not assigned to any seeded cluster.
  std::unordered_set<uint64_t> seededOwned;
  seededOwned.reserve(hitPool.size());
  for (const auto& hm : seededClusters)
    for (const auto& [cid, h] : hm)
      seededOwned.insert(cid);

  ClusterSeedingBase::Hitmap unownedPool;
  for (const auto& [cid, hit] : hitPool) {
    if (!seededOwned.count(cid))
      unownedPool.emplace(cid, hit);
  }

  const int nUnownedAfterSeeded = static_cast<int>(unownedPool.size());

  // Find connected seed groups, then grow them (merging on contest).
  auto unseededClusters =
      grow(ContestStrategy::MergeClusters, unownedPool, connectedSeeds(unownedPool, vnUnseeded, unseededThr, minUHits),
           vnUnseeded, m_growThreshold.value());

  // Count discarded hits (in unownedPool but not in any unseeded cluster)
  int nAssignedUnseeded = 0;
  for (const auto& hm : unseededClusters)
    nAssignedUnseeded += static_cast<int>(hm.size());

  const int nDiscarded = nUnownedAfterSeeded - nAssignedUnseeded;
  debug() << "ClusterSeedGrower: " << nUnownedAfterSeeded << " unowned hits after seeded BFS; "
          << unseededClusters.size() << " unseeded cluster(s) formed; " << nDiscarded
          << " hit(s) discarded after unseeded BFS." << endmsg;

  // ------------------------------------------------------------------
  // Step 5 (optional): Attach remaining unowned hitPool cells to the nearest cluster.
  //
  //   For every hitPool cell that belongs to neither seededClusters nor
  //   unseededClusters:
  //     1. Preselect clusters whose barycenter opening angle to the hit is
  //        < MaxOpeningAngle (snapshot barycenters computed once here; NOT
  //        updated between individual hit attachments).
  //     2. Among the preselected clusters, the one with the smallest
  //        opening-angle distance wins and the hit is appended to it.
  //   Hits for which no cluster passes the angle gate are left unassigned.
  // ------------------------------------------------------------------
  if (m_attachIsolatedHits.value())
    attachIsolatedHits(hitPool, seededClusters, unseededClusters);

  // ------------------------------------------------------------------
  // Step 6: Add companion hits (FieldStringsToInclude) to ALL clusters
  //   (seeded and unseeded).  e.g. Cherenkov hits that correspond to the
  //   same crystal but were excluded by the field selection filter.
  // ------------------------------------------------------------------
  auto addCompanions = [&](ClusterSeedingBase::Hitmap& hm) {
    for (unsigned int idx = 0; idx < m_fieldStringsToInclude.size(); ++idx) {
      const std::string& field = m_fieldStringsToInclude.value()[idx];
      const int64_t value = m_fieldValuesToInclude.value()[idx];
      // Snapshot before the loop so we don't re-process companions just added.
      const ClusterSeedingBase::Hitmap snapshot = hm;
      for (const auto& [origCid, hit] : snapshot) {
        uint64_t companionCid = origCid;
        m_decoder->set(companionCid, field, value); // derive the companion cellID
        const auto it = leftoverHits.find(companionCid);
        if (it != leftoverHits.end()) {
          hm.emplace(companionCid, it->second);
          leftoverHits.erase(it); // avoid duplicates — first come first served
        }
      } // loop over (original) hits
    } // loop over fields to include
  }; // lambda addCompanions

  // ------------------------------------------------------------------
  // Step 7: Build output ClusterCollection.
  //   Seeded clusters first (one per input seed, preserving seed type),
  //   then unseeded clusters.
  // ------------------------------------------------------------------
  for (int ci = 0; ci < static_cast<int>(seededClusters.size()); ++ci) {
    auto& hm = seededClusters[ci];
    addCompanions(hm);

    auto out = output.create();
    out.setType(seedTypes[ci]);

    const ClusterState finalState = calcBarycenter(hm);
    out.setPosition({finalState.x, finalState.y, finalState.z});
    out.setEnergy(finalState.energy);

    for (const auto& [cid, hit] : hm)
      out.addToHits(hit);
  } // loop over seeded clusters

  int nUnseededOut = 0;
  for (auto& hm : unseededClusters) {
    addCompanions(hm);

    auto out = output.create();
    out.setType(static_cast<int32_t>(ClusterSeeding::SeedType::Unseeded)); // unseeded cluster

    const ClusterState finalState = calcBarycenter(hm);
    out.setPosition({finalState.x, finalState.y, finalState.z});
    out.setEnergy(finalState.energy);

    for (const auto& [cid, hit] : hm)
      out.addToHits(hit);

    ++nUnseededOut;
  } // loop over unseeded clusters

  debug() << "ClusterSeedGrower: produced " << output.size() << " clusters total (" << (output.size() - nUnseededOut)
          << " seeded + " << nUnseededOut << " unseeded) from a hit pool of " << hitPool.size() << " cells." << endmsg;

  return std::make_tuple(std::move(output));
}

DECLARE_COMPONENT(ClusterSeedGrower)
