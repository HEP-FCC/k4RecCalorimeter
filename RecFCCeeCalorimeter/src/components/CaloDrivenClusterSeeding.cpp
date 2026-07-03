#include "CaloDrivenClusterSeeding.h"

#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/MutableCluster.h"

#include "DD4hep/Detector.h"
#include "DDSegmentation/Segmentation.h"

#include "ClusterSeedMerging.h" // for enum

#include <cstdint>
#include <queue>
#include <unordered_map>
#include <vector>

// ============================================================
//  CaloDrivenClusterSeeding
// ============================================================

CaloDrivenClusterSeeding::CaloDrivenClusterSeeding(const std::string& name, ISvcLocator* svcLoc)
    : ClusterSeedingBase(name, svcLoc, {KeyValues("InputCaloHitCollections", {})},
                         {KeyValue("OutputSeedsA", "CaloDrivenSeedsA"), KeyValue("OutputSeedsB", "CaloDrivenSeedsB")}) {
}

StatusCode CaloDrivenClusterSeeding::initialize() { return ClusterSeedingBase::initialize(); }

// ------------------------------------------------------------

std::tuple<edm4hep::ClusterCollection, edm4hep::ClusterCollection>
CaloDrivenClusterSeeding::operator()(const std::vector<const edm4hep::CalorimeterHitCollection*>& caloHitColls) const {

  edm4hep::ClusterCollection seedsA;
  edm4hep::ClusterCollection seedsB;

  const float thrA = m_seedEnergyThresholdA.value();
  const float thrB = m_seedEnergyThresholdB.value();

  // ------------------------------------------------------------------
  // Step 1: Accumulate ECAL scintillation depth-0 energy per cellID.
  //         Also keep a representative CalorimeterHit per cellID so we
  //         can attach it to the output cluster.
  // ------------------------------------------------------------------
  std::unordered_map<uint64_t, float> energyMap;
  ClusterSeedingBase::Hitmap hitMap;

  // create lookup maps by looping over all input collections and hits
  for (const auto* coll : caloHitColls) {
    if (!coll)
      continue;

    for (const auto& hit : *coll) {
      const uint64_t cellID = hit.getCellID();

      if (!passSelection(cellID))
        continue;

      energyMap[cellID] = hit.getEnergy();
      hitMap.try_emplace(cellID, hit);
    }
  }

  // ------------------------------------------------------------------
  // Step 2: For each above-threshold crystal, check seeding conditions.
  // ------------------------------------------------------------------

  auto isLocalMax = [&energyMap, this](const uint64_t cellID, const float energy,
                                       const std::set<uint64_t>& nbrs) -> bool {
    for (const uint64_t nb : nbrs) {
      if (nb == cellID)
        continue;

      auto it = energyMap.find(nb);
      if (it != energyMap.end() && it->second > energy) {
        return false;
      }
    }

    return true;
  }; // lambda isLocalMax

  struct SeedCandidate {
    uint64_t cellID;
    std::set<uint64_t> neighbors;
  };

  std::vector<SeedCandidate> seedsTypeA, seedsTypeB; // store cellIDs of seeds

  for (const auto& [cellID, energy] : energyMap) {
    if (!passSelection(cellID))
      continue;

    // --- Type A ---
    // the seed crystal exceeds the threshold
    // and is a local maximum within its Von Neumann neighbourhood
    if (energy >= thrA) {
      const auto nbrs = vonNeumannNeighbors(cellID, m_vnDistance.value());

      if (isLocalMax(cellID, energy, nbrs))
        seedsTypeA.push_back({cellID, nbrs}); // store the whole neighborhood as a seed
    } // end Type A

    // --- Type B ---
    // the seed crystal exceeds the threshold
    // in the same VN neighbourhood, all above-threshold neighbours are connected
    // and the seed is a local maximum within the neighbourhood
    if (energy >= thrB) {
      const auto nbrs = vonNeumannNeighbors(cellID, m_vnDistance.value());

      // Collect above-threshold neighbours (including self)
      std::set<uint64_t> aboveThresh;
      for (const uint64_t nb : nbrs) {
        auto it = energyMap.find(nb);
        if (it != energyMap.end() && it->second >= thrB) {
          aboveThresh.insert(nb);
        }
      } // loop over neighbours

      // check above threshold neighbour count
      // and connectivity of above-threshold neighbours
      if (aboveThresh.size() >= m_minAboveThresholdNeighbours.value() && isConnected(aboveThresh)) {
        // Local maximum
        if (isLocalMax(cellID, energy, nbrs))
          seedsTypeB.push_back({cellID, nbrs}); // store the whole neighborhood as a seed
      } // end if neighbor count and connectivity
    } // end Type B
  } // loop over energyMap

  // ------------------------------------------------------------------
  // Step 3: Create output clusters for the seeds
  for (const auto& seed : seedsTypeA) {
    if (seed.neighbors.empty())
      continue;

    auto cluster = seedsA.create();
    cluster.setType(static_cast<int>(ClusterSeeding::SeedType::CaloDrivenA)); // Type A seed
    // use the seed hit position as the cluster position
    cluster.setPosition(hitMap[seed.cellID].getPosition());
    // attach only hits that have actual energy deposits
    for (const auto& cellID : seed.neighbors) {
      auto it = hitMap.find(cellID);
      if (it != hitMap.end())
        cluster.addToHits(it->second);
    }
  } // loop over seedsTypeA

  for (const auto& seed : seedsTypeB) {
    if (seed.neighbors.empty())
      continue;

    auto cluster = seedsB.create();
    cluster.setType(static_cast<int>(ClusterSeeding::SeedType::CaloDrivenB)); // Type B seed
    // use the seed hit position as the cluster position
    cluster.setPosition(hitMap[seed.cellID].getPosition());
    // attach only hits that have actual energy deposits
    for (const auto& cellID : seed.neighbors) {
      auto it = hitMap.find(cellID);
      if (it != hitMap.end())
        cluster.addToHits(it->second);
    }
  } // loop over seedsTypeB

  debug() << "CaloDrivenClusterSeeding: found " << seedsTypeA.size() << " Type-A seeds and " << seedsTypeB.size()
          << " Type-B seeds." << endmsg;

  return std::make_tuple(std::move(seedsA), std::move(seedsB));
} // operator()

DECLARE_COMPONENT(CaloDrivenClusterSeeding)
