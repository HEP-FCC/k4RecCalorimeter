#include "TrackDrivenClusterSeeding.h"

#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/MutableCluster.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackState.h"

#include "DD4hep/Detector.h"
#include "DDSegmentation/Segmentation.h"

#include <cmath>
#include <limits>
#include <set>
#include <unordered_map>
#include <vector>

// ============================================================
//  TrackDrivenClusterSeeding
// ============================================================

TrackDrivenClusterSeeding::TrackDrivenClusterSeeding(const std::string& name, ISvcLocator* svcLoc)
    : ClusterSeedingBase(
          name, svcLoc,
          {KeyValue("InputTrackCollection", "TracksFromGenParticles"), KeyValues("InputCaloHitCollections", {})},
          {KeyValue("OutputSeedsC", "TrackDrivenSeedsC")}) {}

// ------------------------------------------------------------

StatusCode TrackDrivenClusterSeeding::initialize() { return ClusterSeedingBase::initialize(); }

// ------------------------------------------------------------

std::tuple<edm4hep::ClusterCollection>
TrackDrivenClusterSeeding::operator()(const edm4hep::TrackCollection& trackColl,
                                      const std::vector<const edm4hep::CalorimeterHitCollection*>& caloHitColls) const {

  edm4hep::ClusterCollection seedsC;

  const float thr = m_seedEnergyThreshold.value();
  const float window = m_trackWindow.value();

  // ------------------------------------------------------------------
  // Step 1: Build energy and position maps for above-threshold hits
  //         that pass the cell-ID selection filter.
  // ------------------------------------------------------------------
  std::unordered_map<uint64_t, float> energyMap;
  ClusterSeedingBase::Hitmap hitMap;

  for (const auto* coll : caloHitColls) {
    if (!coll)
      continue;
    for (const auto& hit : *coll) {
      const uint64_t cellID = hit.getCellID();
      if (!passSelection(cellID))
        continue;
      energyMap[cellID] += hit.getEnergy();
      hitMap.try_emplace(cellID, hit);
    }
  }

  // ------------------------------------------------------------------
  // Step 2: Pre-compute (theta, phi) for each hit to avoid repeated
  //         sqrt/atan2 calls during the per-track search.
  // ------------------------------------------------------------------

  std::unordered_map<uint64_t, std::pair<float, float>> posMap; // cellID -> (theta, phi)
  posMap.reserve(hitMap.size());
  for (const auto& [cellID, hit] : hitMap) {
    const auto& p = hit.getPosition();
    posMap[cellID] = ClusterSeeding::toThetaPhi(p.x, p.y, p.z);
  }

  // ------------------------------------------------------------------
  // Step 3: Collect all track states at the calorimeter surface.
  //         edm4hep TrackState location == 4: AtCalorimeter.
  // ------------------------------------------------------------------
  std::vector<std::pair<float, float>> trackImpacts; // (theta, phi)
  for (const auto& track : trackColl) {
    for (const auto& ts : track.getTrackStates()) {
      if (ts.location != edm4hep::TrackState::AtCalorimeter)
        continue; // AtCalorimeter only

      const auto& rp = ts.referencePoint;
      trackImpacts.push_back(ClusterSeeding::toThetaPhi(rp.x, rp.y, rp.z));
      break;
    }
  } // loop over tracks

  // ------------------------------------------------------------------
  // Step 4: For each track state on calo, find the best seed candidate.
  // ------------------------------------------------------------------
  std::set<uint64_t> usedSeeds; // avoid duplicating a seed for two close tracks

  for (const auto& [tth, tph] : trackImpacts) {
    uint64_t bestCell = 0;
    std::set<uint64_t> bestNbrs; // best seed + VN-d1 neighbors that pass selection
    float bestDist = std::numeric_limits<float>::max();

    for (const auto& [cellID, energy] : energyMap) {
      if (energy < thr)
        continue;

      auto posIt = posMap.find(cellID);
      if (posIt == posMap.end())
        continue;

      const auto& [cth, cph] = posIt->second;

      // Angular distance from track impact
      const float dist = ClusterSeeding::angularDist(cth, cph, tth, tph);
      if (dist >= window)
        continue;

      // Local maximum check within VN neighbourhood
      std::set<uint64_t> rawNbrs;
      m_segmentation->neighbours(cellID, rawNbrs);
      // filter neighbors through passSelection
      std::set<uint64_t> filteredNbrs;
      for (const auto& nb : rawNbrs) {
        if (passSelection(nb))
          filteredNbrs.insert(nb);
      }

      bool isLocalMax = true;
      for (const uint64_t nb : filteredNbrs) {
        if (nb == cellID)
          continue;

        auto it = energyMap.find(nb);
        if (it != energyMap.end() && it->second > energy) {
          isLocalMax = false;
          break;
        }
      } // for filtered neighbors

      if (!isLocalMax)
        continue;

      if (dist < bestDist) {
        bestDist = dist;
        bestCell = cellID;
        bestNbrs = std::move(filteredNbrs);
      }
    } // loop over hits

    if (bestCell == 0)
      continue; // no candidate found for this track
    if (usedSeeds.count(bestCell))
      continue; // already seeded by another track
    usedSeeds.insert(bestCell);

    // Build the cluster: seed hit + above-threshold VN-d1 neighbourhood
    auto cluster = seedsC.create();
    cluster.setType(static_cast<int>(ClusterSeeding::SeedType::TrackDrivenC)); // Type C seed (track-driven)
    // use the seed hit position as the cluster position
    cluster.setPosition(hitMap[bestCell].getPosition());

    // Add the seed cell itself. m_segmentation->neighbours() returns only the
    // surrounding cells, so bestCell (the local-max that defines the seed) is NOT in
    // bestNbrs and would otherwise be dropped -- leaving isolated MIP seeds with 0 hits.
    cluster.addToHits(hitMap[bestCell]);

    for (const auto& id : bestNbrs) {
      if (id == bestCell)
        continue; // guard against double-add if neighbours() ever includes the query cell
      auto it = hitMap.find(id);
      if (it != hitMap.end())
        cluster.addToHits(it->second);
    }
  } // loop over track impacts

  debug() << "TrackDrivenClusterSeeding: found " << seedsC.size() << " Type-C seeds from " << trackImpacts.size()
          << " track states." << endmsg;

  return std::make_tuple(std::move(seedsC));
}

DECLARE_COMPONENT(TrackDrivenClusterSeeding)
