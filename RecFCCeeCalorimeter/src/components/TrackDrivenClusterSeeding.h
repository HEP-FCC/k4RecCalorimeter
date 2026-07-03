#ifndef TrackDrivenClusterSeeding_h
#define TrackDrivenClusterSeeding_h 1

/*
 * TrackDrivenClusterSeeding
 *
 * Gaudi MultiTransformer that identifies cluster seeds driven by charged-
 * particle track states at the calorimeter surface (Type-C seeds).
 *
 * Algorithm
 * ---------
 *   For each track state at the calorimeter surface (TrackState::AtCalorimeter,
 *   location == 4):
 *     1. Find all calorimeter hits (filtered by FieldStrings/FieldValues) whose
 *        angular distance from the track impact point is less than TrackWindow.
 *     2. Among those, keep only hits with E > SeedEnergyThreshold that are a
 *        local maximum within their Von Neumann distance-1 neighbourhood.
 *     3. Retain only the single closest hit to the track impact point.
 *
 *   One cluster (containing the single seed hit and its VN-d1 neighbourhood
 *   hits) is produced per track state, at most.  Two track states that resolve
 *   to the same seed crystal produce only one cluster.
 *
 * Input
 * -----
 *   edm4hep::TrackCollection
 *       Reconstructed charged-particle tracks.  All TrackState entries with
 *       location == AtCalorimeter (4) are used.
 *
 *   std::vector<const edm4hep::CalorimeterHitCollection*>
 *       Digitised calorimeter hit collections.  Filtered internally via
 *       FieldStrings / FieldValues (same mechanism as CaloDrivenClusterSeeding).
 *
 * Output
 * ------
 *   edm4hep::ClusterCollection
 *       One cluster per accepted track seed.  Each cluster holds the seed
 *       crystal plus all above-threshold hits from its VN-d1 neighbourhood.
 */

#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/TrackCollection.h"

#include "Gaudi/Property.h"
#include "k4FWCore/Transformer.h"

#include "ClusterSeedMerging.h" // for SeedType enum
#include "ClusterSeedingBase.h"

#include <cmath>
#include <limits>
#include <set>
#include <tuple>
#include <vector>

struct TrackDrivenClusterSeeding final
    : ClusterSeedingBase<std::tuple<edm4hep::ClusterCollection>(
          const edm4hep::TrackCollection&, const std::vector<const edm4hep::CalorimeterHitCollection*>&)> {

  TrackDrivenClusterSeeding(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode finalize() override { return StatusCode::SUCCESS; }

  std::tuple<edm4hep::ClusterCollection>
  operator()(const edm4hep::TrackCollection& trackColl,
             const std::vector<const edm4hep::CalorimeterHitCollection*>& caloHitColls) const override;

  // ---- Algorithm parameters ----
  Gaudi::Property<float> m_seedEnergyThreshold{
      this, "SeedEnergyThreshold", 0.020f, "Minimum crystal energy [GeV] to be considered as a Type-C seed candidate"};
  Gaudi::Property<float> m_trackWindow{this, "TrackWindow", 0.05f,
                                       "Angular search cone half-angle [rad] around the track impact point: "
                                       "sqrt(dTheta^2 + dPhi^2) < TrackWindow"};
};

#endif // TrackDrivenClusterSeeding_h
