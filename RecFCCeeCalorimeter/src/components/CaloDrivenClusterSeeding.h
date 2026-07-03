#ifndef CaloDrivenClusterSeeding_h
#define CaloDrivenClusterSeeding_h 1

/*
 * CaloDrivenClusterSeeding
 *
 * Gaudi MultiTransformer that identifies ECAL cluster seeds purely from
 * calorimeter hit energy deposits (no tracking information).
 *
 * Input
 * -----
 *   std::vector<const edm4hep::CalorimeterHitCollection*>
 *       All digitised calorimeter hit collections.  The algorithm filters
 *       internally for (ECAL) depth-0 hits
 *
 * Output  std::tuple<seedsA, seedsB>
 * ------
 *   seedsA  (edm4hep::ClusterCollection)
 *       Type-A seeds: crystals with E > SeedEnergyThresholdA that are a
 *       local maximum within their Von Neumann distance-k neighbourhood.
 *
 *   seedsB  (edm4hep::ClusterCollection)
 *       Type-B seeds: crystals with E > SeedEnergyThresholdB that are a
 *       local maximum within their Von Neumann distance-k neighbourhood,
 *       have at least MinAboveThresholdNeighbours (including self) above the
 *       above-threshold neighbours in that neighbourhood, AND all
 *       above-threshold crystals in the neighbourhood are topologically
 *       connected via VN-d1 adjacency.
 */

#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ClusterCollection.h"

#include "k4FWCore/Transformer.h"
#include "Gaudi/Property.h"

#include "ClusterSeedingBase.h"

#include <tuple>
#include <vector>

struct CaloDrivenClusterSeeding final
    : ClusterSeedingBase<
          std::tuple<edm4hep::ClusterCollection, edm4hep::ClusterCollection>(
              const std::vector<const edm4hep::CalorimeterHitCollection*>&)> {

  CaloDrivenClusterSeeding(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode finalize() override { return StatusCode::SUCCESS; }

  std::tuple<edm4hep::ClusterCollection, edm4hep::ClusterCollection>
  operator()(const std::vector<const edm4hep::CalorimeterHitCollection*>& caloHitColls) const override;

  // ---- Configurable parameters ----
  Gaudi::Property<float> m_seedEnergyThresholdA{
      this, "SeedEnergyThresholdA", 0.040f,
      "Minimum crystal energy [GeV] for a Type-A seed (local max, Von Neumann neighbourhood)"};
  Gaudi::Property<float> m_seedEnergyThresholdB{
      this, "SeedEnergyThresholdB", 0.020f,
      "Minimum crystal energy [GeV] for a Type-B seed (local max + connectivity)"};
  Gaudi::Property<unsigned int> m_minAboveThresholdNeighbours{
      this, "MinAboveThresholdNeighbours", 2,
      "Minimum number of above-threshold neighbours (including seed) for a Type-B seed"};
  Gaudi::Property<int> m_vnDistance{this, "VonNeumannDistance", 2,
      "Von Neumann distance for neighbourhood definition (applies to both seed types)"};
};

#endif // CaloDrivenClusterSeeding_h
