#ifndef RECCALORIMETER_CORRECTCALOCLUSTERS_H
#define RECCALORIMETER_CORRECTCALOCLUSTERS_H

// Key4HEP
#include "k4FWCore/DataHandle.h"

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/RndmGenerators.h"
class IGeoSvc;
class IRndmGenSvc;
class ITHistSvc;

// ROOT
#include "TF2.h"

// EDM4HEP
namespace edm4hep {
  class Cluster;
  class MutableCluster;
  class ClusterCollection;
  class CalorimeterHitCollection;
  class MCParticleCollection;
  class VertexCollection;
}

// DD4HEP
namespace dd4hep {
  namespace DDSegmentation {
    class FCCSWGridPhiEta;
    class MultiSegmentation;
    class BitFieldCoder;
  }
}

/** @class CorrectCaloClusters
 *
 *  Apply corrections to the clusters reconstructed in ECAL.
 *  * Upstream energy correction corrects for the energy lost by the particles in the material before they enter
 *    active volume of the calorimeter. The correction is parametrized in one (cluster energy) or two variables (cluster
 *    energy, cluster angle)
 *  * Downstream energy correction corrects for the energy lost due to punch trough. This energy is deposited in the
 *    instrumentation behind the active volume of the calorimeter. It is parametrized in one (cluster energy) or two
 *    variables (cluster energy, cluster angle)
 *  * Benchmark method energy correction to be used for ECal+HCal simulation when shooting hadrons. 
 *    This correction brings ECal to HAD scale (HCal is expected to be at HAD scale at the input level) and also corrects 
 *    for the amount of energy lost between ECal and HCal, for the non-compensation of the ECal 
 *
 *  Based on similar corrections by Jana Faltova and Anna Zaborowska.
 *
 *  @author Juraj Smiesko, March 2023: M. Mlynarikova added benchmark method
 */

class CorrectCaloClusters : public GaudiAlgorithm {

public:
  CorrectCaloClusters(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  /**
   * 
   * Initialize output calorimeter cluster collection.
   *
   * @param[in] inClusters  Pointer to the input cluster collection.
   *
   * @return                Pointer to the output cluster collection.
   */
  edm4hep::ClusterCollection* initializeOutputClusters(const edm4hep::ClusterCollection* inClusters);

  /**
   * Initialize vectors of upstream and downstream correction functions.
   *
   * @return                  Status code.
   */
  StatusCode initializeCorrFunctions(std::vector<std::vector<TF1*>>& functions,
                                     std::vector<std::vector<std::string>> formulas,
                                     std::vector<std::vector<double>> parameters,
                                     const std::string& funcNameStem = "upDown");

  /**
   * Initialize vectors of benchmark correction functions.
   *
   * @return                  Status code.
   */
  StatusCode initializeBenchmarkCorrFunctions(std::vector<TF1*>& functions,
                                              std::vector<std::string> formulas,
                                              std::vector<double> parameters,
                                              const std::string& funcNameStem = "benchmark");

  /**
   * Apply upstream correction to the output clusters.
   *
   * @param[in]  inClusters   Pointer to the input cluster collection.
   * @param[out] outClusters  Pointer to the output cluster collection.
   *
   * @return                  Status code.
   */
  StatusCode applyUpstreamCorr(const edm4hep::ClusterCollection* inClusters,
                               edm4hep::ClusterCollection* outClusters);

  /**
   * Apply downstream correction to the output clusters.
   *
   * @param[in]  inClusters   Pointer to the input cluster collection.
   * @param[out] outClusters  Pointer to the output cluster collection.
   *
   * @return                  Status code.
   */
  StatusCode applyDownstreamCorr(const edm4hep::ClusterCollection* inClusters,
                                 edm4hep::ClusterCollection* outClusters);

  /**
   * Apply benchmark correction to the output clusters.
   *
   * @param[in]  inClusters   Pointer to the input cluster collection.
   * @param[out] outClusters  Pointer to the output cluster collection.
   *
   * @return                  Status code.
   */
  StatusCode applyBenchmarkCorr(const edm4hep::ClusterCollection* inClusters,
                                 edm4hep::ClusterCollection* outClusters);

  /**
   * Get sum of energy from cells in specified layer.
   * This energy is not calibrated.
   *
   * @param[in]  cluster       Pointer to cluster of interest.
   * @param[in]  readoutName   Name of the readout.
   * @param[in]  layerID       ID of the layer of the readout.
   *
   * @return                   Energy in layer.
   */
  double getEnergyInLayer(edm4hep::Cluster cluster,
                          const std::string& readoutName,
                          size_t systemID,
                          size_t layerID);

  /**
   * Get the theta angle of the specified cluster.
   *
   * @param[in]  cluster  Pointer to cluster of interest.
   *
   * @return              theta angle value.
   */
  double getClusterTheta(edm4hep::Cluster cluster);

  /// Handle for input calorimeter clusters collection
  DataHandle<edm4hep::ClusterCollection> m_inClusters {
    "inClusters", Gaudi::DataHandle::Reader, this
  };
  /// Handle for corrected (output) calorimeter clusters collection
  DataHandle<edm4hep::ClusterCollection> m_outClusters {
    "outClusters", Gaudi::DataHandle::Writer, this
  };

  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;
  /// Pointers to upstream correction functions
  std::vector<std::vector<TF1*>> m_upstreamFunctions;
  /// Pointers to downstream correction functions
  std::vector<std::vector<TF1*>> m_downstreamFunctions;
  /// Pointers to benchmark method correction functions
  std::vector<TF1*> m_benchmarkFunctions;

  /// IDs of the detectors
  Gaudi::Property<std::vector<size_t>> m_systemIDs {
      this, "systemIDs", {4,8}, "IDs of systems"
  };
  Gaudi::Property<uint> m_ECalSystemID{this, "ECalSystemID", 4, "ID of ECal system"};
  Gaudi::Property<uint> m_HCalSystemID{this, "HCalSystemID", 8, "ID of the HCal system"};
  /// Names of the detector readouts, corresponding to system IDs
  Gaudi::Property<std::vector<std::string>> m_readoutNames {
      this, "readoutNames", {"ECalBarrelPhiEta","HCalBarrelReadout"},
      "Names of the detector readout, corresponding to systemID"
  };
  /// Numbers of layers of the detectors
  Gaudi::Property<std::vector<size_t>> m_numLayers {
      this, "numLayers", {12,13}, "Numbers of layers of the systems"};
  /// IDs of the first layers of the detectors
  Gaudi::Property<std::vector<size_t>> m_firstLayerIDs {
      this, "firstLayerIDs", {0,0}, "IDs of first layers in the systems"
  };
  /// IDs of the last layers of the detectors
  Gaudi::Property<std::vector<size_t>> m_lastLayerIDs {
      this, "lastLayerIDs", {11,12}, "IDs of last layers in the systems"
  };

  /// Upstream correction parameters (a, b, c, ...)
  Gaudi::Property<std::vector<std::vector<double>>> m_upstreamParams {
      this, "upstreamParameters", {}, "Upstream parameters"};
  /// Upstream correction formulas in ROOT notation
  Gaudi::Property<std::vector<std::vector<std::string>>> m_upstreamFormulas {
      this, "upstreamFormulas", {}, "Upstream formulas (in ROOT notation)"};

  /// Downstream correction parameters (a, b, c, ...)
  Gaudi::Property<std::vector<std::vector<double>>> m_downstreamParams {
      this, "downstreamParameters", {}, "Downstream parameters"};
  /// Downstream correction formulas in ROOT notation
  Gaudi::Property<std::vector<std::vector<std::string>>> m_downstreamFormulas {
      this, "downstreamFormulas", {}, "Downstream formulas (in ROOT notation)"};

  /// Benchmark method correction parameters (a, b, c, ...)
  Gaudi::Property<std::vector<double>> m_benchmarkParams {
      this, "benchmarkParameters", {}, "Benchmark parameters"};
  /// Downstream correction formulas in ROOT notation
  Gaudi::Property<std::vector<std::string>> m_benchmarkFormulas {
      this, "benchmarkFormulas", {}, "Benchmark formulas (in ROOT notation)"};

  /// benchmarkParamsApprox are used for the first estimate of the energy, below parameters for 100 GeV pion
  Gaudi::Property<std::vector<double>> m_benchmarkParamsApprox {
      this, "benchmarkParamsApprox", {1.2909, 1, 0.91, -0.0019}, "Approximate benchmark parameters"};

  /// Flag if upstream correction should be applied
  Gaudi::Property<bool> m_upstreamCorr{this, "upstreamCorr", true};
  /// Flag if downstream correction should be applied
  Gaudi::Property<bool> m_downstreamCorr{this, "downstreamCorr", true};
  /// Flag if benchmark correction should be applied
  Gaudi::Property<bool> m_benchmarkCorr{this, "benchmarkCorr", true};
};

#endif /* RECCALORIMETER_CORRECTCALOCLUSTERS_H */
