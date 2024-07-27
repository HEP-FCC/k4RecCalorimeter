#ifndef RECCALORIMETER_CORRECTCALOCLUSTERS_H
#define RECCALORIMETER_CORRECTCALOCLUSTERS_H

// Key4HEP
#include "k4FWCore/DataHandle.h"

// Gaudi
#include "Gaudi/Algorithm.h"
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
    class FCCSWGridPhiEta_k4geo;
    class FCCSWGridModuleThetaMerged_k4geo;
    class FCCSWGridPhiTheta_k4geo; 
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
 * * Benchmark calibration should be used for the combined simulation of ECal and HCal when using charged pions. 
 *    At the input level, the ECal should be calibrated to EM scale and HCal should be calibrated to HAD scale. 
 *    The aim of the benchmark calibration is to bring ECal to HAD scale and also to take into account
 *    the energy loss between the ECal and HCal (e.g. in cryostat) - for this, the energy from the last ECal layer and the first HCal layer is used. 
 *    While, downstream correction is part of the benchmark method (energy lost between ECal and HCal), as well as the upstream correction for ECal
 *    For standalone ECal, include upstream and downstream correction; for ECal+HCal simulation apply only benchmark correction should be applied
 *    To obtain the actual parameters run RecCalorimeter/tests/options/fcc_ee_caloBenchmarkCalibration.py which calls CalibrateBenchmarkMethod
 *
 *  Based on similar corrections by Jana Faltova and Anna Zaborowska.
 *
 *  @author Juraj Smiesko, benchmark calibration added by Michaela Mlynarikova
 */

class CorrectCaloClusters : public Gaudi::Algorithm {

public:
  CorrectCaloClusters(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize();

  StatusCode execute(const EventContext&) const;

  StatusCode finalize();

private:
  /**
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
                                     const std::string& funcNameStem = "upDownBenchmark");

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
   * @param[in]  systemID      ID of the system.
   * @param[in]  layerID       ID of the layer of the readout.
   *
   * @return                   Energy in layer.
   */
  double getEnergyInLayer(edm4hep::Cluster cluster,
                          const std::string& readoutName,
                          int systemID,
                          int layerID);

  /**
   * Get sum of energy from cells in the whole calorimeter.
   * This energy is not calibrated.
   *
   * @param[in]  cluster       Pointer to cluster of interest.
   * @param[in]  readoutName   Name of the readout.
   * @param[in]  systemID      ID of the system.
   *
   * @return                   Total energy.
   */
  double getTotalEnergy(edm4hep::Cluster cluster,
                        const std::string& readoutName,
                        int systemID); 

  /**
   * Get the theta angle of the specified cluster.
   *
   * @param[in]  cluster  Pointer to cluster of interest.
   *
   * @return              theta angle value.
   */
  double getClusterTheta(edm4hep::Cluster cluster);

  /// Handle for input calorimeter clusters collection
  mutable DataHandle<edm4hep::ClusterCollection> m_inClusters {
    "inClusters", Gaudi::DataHandle::Reader, this
  };
  /// Handle for corrected (output) calorimeter clusters collection
  mutable DataHandle<edm4hep::ClusterCollection> m_outClusters {
    "outClusters", Gaudi::DataHandle::Writer, this
  };

  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;
  /// Pointers to upstream correction functions
  std::vector<std::vector<TF1*>> m_upstreamFunctions;
  /// Pointers to downstream correction functions
  std::vector<std::vector<TF1*>> m_downstreamFunctions;
  /// Pointers to benchmark method correction functions
  std::vector<std::vector<TF1*>> m_benchmarkFunctions;

  /// IDs of the detectors
  Gaudi::Property<std::vector<int>> m_systemIDs {
      this, "systemIDs", {4,8}, "IDs of systems"
  };
  Gaudi::Property<uint> m_systemIDECal{this, "systemIDECal", 4, "ID of ECal system"};
  Gaudi::Property<uint> m_systemIDHCal{this, "systemIDHCal", 8, "ID of the HCal system"};
  /// Names of the detector readouts, corresponding to system IDs
  Gaudi::Property<std::vector<std::string>> m_readoutNames {
      this, "readoutNames", {"ECalBarrelModuleThetaMerged","BarHCal_Readout_phitheta"},
      "Names of the detector readout, corresponding to systemID"
  };
  /// Numbers of layers of the detectors
  Gaudi::Property<std::vector<size_t>> m_numLayers {
      this, "numLayers", {12,13}, "Numbers of layers of the systems"};
  /// IDs of the first layers of the detectors
  Gaudi::Property<std::vector<int>> m_firstLayerIDs {
      this, "firstLayerIDs", {0,0}, "IDs of first layers in the systems"
  };
  /// IDs of the last layers of the detectors
  Gaudi::Property<std::vector<int>> m_lastLayerIDs {
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
  Gaudi::Property<std::vector<std::vector<double>>> m_benchmarkParametrization {
      this, "benchmarkParametrization", {}, "Benchmark formula parameters"};
  /// Benchmark correction formulas in ROOT notation
  Gaudi::Property<std::vector<std::vector<std::string>>> m_benchmarkFormulas {
      this, "benchmarkFormulas", {}, "Benchmark formulas (in ROOT notation)"};

  /// Possibility to use different parametrization for very low energies (<10GeV) for benchmark method 
  /// For the energy (GeV) below this given threshold, the second energy formula in benchmarkFormulas is used 
  Gaudi::Property<double> m_benchmarkEneSwitch {
    this, "benchmarkEneSwitch", -1., 
    "Energy threshold in GeV to use the second formula. Set to a negative value to disable using the second formula."
};

  /// benchmarkParamsApprox are used for the first estimate of the energy 
  /// which is then used to obtain energy dependent benchmark parameters
  /// Parameters below were obtained by running a single benchmark calibration for 100 GeV pion
  Gaudi::Property<std::vector<double>> m_benchmarkParamsApprox {
      this, "benchmarkParamsApprox", {1.2281, 1, 1.07731, -0.00191622, 18.6772, 0.}, "Approximate benchmark parameters"
  };
  /// Flag if upstream correction should be applied
  Gaudi::Property<bool> m_upstreamCorr{this, "upstreamCorr", true};
  /// Flag if downstream correction should be applied
  Gaudi::Property<bool> m_downstreamCorr{this, "downstreamCorr", true};
  /// Flag if benchmark correction should be applied
  Gaudi::Property<bool> m_benchmarkCorr{this, "benchmarkCorr", false};
};

#endif /* RECCALORIMETER_CORRECTCALOCLUSTERS_H */
