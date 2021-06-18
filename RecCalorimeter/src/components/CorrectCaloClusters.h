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
 *
 *  Based on similar corrections by Jana Faltova and Anna Zaborowska.
 *
 *  @author Juraj Smiesko
 */

class CorrectCaloClusters : public GaudiAlgorithm {

public:
  CorrectCaloClusters(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize();

  StatusCode execute();

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
                                     const std::string& funcNameStem = "upDown");

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
   * Get sum of energy from cells in specified layer.
   * This energy is not calibrated.
   *
   * @param[in]  cluster       Pointer to cluster of interest.
   * @param[in]  readoutName   Name of the readout.
   * @param[in]  layerID       ID of the layer of the readout.
   *
   * @return                   Energy in layer.
   */
  double getEnergyInLayer(const edm4hep::Cluster& cluster,
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
  double getClusterTheta(const edm4hep::Cluster& cluster);

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

  /// IDs of the detectors
  Gaudi::Property<std::vector<size_t>> m_systemIDs {
      this, "systemIDs", {4}, "IDs of systems"
  };
  /// Names of the detector readouts, corresponding to system IDs
  Gaudi::Property<std::vector<std::string>> m_readoutNames {
      this, "readoutNames", {"ECalBarrelPhiEta"},
      "Names of the detector readout, corresponding to systemID"
  };
  /// Numbers of layers of the detectors
  Gaudi::Property<std::vector<size_t>> m_numLayers {
      this, "numLayers", {12}, "Numbers of layers of the systems"};
  /// IDs of the first layers of the detectors
  Gaudi::Property<std::vector<size_t>> m_firstLayerIDs {
      this, "firstLayerIDs", {0}, "IDs of first layers in the systems"
  };
  /// IDs of the last layers of the detectors
  Gaudi::Property<std::vector<size_t>> m_lastLayerIDs {
      this, "lastLayerIDs", {7}, "IDs of last layers in the systems"
  };
  /// Values of sampling fractions used for energy calibration of the systems
  Gaudi::Property<std::vector<std::vector<double>>> m_samplingFractions {
      this, "samplingFractions",
      {{0.299041341789, 0.1306220735, 0.163243999965, 0.186360269398,
        0.203778124831, 0.216211280314, 0.227140796653, 0.243315422934}},
      "Values of sampling fractions used in energy calibration of the systems"};

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
};

#endif /* RECCALORIMETER_CORRECTCALOCLUSTERS_H */
