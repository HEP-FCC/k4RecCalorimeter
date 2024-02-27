#ifndef RECFCCEECALORIMETER_CALIBRATECALOCLUSTERS_H
#define RECFCCEECALORIMETER_CALIBRATECALOCLUSTERS_H

// Key4HEP
#include "k4FWCore/DataHandle.h"

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/MsgStream.h"
class IGeoSvc;

// EDM4HEP
namespace edm4hep {
  class Cluster;
  class ClusterCollection;
  class CalorimeterHitCollection;
}

// DD4HEP
namespace dd4hep {
  namespace DDSegmentation {
    class BitFieldCoder;
  }
}

// ONNX
#include "onnxruntime/core/session/experimental_onnxruntime_cxx_api.h"

/** @class CalibrateCaloClusters
 *
 *  Apply an MVA energy calibration to the clusters reconstructed in the calorimeter.
 *
 *  @author Giovanni Marchiori
 */

class CalibrateCaloClusters : public GaudiAlgorithm {

public:
  CalibrateCaloClusters(const std::string& name, ISvcLocator* svcLoc);

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
  StatusCode readCalibrationFile(const std::string& calibrationFileName);

  /**
   * Apply the calibration to the input clusters.
   *
   * @param[in]  inClusters   Pointer to the input cluster collection.
   * @param[out] outClusters  Pointer to the output cluster collection.
   *
   * @return                  Status code.
   */
  StatusCode calibrateClusters(const edm4hep::ClusterCollection* inClusters,
                               edm4hep::ClusterCollection* outClusters);

  /**
   * Get sum of energy from cells in each layer.
   * This energy is not calibrated.
   *
   * @param[in]  cluster          Pointer to cluster of interest.
   * @param[out] energiesInLayer  Reference to vector that will contain the energies
   */
  void calcEnergiesInLayers(edm4hep::Cluster cluster,
                            std::vector<float>& energiesInLayer);


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

  /// IDs of the detectors
  Gaudi::Property<std::vector<int>> m_systemIDs {
      this, "systemIDs", {4}, "IDs of systems"
  };
  /// Names of the detector readouts, corresponding to system IDs
  Gaudi::Property<std::vector<std::string>> m_readoutNames {
      this, "readoutNames", {"ECalBarrelModuleThetaMerged"},
      "Names of the detector readout, corresponding to systemID"
  };
  /// Name of the layer field
  Gaudi::Property<std::vector<std::string>> m_layerFieldNames {
      this, "layerFieldNames", {"layer"},
      "Identifier of layers, corresponding to systemID"
  }; 
  /// Numbers of layers of the detectors
  Gaudi::Property<std::vector<size_t>> m_numLayers {
      this, "numLayers", {12}, "Numbers of layers of the systems"
  };
  /// IDs of the first layers of the detectors
  Gaudi::Property<std::vector<int>> m_firstLayerIDs {
      this, "firstLayerIDs", {0}, "IDs of first layers in the systems"
  };

  /// File with the calibration model 
  // Gaudi::Property<std::vector<std::string>> m_calibrationFiles {
  //    this, "calibrationFiles", {}, "Files with the calibration parameters"};
  Gaudi::Property<std::string> m_calibrationFile {
      this, "calibrationFile", {}, "File with the calibration parameters"};

  // total number of layers summed over the various subsystems
  // should be equal to the number of input features of the MVA
  size_t m_numLayersTotal;

  // the ONNX runtime session for applying the calibration,
  // the environment, and the input and output shapes and names
  Ort::Experimental::Session* ortSession = nullptr;
  Ort::Env* ortEnv = nullptr;
  std::vector<std::int64_t> input_shapes;
  std::vector<std::int64_t> output_shapes;
  std::vector<std::string> input_names;
  std::vector<std::string> output_names;
};

#endif /* RECFCCEECALORIMETER_CALIBRATECALOCLUSTERS_H */
