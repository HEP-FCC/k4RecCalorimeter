#ifndef RECFCCEECALORIMETER_CALIBRATECALOCLUSTERS_H
#define RECFCCEECALORIMETER_CALIBRATECALOCLUSTERS_H

// Key4HEP
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/MetaDataHandle.h"

// Gaudi
#include "GaudiKernel/Algorithm.h"
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

class CalibrateCaloClusters : public Gaudi::Algorithm {

public:
  CalibrateCaloClusters(const std::string& name, ISvcLocator* svcLoc);

  virtual StatusCode initialize();

  virtual StatusCode execute(const EventContext& evtCtx) const;

  virtual StatusCode finalize();

private:
  /**
   * Initialize output calorimeter cluster collection.
   *
   * @param[in] inClusters  Pointer to the input cluster collection.
   *
   * @return                Pointer to the output cluster collection.
   */
  edm4hep::ClusterCollection* initializeOutputClusters(const edm4hep::ClusterCollection* inClusters) const;

  /**
   * Load file with MVA model into memory.
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
                               edm4hep::ClusterCollection* outClusters) const;

  /**
   * Get sum of energy from cells in each layer.
   * This energy is not calibrated.
   *
   * @param[in]  cluster          Pointer to cluster of interest.
   * @param[out] energiesInLayer  Reference to vector that will contain the energies
   */
  void calcEnergiesInLayers(edm4hep::Cluster cluster,
                            std::vector<float>& energiesInLayer) const;


  /// Handle for input calorimeter clusters collection
  mutable DataHandle<edm4hep::ClusterCollection> m_inClusters {
    "inClusters", Gaudi::DataHandle::Reader, this
  };

  /// Handle for corrected (output) calorimeter clusters collection
  mutable DataHandle<edm4hep::ClusterCollection> m_outClusters {
    "outClusters", Gaudi::DataHandle::Writer, this
  };

  /// Handles for the cluster shower shape metadata to read and to write
  MetaDataHandle<std::vector<std::string>> m_inShapeParameterHandle{
    m_inClusters,
    edm4hep::labels::ShapeParameterNames,
    Gaudi::DataHandle::Reader};
  MetaDataHandle<std::vector<std::string>> m_outShapeParameterHandle{
    m_outClusters,
    edm4hep::labels::ShapeParameterNames,
    Gaudi::DataHandle::Writer};

  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;

  /// IDs of the detectors
  Gaudi::Property<std::vector<int>> m_systemIDs {
      this, "systemIDs", {4}, "IDs of systems"
  };
  /// Name of the detectors (for the metadata)
  /// If the calibration inputs are saved in the cluster shapeParameters
  Gaudi::Property<std::vector<std::string>> m_detectorNames{
      this, "systemNames", {"EMB"}, "Names of the detectors, corresponding to systemIDs"};
  /// Names of the detector readouts, corresponding to system IDs
  /// Needed to calculate calibration inputs from cells if not
  /// present in cluster shapeParameters
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
  Gaudi::Property<std::vector<unsigned short int>> m_numLayers {
      this, "numLayers", {12}, "Numbers of layers of the systems"
  };
  /// IDs of the first layers of the detectors
  Gaudi::Property<std::vector<unsigned short int>> m_firstLayerIDs {
      this, "firstLayerIDs", {0}, "IDs of first layers in the systems"
  };

  /// File with the calibration model 
  // Gaudi::Property<std::vector<std::string>> m_calibrationFiles {
  //    this, "calibrationFiles", {}, "Files with the calibration parameters"};
  Gaudi::Property<std::string> m_calibrationFile {
      this, "calibrationFile", {}, "File with the calibration parameters"};

  // total number of layers summed over the various subsystems
  // should be equal to the number of input features of the MVA
  unsigned short int m_numLayersTotal;

  // the ONNX runtime session for applying the calibration,
  // the environment, and the input and output shapes and names
  Ort::Experimental::Session* m_ortSession = nullptr;
  Ort::Env* m_ortEnv = nullptr;
  std::vector<std::int64_t> m_input_shapes;
  std::vector<std::int64_t> m_output_shapes;
  std::vector<std::string> m_input_names;
  std::vector<std::string> m_output_names;

  // the indices of the shapeParameters containing the inputs to the model (if they exist)
  std::vector<unsigned short int> m_inputPositionsInShapeParameters;
};

#endif /* RECFCCEECALORIMETER_CALIBRATECALOCLUSTERS_H */
