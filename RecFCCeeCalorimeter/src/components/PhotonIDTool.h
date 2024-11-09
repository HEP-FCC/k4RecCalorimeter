#ifndef RECFCCEECALORIMETER_PHOTONIDTOOL_H
#define RECFCCEECALORIMETER_PHOTONIDTOOL_H


#include "edm4hep/Constants.h"

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
}

// ONNX
#include "onnxruntime_cxx_api.h"

/** @class PhotonIDTool
 *
 *  Apply a binary MVA classifier to discriminate between photons and pi0s.
 *  It takes a cluster collection in inputs, runs the inference using as inputs
 *  the variables in the shapeParameters of the input clusters, decorates the
 *  cluster with the photon probability (appended to the shapeParameters vector)
 *  and saves the cluster in a new output collection.
 *
 *  @author Giovanni Marchiori
 */

class PhotonIDTool : public Gaudi::Algorithm {

public:
  PhotonIDTool(const std::string& name, ISvcLocator* svcLoc);

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
  StatusCode readMVAFiles(const std::string& mvaInputsFileName,
			  const std::string& mvaModelFileName);

  /**
   * Calculate the MVA score for the input clusters and adds it
   * as a new shapeParameter of the output clusters
   *
   * @param[in]  inClusters   Pointer to the input cluster collection.
   * @param[out] outClusters  Pointer to the output cluster collection.
   *
   * @return                  Status code.
   */
  StatusCode applyMVAtoClusters(const edm4hep::ClusterCollection* inClusters,
				edm4hep::ClusterCollection* outClusters) const;

  /// Handle for input calorimeter clusters collection
  mutable DataHandle<edm4hep::ClusterCollection> m_inClusters {
    "inClusters", Gaudi::DataHandle::Reader, this
  };

  /// Handle for output calorimeter clusters collection
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

  /// Files with the MVA model and list of inputs
  Gaudi::Property<std::string> m_mvaModelFile {
      this, "mvaModelFile", {}, "ONNX file with the mva model"};
  Gaudi::Property<std::string> m_mvaInputsFile {
      this, "mvaInputsFile", {}, "JSON file with the mva inputs"};

  // the ONNX runtime session for running the inference,
  // the environment, and the input and output shapes and names
  Ort::Session* m_ortSession = nullptr;
  Ort::Env* m_ortEnv = nullptr;
  Ort::MemoryInfo m_ortMemInfo;
  std::vector<std::int64_t> m_input_shapes;
  std::vector<std::int64_t> m_output_shapes;
  std::vector<const char*> m_input_names;
  std::vector<const char*> m_output_names;
  std::vector<std::string> m_internal_input_names;
  
  // the indices of the shapeParameters containing the inputs to the model (-1 if not found)
  std::vector<short int> m_inputPositionsInShapeParameters;
};

#endif /* RECFCCEECALORIMETER_PHOTONIDTOOL_H */
