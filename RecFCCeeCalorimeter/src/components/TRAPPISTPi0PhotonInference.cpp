// k4FWCore
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/MetadataUtils.h"
#include "k4FWCore/Transformer.h"

// Gaudi
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"
#include "RecCaloCommon/ICaloCellConstantsSvc.h"

// EDM4hep
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ClusterCollection.h"

// ONNX
#include "OnnxruntimeUtilities.h"
#include "onnxruntime_cxx_api.h"
#include "onnxruntime_run_options_config_keys.h"

// STL
#include <memory>
#include <vector>

namespace {
/// An ONNX Runtime environment together with a session created from it.
/// ONNX Runtime requires the environment to outlive every session made from it, so
/// the two are kept together and shared as a pair. Members are destroyed in reverse
/// declaration order, so the session is always destroyed before the environment.
struct OnnxModel {
  Ort::Env env{ORT_LOGGING_LEVEL_WARNING, "ONNX_Runtime"};
  Ort::Session session;
  OnnxModel(const std::string& path, const Ort::SessionOptions& opts) : session(env, path.c_str(), opts) {}
};
} // namespace

/** @struct TRAPPISTPi0PhotonInference
 *
 * Gaudi MultiTransformer that produces produces a classification score for a given cluster using a
 * pre-trained ONNX model based on the GATr architecture (https://arxiv.org/abs/2305.18415). For each
 * input cluster, the tool extracts the positions and energies of its hits, formats them into tensors,
 * and feeds them into the ONNX model to obtain a classification score. The score is then appended
 * to the cluster's shapeParameters, and the new cluster is added to the output collection.
 *
 * input: edm4hep::ClusterCollection
 * output: edm4hep::ClusterCollection
 *
 *  @author Jacopo Fanini
 *  @date   2026-07
 *
 */

struct TRAPPISTPi0PhotonInference final
    : k4FWCore::MultiTransformer<std::tuple<edm4hep::ClusterCollection>(const edm4hep::ClusterCollection&)> {
public:
  TRAPPISTPi0PhotonInference(const std::string& name, ISvcLocator* svcLoc)
      : MultiTransformer(name, svcLoc, {KeyValues("inClusters", {"unpairedClusters"})},
                         {KeyValues("outClusters", {"clustersWithScore"})}) {}

  StatusCode initialize() override {
    // Initialize ONNX memory allocation
    m_memoryInfo = std::make_unique<Ort::MemoryInfo>(Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault));

    // Retrieve the shared constants service
    if (m_constantsSvc.retrieve().isFailure()) {
      error() << "Failed to retrieve CaloCellConstantsSvc" << endmsg;
      return StatusCode::FAILURE;
    }

    // Share the ONNX environment and session across instances using the same model file.
    // The first instance to initialize creates them; subsequent ones reuse them.
    const std::string svcKey = "ort_model:" + m_modelPath.value();
    if (m_constantsSvc->getObj<std::shared_ptr<OnnxModel>>(svcKey) == nullptr) {
      Ort::SessionOptions session_opts;
      session_opts.SetIntraOpNumThreads(1);
      session_opts.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_DISABLE_ALL);
      m_constantsSvc->putObj(svcKey, std::make_shared<OnnxModel>(m_modelPath.value(), session_opts));
    }
    m_ortModel = *m_constantsSvc->getObj<std::shared_ptr<OnnxModel>>(svcKey);
    Ort::Session& session = m_ortModel->session;
    // Create ONNX allocator to manage memory allocations during runtime
    Ort::AllocatorWithDefaultOptions allocator;
    // Retrieve input and output names from the ONNX model for inference
    for (std::size_t i = 0; i < session.GetInputCount(); i++) {
      m_input_names.emplace_back(session.GetInputNameAllocated(i, allocator).release());
    }
    for (std::size_t i = 0; i < session.GetOutputCount(); i++) {
      m_output_names.emplace_back(session.GetOutputNameAllocated(i, allocator).release());
    }
    // Hooks for handling collection metadata
    auto inputKey = this->inputLocations("inClusters")[0];
    auto outputKey = this->outputLocations("outClusters")[0];
    auto shapeParameterNames =
        k4FWCore::getCollectionParameter<std::vector<std::string>>(inputKey, edm4hep::labels::ShapeParameterNames, this)
            .value_or(std::vector<std::string>{});
    debug() << "Input cluster has " << shapeParameterNames.size() << " names in metadata" << endmsg;
    // Append ClusterScore to the shape parameter names
    shapeParameterNames.push_back("ClusterTRAPPISTScore");
    k4FWCore::putCollectionParameter(outputKey, edm4hep::labels::ShapeParameterNames, shapeParameterNames, this);

    return StatusCode::SUCCESS;
  }

  std::tuple<edm4hep::ClusterCollection> operator()(const edm4hep::ClusterCollection& unpairedClusters) const override {
    edm4hep::ClusterCollection outputClusters;
    debug() << "Processing " << unpairedClusters.size() << " unpaired clusters" << endmsg;
    // Loop over each cluster in the input collection and extract hits
    for (const auto& cluster : unpairedClusters) {
      const auto& hits = cluster.getHits();
      const std::size_t numHits = hits.size();
      debug() << "Processing cluster with " << numHits << " hits" << endmsg;
      std::vector<float> positions;
      std::vector<float> energies;
      positions.reserve(numHits * 3);
      energies.reserve(numHits);
      for (const auto& hit : hits) {
        const auto& pos = hit.getPosition();
        positions.push_back(pos.x);
        positions.push_back(pos.y);
        positions.push_back(pos.z);
        energies.push_back(hit.getEnergy());
      }
      // Define shape of input tensors and convert cluster data to ONNX tensors
      std::vector<int64_t> pos_shape = {static_cast<int64_t>(numHits), 3};
      std::vector<int64_t> en_shape = {static_cast<int64_t>(numHits), 1};
      std::vector<Ort::Value> input_tensors;
      input_tensors.emplace_back(vec_to_tensor<float>(positions, pos_shape, *m_memoryInfo));
      input_tensors.emplace_back(vec_to_tensor<float>(energies, en_shape, *m_memoryInfo));
      // Run ONNX inference
      auto output_tensors =
          m_ortModel->session.Run(Ort::RunOptions{nullptr}, m_input_names.data(), input_tensors.data(),
                                  input_tensors.size(), m_output_names.data(), m_output_names.size());
      if (output_tensors.empty() || !output_tensors[0].IsTensor()) {
        error() << "Failed to get a valid output tensor from ONNX session." << endmsg;
        throw std::runtime_error("ONNX output tensor is empty or not a tensor");
      }
      // Extract the two raw logits and compute sigmoid to get the probability for the pi0 class
      float* logits = output_tensors[0].GetTensorMutableData<float>();
      float logit_gamma = logits[0];
      float logit_pi0 = logits[1];
      float logits_diff = logit_pi0 - logit_gamma;
      // Sigmoid, but protected agaist overflows by using std::tanh
      float sigmoid_output = 0.5f * (1.0f + std::tanh(0.5f * logits_diff));
      debug() << "Logits: " << logit_gamma << " " << logit_pi0 << "\n"
              << "Sigmoid output: " << sigmoid_output << endmsg;
      // Create  new cluster with the computed score added to its shapeParameters
      auto outputCluster = cluster.clone();
      outputCluster.addToShapeParameters(sigmoid_output);
      outputClusters.push_back(std::move(outputCluster));
    }
    return std::make_tuple(std::move(outputClusters));
  }

  StatusCode finalize() override { return StatusCode::SUCCESS; }

private:
  // ONNX memory info: manages memory allocation for tensors during inference
  std::unique_ptr<Ort::MemoryInfo> m_memoryInfo{nullptr};
  // ONNX environment and session (loaded model used for inference).
  // Shared across instances using the same model file via ICaloCellConstantsSvc.
  std::shared_ptr<OnnxModel> m_ortModel{nullptr};
  // Handle to the shared constants service used to cache the ONNX model.
  ServiceHandle<k4::recCalo::ICaloCellConstantsSvc> m_constantsSvc{this, "CaloCellConstantsSvc",
                                                                   "k4::recCalo::CaloCellConstantsSvc", ""};
  // Input and output names for the ONNX model, used to identify tensors during inference
  std::vector<const char*> m_input_names;
  std::vector<const char*> m_output_names;
  // Path to ONNX model to be used for inference
  Gaudi::Property<std::string> m_modelPath{this, "ONNXModelPath", "",
                                           "Path to the ONNX model to be used for inference"};
};

DECLARE_COMPONENT(TRAPPISTPi0PhotonInference)