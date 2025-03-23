#include "PhotonIDTool.h"

// our EDM
#include "edm4hep/Cluster.h"
#include "edm4hep/ClusterCollection.h"

#include <fstream>

#include "nlohmann/json.hpp"

#include "OnnxruntimeUtilities.h"

using json = nlohmann::json;

DECLARE_COMPONENT(PhotonIDTool)

PhotonIDTool::PhotonIDTool(const std::string& name, ISvcLocator* svcLoc)
    : Gaudi::Algorithm(name, svcLoc),
      m_ortMemInfo(Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault)) {
  declareProperty("inClusters", m_inClusters, "Input cluster collection");
  declareProperty("outClusters", m_outClusters, "Output cluster collection");
}

StatusCode PhotonIDTool::initialize() {
  // Initialize base class
  {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure()) {
      return sc;
    }
  }

  // read the files defining the model
  StatusCode sc = readMVAFiles(m_mvaInputsFile, m_mvaModelFile);
  if (sc.isFailure()) {
    error() << "Initialization of photon ID tool config files not successful!" << endmsg;
    return sc;
  }

  // read from the metadata the names of the shape parameters in the input clusters
  std::vector<std::string> shapeParameters = m_inShapeParameterHandle.get({});
  debug() << "Variables in shapeParameters of input clusters:" << endmsg;
  for (const auto& str : shapeParameters) {
    debug() << str << endmsg;
  }

  // check if the shape parameters contain the inputs needed for the inference
  m_inputPositionsInShapeParameters.clear();
  for (const auto& feature : m_internal_input_names) {

    if (feature == "ecl") {
      // for the cluster energy, check if we have rawE in decorations
      // this is for cluster that have been passed through the MVA calibration
      // otherwise, we will use the energy of the cluster object
      auto it = std::find(shapeParameters.begin(), shapeParameters.end(), "rawE");
      if (it != shapeParameters.end()) {
        int position = std::distance(shapeParameters.begin(), it);
        m_inputPositionsInShapeParameters.push_back(position);
        info() << "Feature " << feature << " found in position " << position << " of shapeParameters" << endmsg;
      } else {
        m_inputPositionsInShapeParameters.push_back(-1);
      }
    } else {
      // for the other features, check if they are in the shape parameters
      auto it = std::find(shapeParameters.begin(), shapeParameters.end(), feature);
      if (it != shapeParameters.end()) {
        int position = std::distance(shapeParameters.begin(), it);
        m_inputPositionsInShapeParameters.push_back(position);
        info() << "Feature " << feature << " found in position " << position << " of shapeParameters" << endmsg;
      } else {
        // at least one of the inputs of the MVA was not found in the shapeParameters
        // so we can stop checking the others
        m_inputPositionsInShapeParameters.clear();
        error() << "Feature " << feature << " not found, aborting..." << endmsg;
        return StatusCode::FAILURE;
      }
    }
  }

  // append the MVA score to the output shape parameters
  shapeParameters.push_back("photonIDscore");
  m_outShapeParameterHandle.put(shapeParameters);

  info() << "Initialized the photonID MVA tool" << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode PhotonIDTool::execute([[maybe_unused]] const EventContext& evtCtx) const {
  verbose() << "-------------------------------------------" << endmsg;

  // Get the input collection with clusters
  const edm4hep::ClusterCollection* inClusters = m_inClusters.get();

  // Initialize output clusters
  edm4hep::ClusterCollection* outClusters = initializeOutputClusters(inClusters);
  if (!outClusters) {
    error() << "Something went wrong in initialization of the output cluster collection, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (inClusters->size() != outClusters->size()) {
    error() << "Sizes of input and output cluster collections does not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }

  // Run inference
  {
    StatusCode sc = applyMVAtoClusters(inClusters, outClusters);
    if (sc.isFailure()) {
      return sc;
    }
  }

  return StatusCode::SUCCESS;
}

StatusCode PhotonIDTool::finalize() {
  if (m_ortSession)
    delete m_ortSession;

  if (m_ortEnv)
    delete m_ortEnv;

  for (auto& name : m_input_names) {
    delete name;
  }

  for (auto& name : m_output_names) {
    delete name;
  }

  return Gaudi::Algorithm::finalize();
}

edm4hep::ClusterCollection* PhotonIDTool::initializeOutputClusters(const edm4hep::ClusterCollection* inClusters) const {
  edm4hep::ClusterCollection* outClusters = m_outClusters.createAndPut();

  for (auto const& inCluster : *inClusters) {
    auto outCluster = inCluster.clone();
    outClusters->push_back(outCluster);
  }

  return outClusters;
}

StatusCode PhotonIDTool::readMVAFiles(const std::string& mvaInputsFileName, const std::string& mvaModelFileName) {
  // 1. read the file with the list of input features
  // Open the JSON file
  std::ifstream file(mvaInputsFileName);
  if (!file.is_open()) {
    error() << "Error opening file: " << mvaInputsFileName << endmsg;
    return StatusCode::FAILURE;
  }

  // Parse the JSON file
  json j;
  try {
    file >> j;
  } catch (const nlohmann::json::exception& e) {
    error() << "Error parsing JSON: " << e.what() << endmsg;
    return StatusCode::FAILURE;
  }
  file.close();

  // Access the data and print to screen
  std::string timeStamp;
  if (!j.contains("timeStamp")) {
    error() << "Error: timeStamp key not found in JSON" << endmsg;
    return StatusCode::FAILURE;
  } else {
    timeStamp = j["timeStamp"];
  }

  std::string clusterCollection;
  if (!j.contains("clusterCollection")) {
    error() << "Error: clusterCollection key not found in JSON" << endmsg;
    return StatusCode::FAILURE;
  } else {
    clusterCollection = j["clusterCollection"];
  }

  std::string trainingTool;
  if (!j.contains("trainingTool")) {
    error() << "Error: trainingTool key not found in JSON" << endmsg;
    return StatusCode::FAILURE;
  } else {
    trainingTool = j["trainingTool"];
  }

  info() << "Using the following photon-ID training:" << endmsg;
  info() << "   Timestamp: " << timeStamp << endmsg;
  info() << "   Training tool used: " << trainingTool << endmsg;
  info() << "   Input cluster collection: " << clusterCollection << endmsg;
  if (!j.contains("shapeParameters")) {
    error() << "Error: shapeParameters key not found in JSON" << endmsg;
    return StatusCode::FAILURE;
  } else {
    try {
      const auto& shape_params = j["shapeParameters"];
      if (!shape_params.is_array()) {
        throw std::runtime_error("shapeParameters is not an array");
      }
      for (const auto& param : shape_params) {
        if (!param.is_string()) {
          throw std::runtime_error("shapeParameters contains non-string values");
        }
        m_internal_input_names.push_back(param.get<std::string>());
      }
    } catch (const std::exception& e) {
      error() << "Error: " << e.what() << endmsg;
      return StatusCode::FAILURE;
    }
  }
  info() << "   Input shape parameters:" << endmsg;
  for (const auto& str : m_internal_input_names) {
    info() << "      " << str << endmsg;
  }
  if (!j.contains("trainingParameters")) {
    error() << "Error: trainingParameters key not found in JSON" << endmsg;
    return StatusCode::FAILURE;
  } else {
    info() << "   Training parameters:" << endmsg;
    for (const auto& param : j["trainingParameters"].items()) {
      std::string key = param.key();
      std::string value;
      if (param.value().is_string()) {
        value = param.value().get<std::string>();
      } else if (param.value().is_number()) {
        value = std::to_string(param.value().get<double>());
      } else if (param.value().is_null()) {
        value = "null";
      } else {
        value = "invalid";
      }
      info() << "      " << key << " : " << value << endmsg;
    }
  }

  // 2. - read the file with the MVA model and setup the ONNX runtime
  // set ONNX logging level based on output level of this alg
  OrtLoggingLevel loggingLevel = ORT_LOGGING_LEVEL_WARNING;
  MSG::Level outputLevel = this->msgStream().level();
  switch (outputLevel) {
  case MSG::Level::FATAL:                   // 6
    loggingLevel = ORT_LOGGING_LEVEL_FATAL; // 4
    break;
  case MSG::Level::ERROR:                   // 5
    loggingLevel = ORT_LOGGING_LEVEL_ERROR; // 3
    break;
  case MSG::Level::WARNING:                   // 4
    loggingLevel = ORT_LOGGING_LEVEL_WARNING; // 2
    break;
  case MSG::Level::INFO:                      // 3
    loggingLevel = ORT_LOGGING_LEVEL_WARNING; // 2 (ORT_LOGGING_LEVEL_INFO too verbose..)
    break;
  case MSG::Level::DEBUG:                  // 2
    loggingLevel = ORT_LOGGING_LEVEL_INFO; // 1
    break;
  case MSG::Level::VERBOSE:                   // 1
    loggingLevel = ORT_LOGGING_LEVEL_VERBOSE; // 0
    break;
  default:
    break;
  }
  try {
    m_ortEnv = new Ort::Env(loggingLevel, "ONNX runtime environment for photonID");
    Ort::SessionOptions session_options;
    session_options.SetIntraOpNumThreads(1);
    m_ortSession = new Ort::Session(*m_ortEnv, mvaModelFileName.data(), session_options);
  } catch (const Ort::Exception& exception) {
    error() << "ERROR setting up ONNX runtime environment: " << exception.what() << endmsg;
    return StatusCode::FAILURE;
  }

  // print name/shape of inputs
  // use default allocator (CPU)
  Ort::AllocatorWithDefaultOptions allocator;
#if ORT_API_VERSION < 13
  // Before 1.13 we have to roll our own unique_ptr wrapper here
  auto allocDeleter = [&allocator](char* p) { allocator.Free(p); };
  using AllocatedStringPtr = std::unique_ptr<char, decltype(allocDeleter)>;
#endif

  debug() << "Input Node Name/Shape (" << m_ortSession->GetInputCount() << "):" << endmsg;
  for (std::size_t i = 0; i < m_ortSession->GetInputCount(); i++) {
#if ORT_API_VERSION < 13
    m_input_names.emplace_back(AllocatedStringPtr(m_ortSession->GetInputName(i, allocator), allocDeleter).release());
#else
    m_input_names.emplace_back(m_ortSession->GetInputNameAllocated(i, allocator).release());
#endif

    m_input_shapes = m_ortSession->GetInputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape();
    debug() << "\t" << m_input_names.at(i) << " : ";
    for (std::size_t k = 0; k < m_input_shapes.size() - 1; k++) {
      debug() << m_input_shapes[k] << "x";
    }
    debug() << m_input_shapes[m_input_shapes.size() - 1] << endmsg;
  }
  // some models might have negative shape values to indicate dynamic shape, e.g., for variable batch size.
  for (auto& s : m_input_shapes) {
    if (s < 0) {
      s = 1;
    }
  }

  // print name/shape of outputs
  debug() << "Output Node Name/Shape (" << m_ortSession->GetOutputCount() << "):" << endmsg;
  for (std::size_t i = 0; i < m_ortSession->GetOutputCount(); i++) {
#if ORT_API_VERSION < 13
    m_output_names.emplace_back(AllocatedStringPtr(m_ortSession->GetOutputName(i, allocator), allocDeleter).release());
#else
    m_output_names.emplace_back(m_ortSession->GetOutputNameAllocated(i, allocator).release());
#endif

    m_output_shapes = m_ortSession->GetOutputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape();
    debug() << m_output_shapes.size() << endmsg;
    debug() << "\t" << m_output_names.at(i) << " : ";
    for (std::size_t k = 0; k < m_output_shapes.size() - 1; k++) {
      debug() << m_output_shapes[k] << "x";
    }
    debug() << m_output_shapes[m_output_shapes.size() - 1] << endmsg;
  }

  debug() << "PhotonID config files read out successfully" << endmsg;

  return StatusCode::SUCCESS;
}

StatusCode PhotonIDTool::applyMVAtoClusters(const edm4hep::ClusterCollection* inClusters,
                                            edm4hep::ClusterCollection* outClusters) const {
  size_t numShapeVars = m_internal_input_names.size();
  std::vector<float> mvaInputs(numShapeVars);

  // loop over the input clusters and perform the inference
  for (unsigned int j = 0; j < inClusters->size(); ++j) {
    // read the values of the input features
    for (unsigned int i = 0; i < m_inputPositionsInShapeParameters.size(); i++) {
      int position = m_inputPositionsInShapeParameters[i];
      if (position == -1)
        mvaInputs[i] = (inClusters->at(j)).getEnergy();
      else
        mvaInputs[i] = (inClusters->at(j)).getShapeParameters(position);
    }

    // print the values of the input features
    verbose() << "MVA inputs:" << endmsg;
    for (unsigned short int k = 0; k < numShapeVars; ++k) {
      verbose() << "var " << k << " : " << mvaInputs[k] << endmsg;
    }

    // run the MVA and save the output score in output
    float score = -1.0;
    // Create a single Ort tensor
    std::vector<Ort::Value> input_tensors;
    input_tensors.emplace_back(vec_to_tensor<float>(mvaInputs, m_input_shapes, m_ortMemInfo));

    // pass data through model
    try {
      auto output_tensors = m_ortSession->Run(Ort::RunOptions{nullptr}, m_input_names.data(), input_tensors.data(),
                                              input_tensors.size(), m_output_names.data(), m_output_names.size());

      // double-check the dimensions of the output tensors
      // NOTE: the number of output tensors is equal to the number of output nodes specified in the Run() call
      // assert(output_tensors.size() == output_names.size() && output_tensors[0].IsTensor());
      // the probabilities are in the 2nd entry of the output
      debug() << output_tensors.size() << endmsg;
      debug() << output_tensors[1].GetTensorTypeAndShapeInfo().GetShape() << endmsg;
      float* outputData = output_tensors[1].GetTensorMutableData<float>();
      for (int i = 0; i < 2; i++)
        debug() << i << " " << outputData[i] << endmsg;
      score = outputData[1];
    } catch (const Ort::Exception& exception) {
      error() << "ERROR running model inference: " << exception.what() << endmsg;
      return StatusCode::FAILURE;
    }

    verbose() << "Photon ID score: " << score << endmsg;
    outClusters->at(j).addToShapeParameters(score);
  }

  return StatusCode::SUCCESS;
}
