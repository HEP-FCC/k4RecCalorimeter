#include "PhotonIDTool.h"

// our EDM
#include "edm4hep/Cluster.h"
#include "edm4hep/ClusterCollection.h"

#include <fstream>

DECLARE_COMPONENT(PhotonIDTool)

// convert vector data with given shape into ONNX runtime tensor
template <typename T>
Ort::Value vec_to_tensor(std::vector<T> &data, const std::vector<std::int64_t> &shape)
{
  Ort::MemoryInfo mem_info =
      Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault);
  auto tensor = Ort::Value::CreateTensor<T>(mem_info, data.data(), data.size(), shape.data(), shape.size());
  return tensor;
}

PhotonIDTool::PhotonIDTool(const std::string &name,
                           ISvcLocator *svcLoc)
    : Gaudi::Algorithm(name, svcLoc)
{
  declareProperty("inClusters", m_inClusters, "Input cluster collection");
  declareProperty("outClusters", m_outClusters, "Output cluster collection");
}

StatusCode PhotonIDTool::initialize()
{
  // Initialize base class
  {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure())
    {
      return sc;
    }
  }

  // read the files defining the model
  StatusCode sc = readMVAFiles(m_mvaInputsFile, m_mvaModelFile);
  if (sc.isFailure())
  {
    error() << "Initialization of photon ID tool config files not successful!" << endmsg;
    return sc;
  }

  // read from the metadata the names of the shape parameters in the input clusters
  std::vector<std::string> shapeParameters = m_inShapeParameterHandle.get({});
  debug() << "Variables in shapeParameters of input clusters:" << endmsg;
  for (const auto &str : shapeParameters) {
    debug() << str << endmsg;
  }

  // check if the shape parameters contain the inputs needed for the inference
  m_inputPositionsInShapeParameters.clear();
  for (const auto &feature : m_internal_input_names) {

    if (feature == "ecl") {
      // for the cluster energy, check if we have rawE in decorations
      // this is for cluster that have been passed through the MVA calibration
      // otherwise, we will use the energy of the cluster object
      auto it = std::find(shapeParameters.begin(), shapeParameters.end(), "rawE");
      if (it != shapeParameters.end())
      {
        int position = std::distance(shapeParameters.begin(), it);
        m_inputPositionsInShapeParameters.push_back(position);
        info() << "Feature " << feature << " found in position " << position << " of shapeParameters" << endmsg;
      }
      else {
        m_inputPositionsInShapeParameters.push_back(-1);
      }
    }
    else {
      // for the other features, check if they are in the shape parameters
      auto it = std::find(shapeParameters.begin(), shapeParameters.end(), feature);
      if (it != shapeParameters.end())
      {
        int position = std::distance(shapeParameters.begin(), it);
        m_inputPositionsInShapeParameters.push_back(position);
        info() << "Feature " << feature << " found in position " << position << " of shapeParameters" << endmsg;
      }
      else
      {
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

StatusCode PhotonIDTool::execute([[maybe_unused]] const EventContext &evtCtx) const
{
  verbose() << "-------------------------------------------" << endmsg;

  // Get the input collection with clusters
  const edm4hep::ClusterCollection *inClusters = m_inClusters.get();

  // Initialize output clusters
  edm4hep::ClusterCollection *outClusters = initializeOutputClusters(inClusters);
  if (!outClusters)
  {
    error() << "Something went wrong in initialization of the output cluster collection, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (inClusters->size() != outClusters->size())
  {
    error() << "Sizes of input and output cluster collections does not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }

  // Run inference
  {
    StatusCode sc = applyMVAtoClusters(inClusters, outClusters);
    if (sc.isFailure())
    {
      return sc;
    }
  }

  return StatusCode::SUCCESS;
}

StatusCode PhotonIDTool::finalize()
{
  if (m_ortSession)
    delete m_ortSession;
  if (m_ortEnv)
    delete m_ortEnv;

  return Gaudi::Algorithm::finalize();
}

edm4hep::ClusterCollection *PhotonIDTool::initializeOutputClusters(
    const edm4hep::ClusterCollection *inClusters) const
{
  edm4hep::ClusterCollection *outClusters = m_outClusters.createAndPut();

  for (auto const &inCluster : *inClusters)
  {
    auto outCluster = inCluster.clone();
    outClusters->push_back(outCluster);
  }

  return outClusters;
}

StatusCode PhotonIDTool::readMVAFiles(const std::string& mvaInputsFileName,
                                      const std::string& mvaModelFileName)
{
  // 1. read the file with the list of input features
  // TODO: figure how to save BDT to ONNX such that input is
  // the list of input features with their names, rather
  // than single tensor "X"
  std::ifstream file(mvaInputsFileName);
  if (!file.is_open()) {
    error() << "Error opening file: " << mvaInputsFileName << endmsg;
    return StatusCode::FAILURE;
  }
  
  std::string line;
  if (std::getline(file, line)) {
    // Remove square brackets
    line.erase(line.begin());
    line.erase(line.end() - 1);

    // Replace double quotes with single quotes to simplify parsing
    std::replace(line.begin(), line.end(), '"', '\'');

    std::istringstream iss(line);
    std::string token;
    while (std::getline(iss, token, ',')) {
      // Remove single quotes and leading/trailing whitespaces
      token.erase(std::remove(token.begin(), token.end(), '\''), token.end());
      token.erase(token.begin(), std::find_if(token.begin(), token.end(), [](unsigned char ch) { return !std::isspace(ch); }));
      token.erase(std::find_if(token.rbegin(), token.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(), token.end());
      m_internal_input_names.push_back(token);
    }
  }
  else {
    error() << "Could not read text from file: " << mvaInputsFileName << endmsg;
    return StatusCode::FAILURE;
  }
  debug() << "Variables needed for the photon ID MVA:" << endmsg;
  for (const auto &str : m_internal_input_names) {
    debug() << str << endmsg;
  }
  file.close();
  
  // 2. - read the file with the MVA model and setup the ONNX runtime  
  // set ONNX logging level based on output level of this alg
  OrtLoggingLevel loggingLevel = ORT_LOGGING_LEVEL_WARNING;
  MSG::Level outputLevel = this->msgStream().level();
  switch (outputLevel)
  {
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
  try
  {
    m_ortEnv = new Ort::Env(loggingLevel, "ONNX runtime environment for photonID");
    Ort::SessionOptions session_options;
    session_options.SetIntraOpNumThreads(1);
    m_ortSession = new Ort::Experimental::Session(*m_ortEnv, const_cast<std::string &>(mvaModelFileName), session_options);
    // m_ortSession = new Ort::Session(*m_ortEnv, const_cast<std::string &>(mvaModelFileName), session_options);
  }
  catch (const Ort::Exception &exception)
  {
    error() << "ERROR setting up ONNX runtime environment: " << exception.what() << endmsg;
    return StatusCode::FAILURE;
  }

  // print name/shape of inputs
  // use default allocator (CPU)
  Ort::AllocatorWithDefaultOptions allocator;
  debug() << "Input Node Name/Shape (" << m_ortSession->GetInputCount() << "):" << endmsg;
  for (std::size_t i = 0; i < m_ortSession->GetInputCount(); i++)
  {
    // for old ONNX runtime version
    // m_input_names.emplace_back(m_ortSession->GetInputName(i, allocator));
    // for new runtime version
    m_input_names.emplace_back(m_ortSession->GetInputNameAllocated(i, allocator).get());
    m_input_shapes = m_ortSession->GetInputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape(); 
    debug() << "\t" << m_input_names.at(i) << " : ";
    for (std::size_t k = 0; k < m_input_shapes.size() - 1; k++)
    {
      debug() << m_input_shapes[k] << "x";
    }
    debug() << m_input_shapes[m_input_shapes.size() - 1] << endmsg;
  }
  // some models might have negative shape values to indicate dynamic shape, e.g., for variable batch size.
  for (auto &s : m_input_shapes)
  {
    if (s < 0)
    {
      s = 1;
    }
  }

  // print name/shape of outputs
  debug() << "Output Node Name/Shape (" << m_ortSession->GetOutputCount() << "):" << endmsg;
  for (std::size_t i = 0; i < m_ortSession->GetOutputCount(); i++)
  {
    // for old ONNX runtime version
    // m_output_names.emplace_back(m_ortSession->GetOutputName(i, allocator));
    // for new runtime version
    m_output_names.emplace_back(m_ortSession->GetOutputNameAllocated(i, allocator).get());
    m_output_shapes = m_ortSession->GetOutputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape();
    debug() << m_output_shapes.size() << endmsg;
    debug() << "\t" << m_output_names.at(i) << " : ";
    for (std::size_t k = 0; k < m_output_shapes.size() - 1; k++)
    {
      debug() << m_output_shapes[k] << "x";
    }
    debug() << m_output_shapes[m_output_shapes.size() - 1] << endmsg;
  }

  debug() << "PhotonID config files read out successfully" << endmsg;
  
  return StatusCode::SUCCESS;
}

StatusCode PhotonIDTool::applyMVAtoClusters(const edm4hep::ClusterCollection *inClusters,
                                            edm4hep::ClusterCollection *outClusters) const
{
  size_t numShapeVars = m_internal_input_names.size();
  std::vector<float> mvaInputs(numShapeVars);
  
  // loop over the input clusters and perform the inference
  for (unsigned int j = 0; j < inClusters->size(); ++j)
  {
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
    for (unsigned short int k = 0; k < numShapeVars; ++k)
    {
      verbose() << "var " << k << " : " << mvaInputs[k] << endmsg;
    }

    // run the MVA and save the output score in output
    float score= -1.0;
    // Create a single Ort tensor
    std::vector<Ort::Value> input_tensors;
    input_tensors.emplace_back(vec_to_tensor<float>(mvaInputs, m_input_shapes));

    // pass data through model
    try
    {
      std::vector<Ort::Value> output_tensors = m_ortSession->Run(m_input_names,
                                                                 input_tensors,
                                                                 m_output_names,
                                                                 Ort::RunOptions{nullptr});

      // double-check the dimensions of the output tensors
      // NOTE: the number of output tensors is equal to the number of output nodes specified in the Run() call
      // assert(output_tensors.size() == output_names.size() && output_tensors[0].IsTensor());
      // the probabilities are in the 2nd entry of the output
      debug() << output_tensors.size() << endmsg;
      debug() << output_tensors[1].GetTensorTypeAndShapeInfo().GetShape() << endmsg;
      float *outputData = output_tensors[1].GetTensorMutableData<float>();
      for (int i=0; i<2; i++)
        debug() << i << " " << outputData[i] << endmsg;
      score = outputData[1];
    }
    catch (const Ort::Exception &exception)
    {
      error() << "ERROR running model inference: " << exception.what() << endmsg;
      return StatusCode::FAILURE;
    }

    verbose() << "Photon ID score: " << score << endmsg;
    outClusters->at(j).addToShapeParameters(score);
  }

  return StatusCode::SUCCESS;
}
