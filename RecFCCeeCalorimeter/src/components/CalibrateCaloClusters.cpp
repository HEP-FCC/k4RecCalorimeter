#include "CalibrateCaloClusters.h"

// Key4HEP
#include "k4Interface/IGeoSvc.h"

// FCC Detectors
#include "detectorCommon/DetUtils_k4geo.h"

// DD4hep
#include "DD4hep/Detector.h"

// our EDM
#include "edm4hep/Cluster.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"

DECLARE_COMPONENT(CalibrateCaloClusters)

// convert vector data with given shape into ONNX runtime tensor
template <typename T>
Ort::Value vec_to_tensor(std::vector<T> &data, const std::vector<std::int64_t> &shape)
{
  Ort::MemoryInfo mem_info =
      Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault);
  auto tensor = Ort::Value::CreateTensor<T>(mem_info, data.data(), data.size(), shape.data(), shape.size());
  return tensor;
}

CalibrateCaloClusters::CalibrateCaloClusters(const std::string &name,
                                             ISvcLocator *svcLoc)
    : Gaudi::Algorithm(name, svcLoc),
      m_geoSvc("GeoSvc", "CalibrateCaloClusters")
{
  declareProperty("inClusters", m_inClusters,
                  "Input cluster collection");
  declareProperty("outClusters", m_outClusters,
                  "Calibrated (output) cluster collection");
}

StatusCode CalibrateCaloClusters::initialize()
{
  // Initialize base class
  {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure())
    {
      return sc;
    }
  }

  // Check if readouts exist
  for (unsigned short int i = 0; i < m_readoutNames.size(); ++i)
  {
    auto readouts = m_geoSvc->getDetector()->readouts();
    if (readouts.find(m_readoutNames.value().at(i)) == readouts.end())
    {
      error() << "Missing readout, exiting!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // Check if readout related variables have the same size
  if (m_systemIDs.size() != m_readoutNames.size())
  {
    error() << "Sizes of the systemIDs vector and readoutNames vector do not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_numLayers.size())
  {
    error() << "Sizes of systemIDs vector and numLayers vector do not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_layerFieldNames.size())
  {
    error() << "Sizes of systemIDs vector and layerFieldNames vector do not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_firstLayerIDs.size())
  {
    error() << "Sizes of systemIDs vector and firstLayerIDs vector do not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }

  // calculate total number of layers summed over the various subsystems
  m_numLayersTotal = 0;
  for (unsigned short int i = 0; i < m_numLayers.size(); ++i)
  {
    m_numLayersTotal += m_numLayers[i];
  }

  // Initialize the calibration tool
  StatusCode sc = readCalibrationFile(m_calibrationFile);
  if (sc.isFailure())
  {
    error() << "Initialization of calibration tool correction functions not successful!" << endmsg;
    return sc;
  }

  // read from the metadata the names of the shape parameters in the input clusters and append the total raw energy to the output
  std::vector<std::string> shapeParameters = m_inShapeParameterHandle.get();
  shapeParameters.push_back("rawE");
  m_outShapeParameterHandle.put(shapeParameters);

  // check if the shape parameters contain the inputs needed for the calibration
  m_inputPositionsInShapeParameters.clear();
  for (unsigned short int iSystem = 0; iSystem < m_systemIDs.size(); iSystem++)
  {
    std::string detector = m_detectorNames[iSystem];
    for (unsigned short int iLayer=0; iLayer<m_numLayers[iSystem]; iLayer++)
    {
      std::string decoration = Form("energy_fraction_%s_layer_%hu", detector.c_str(), iLayer);
      auto it = std::find(shapeParameters.begin(), shapeParameters.end(),
                          decoration);
      if (it != shapeParameters.end())
      {
        int position = std::distance(shapeParameters.begin(), it);
        m_inputPositionsInShapeParameters.push_back(position);
        info() << "Decoration " << decoration << " found in position " << position << " of shapeParameters" << endmsg;
      }
      else
      {
        // at least one of the inputs of the MVA was not found in the shapeParameters
        // so we can stop checking the others 
        m_inputPositionsInShapeParameters.clear();
        info() << "Decoration " << decoration << " not found. Will recalculate energy fractions from cells" << endmsg;
        break;
      }
    }
    if ((int) m_inputPositionsInShapeParameters.size() == m_numLayersTotal)
    {
      info() << "Decorations found in shapeParamers for all BDT inputs. Will not recalculate energy fractions from cells" << endmsg;
    }
  }

  info() << "Initialized the calibration" << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode CalibrateCaloClusters::execute(const EventContext& evtCtx) const
{
  (void) evtCtx;  // event context not used

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

  // Apply calibration
  {
    StatusCode sc = calibrateClusters(inClusters, outClusters);
    if (sc.isFailure())
    {
      return sc;
    }
  }

  return StatusCode::SUCCESS;
}

StatusCode CalibrateCaloClusters::finalize()
{
  if (m_ortSession)
    delete m_ortSession;
  if (m_ortEnv)
    delete m_ortEnv;

  return Gaudi::Algorithm::finalize();
}

edm4hep::ClusterCollection *CalibrateCaloClusters::initializeOutputClusters(
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

StatusCode CalibrateCaloClusters::readCalibrationFile(const std::string &calibrationFile)
{
  // set logging level based on output level of this alg
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
    m_ortEnv = new Ort::Env(loggingLevel, "ONNX runtime environment");
    Ort::SessionOptions session_options;
    session_options.SetIntraOpNumThreads(1);
    m_ortSession = new Ort::Experimental::Session(*m_ortEnv, const_cast<std::string &>(calibrationFile), session_options);
  }
  catch (const Ort::Exception &exception)
  {
    error() << "ERROR setting up ONNX runtime environment: " << exception.what() << endmsg;
    return StatusCode::FAILURE;
  }

  // print name/shape of inputs
  // use default allocator (CPU)
  Ort::AllocatorWithDefaultOptions allocator;
  debug() << "Input Node Name/Shape (" << m_input_names.size() << "):" << endmsg;
  for (std::size_t i = 0; i < m_ortSession->GetInputCount(); i++)
  {
    m_input_names.emplace_back(m_ortSession->GetInputName(i, allocator));
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
  debug() << "Output Node Name/Shape (" << m_output_names.size() << "):" << endmsg;
  for (std::size_t i = 0; i < m_ortSession->GetOutputCount(); i++)
  {
    m_output_names.emplace_back(m_ortSession->GetOutputName(i, allocator));
    m_output_shapes = m_ortSession->GetOutputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape();
    debug() << "\t" << m_output_names.at(i) << " : ";
    for (std::size_t k = 0; k < m_output_shapes.size() - 1; k++)
    {
      debug() << m_output_shapes[k] << "x";
    }
    debug() << m_output_shapes[m_output_shapes.size() - 1] << endmsg;
  }

  // the output should be a single value (the correction)
  // and the inputs should be n(layers)+1 (fractions + total E)
  // the first dimension of the tensors are the number of clusters
  // to be calibrated simultaneously (-1 = dynamic)
  // we will calibrate once at a time
  if (m_input_shapes.size() != 2 ||
      m_output_shapes.size() != 2 ||
      m_input_shapes[1] != (m_numLayersTotal + 1) ||
      m_output_shapes[1] != 1)
  {
    error() << "The input or output shapes in the calibration files do not match the expected architecture" << endmsg;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

StatusCode CalibrateCaloClusters::calibrateClusters(const edm4hep::ClusterCollection *inClusters,
                                                    edm4hep::ClusterCollection *outClusters) const
{

  // this vector will contain the input features for the calibration
  // i.e. the fraction of energy in each layer and the total energy
  std::vector<float> energiesInLayers(m_numLayersTotal + 1);

  // loop over the input clusters and perform the calibration
  for (unsigned int j = 0; j < inClusters->size(); ++j)
  {

    // retrieve total cluster energy
    float ecl = (inClusters->at(j)).getEnergy();
    if (ecl <= 0.)
    {
      warning() << "Energy in calorimeter <= 0, ignoring energy correction!" << endmsg;
      continue;
    }
    verbose() << "Cluster energy before calibration: " << ecl << endmsg;

    // calculate cluster energy in each layer and normalize by total cluster energy
    calcEnergiesInLayers(inClusters->at(j), energiesInLayers);
    verbose() << "Calibration inputs:" << endmsg;
    for (unsigned short int k = 0; k < energiesInLayers.size(); ++k)
    {
      verbose() << "    f" << k << " : " << energiesInLayers[k] << endmsg;
    }

    // run the MVA calibration and correct the cluster energy
    float corr = 1.0;
    // Create a single Ort tensor
    std::vector<Ort::Value> input_tensors;
    input_tensors.emplace_back(vec_to_tensor<float>(energiesInLayers, m_input_shapes));

    // double-check the dimensions of the input tensor
    // assert(input_tensors[0].IsTensor() && input_tensors[0].GetTensorTypeAndShapeInfo().GetShape() == m_input_shapes);

    // pass data through model
    try
    {
      std::vector<Ort::Value> output_tensors = m_ortSession->Run(m_input_names,
                                                                 input_tensors,
                                                                 m_output_names,
                                                                 Ort::RunOptions{nullptr});

      // double-check the dimensions of the output tensors
      // NOTE: the number of output tensors is equal to the number of output nodes specifed in the Run() call
      // assert(output_tensors.size() == output_names.size() && output_tensors[0].IsTensor());

      float *outputData = output_tensors[0].GetTensorMutableData<float>();
      corr = outputData[0];
    }
    catch (const Ort::Exception &exception)
    {
      error() << "ERROR running model inference: " << exception.what() << endmsg;
      return StatusCode::FAILURE;
    }

    verbose() << "Calibration output: " << corr << endmsg;
    outClusters->at(j).setEnergy(ecl * corr);
    outClusters->at(j).addToShapeParameters(ecl);
    verbose() << "Corrected cluster energy: " << ecl * corr << endmsg;
  }

  return StatusCode::SUCCESS;
}

void CalibrateCaloClusters::calcEnergiesInLayers(edm4hep::Cluster cluster,
                                                 std::vector<float> &energiesInLayers) const
{
  // reset vector with energies per layer
  std::fill(energiesInLayers.begin(), energiesInLayers.end(), 0.0);

  // get total cluster energy
  double ecl = cluster.getEnergy();

  // if all inputs are available as shapeParameters, use them
  if (m_inputPositionsInShapeParameters.size() == m_numLayersTotal)
  {
    // the shapeParameters already contain energy fractions so we
    // do not divide by cluster energy
    for (unsigned short int i=0; i < m_numLayersTotal; i++)
    {
      energiesInLayers[i] = cluster.getShapeParameters(m_inputPositionsInShapeParameters[i]);
    }
    // add as last input the total raw cluster energy
    energiesInLayers[m_numLayersTotal] = ecl;
  }
  else
  {
    // calculate the energy fractions from the cells
    // loop over the various readouts/subsystems/...
    unsigned short int startPositionToFill = 0;
    for (unsigned short int k = 0; k < m_readoutNames.size(); k++)
    {
      if (k > 0)
      {
        startPositionToFill += m_numLayers[k - 1];
      }
      int systemID = m_systemIDs[k];
      unsigned short int firstLayer = m_firstLayerIDs[k];
      std::string layerField = m_layerFieldNames[k];
      std::string readoutName = m_readoutNames[k];
      dd4hep::DDSegmentation::BitFieldCoder *decoder = m_geoSvc->getDetector()->readout(readoutName).idSpec().decoder();

      // for each readout, loop over the cells and calculate the energies in the layers
      // of that subsystem
      for (auto cell = cluster.hits_begin(); cell != cluster.hits_end(); ++cell)
      {
        dd4hep::DDSegmentation::CellID cellID = cell->getCellID();
        if (decoder->get(cellID, "system") != systemID)
        {
          continue;
        }
        int layer = decoder->get(cellID, layerField);
        energiesInLayers[startPositionToFill + layer - firstLayer] += cell->getEnergy();
      }
    }
    // divide by the cluster energy to prepare the inputs for the MVA
    for (unsigned short int k = 0; k < m_numLayersTotal; ++k)
    {
      energiesInLayers[k] /= ecl;
    }
    // add as last input the total cluster energy
    energiesInLayers[m_numLayersTotal] = ecl;
  }
}
