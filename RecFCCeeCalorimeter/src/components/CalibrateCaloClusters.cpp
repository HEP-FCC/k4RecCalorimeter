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
    : GaudiAlgorithm(name, svcLoc),
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
    StatusCode sc = GaudiAlgorithm::initialize();
    if (sc.isFailure())
    {
      return sc;
    }
  }

  // Check if readouts exist
  {
    bool readoutMissing = false;
    for (size_t i = 0; i < m_readoutNames.size(); ++i)
    {
      auto readouts = m_geoSvc->getDetector()->readouts();
      if (readouts.find(m_readoutNames.value().at(i)) == readouts.end())
      {
        readoutMissing = true;
      }
    }
    if (readoutMissing)
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
  for (size_t i = 0; i < m_numLayers.size(); ++i)
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
  info() << "Initialized the calibration" << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode CalibrateCaloClusters::execute()
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
  if (ortSession)
    delete ortSession;
  if (ortEnv)
    delete ortEnv;

  return GaudiAlgorithm::finalize();
}

edm4hep::ClusterCollection *CalibrateCaloClusters::initializeOutputClusters(
    const edm4hep::ClusterCollection *inClusters)
{
  edm4hep::ClusterCollection *outClusters = m_outClusters.createAndPut();

  for (auto const &inCluster : *inClusters)
  {
    auto outCluster = inCluster.clone();
    // verbose() << "Cluster energy before calibration:" << outCluster.getEnergy() << endmsg;
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
  case MSG::Level::FATAL: // 6
    loggingLevel = ORT_LOGGING_LEVEL_FATAL; // 4
    break;
  case MSG::Level::ERROR: // 5
    loggingLevel = ORT_LOGGING_LEVEL_ERROR; // 3
    break;
  case MSG::Level::WARNING: // 4
    loggingLevel = ORT_LOGGING_LEVEL_WARNING; // 2
    break;
  case MSG::Level::INFO: // 3
    loggingLevel = ORT_LOGGING_LEVEL_WARNING; // 2 (ORT_LOGGING_LEVEL_INFO too verbose..)
    break;
  case MSG::Level::DEBUG: // 2
    loggingLevel = ORT_LOGGING_LEVEL_INFO; // 1
    break;
  case MSG::Level::VERBOSE: // 1
    loggingLevel = ORT_LOGGING_LEVEL_VERBOSE; // 0
    break;
  }
  try
  {
    ortEnv = new Ort::Env(loggingLevel, "ONNX runtime environment");
    Ort::SessionOptions session_options;
    session_options.SetIntraOpNumThreads(1);
    ortSession = new Ort::Experimental::Session(*ortEnv, const_cast<std::string &>(calibrationFile), session_options);
  }
  catch (const Ort::Exception &exception)
  {
    error() << "ERROR setting up ONNX runtime environment: " << exception.what() << endmsg;
    return StatusCode::FAILURE;
  }

  // print name/shape of inputs
  // use default allocator (CPU)
  Ort::AllocatorWithDefaultOptions allocator;
  debug() << "Input Node Name/Shape (" << input_names.size() << "):" << endmsg;
  for (std::size_t i = 0; i < ortSession->GetInputCount(); i++)
  {
    // input_names.emplace_back(ortSession->GetInputNameAllocated(i, allocator).get());
    input_names.emplace_back(ortSession->GetInputName(i, allocator));
    input_shapes = ortSession->GetInputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape();
    debug() << "\t" << input_names.at(i) << " : ";
    for (std::size_t k = 0; k < input_shapes.size() - 1; k++)
    {
      debug() << input_shapes[k] << "x";
    }
    debug() << input_shapes[input_shapes.size() - 1] << endmsg;
  }
  // some models might have negative shape values to indicate dynamic shape, e.g., for variable batch size.
  for (auto &s : input_shapes)
  {
    if (s < 0)
    {
      s = 1;
    }
  }

  // print name/shape of outputs
  // std::vector<std::string> output_names;
  debug() << "Output Node Name/Shape (" << output_names.size() << "):" << endmsg;
  for (std::size_t i = 0; i < ortSession->GetOutputCount(); i++)
  {
    // output_names.emplace_back(ortSession->GetOutputNameAllocated(i, allocator).get());
    output_names.emplace_back(ortSession->GetOutputName(i, allocator));
    output_shapes = ortSession->GetOutputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape();
    debug() << "\t" << output_names.at(i) << " : ";
    for (std::size_t k = 0; k < output_shapes.size() - 1; k++)
    {
      debug() << output_shapes[k] << "x";
    }
    debug() << output_shapes[output_shapes.size() - 1] << endmsg;
  }

  // the output should be a single value (the correction)
  // and the inputs should be n(layers)+1 (fractions + total E)
  // the first dimension of the tensors are the number of clusters
  // to be calibrated simultaneously (-1 = dynamic)
  // we will calibrate once at a time
  if (input_shapes.size() != 2 ||
      output_shapes.size() != 2 ||
      input_shapes[1] != (m_numLayersTotal + 1) ||
      output_shapes[1] != 1)
  {
    error() << "The input or output shapes in the calibration files do not match the expected architecture" << endmsg;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

StatusCode CalibrateCaloClusters::calibrateClusters(const edm4hep::ClusterCollection *inClusters,
                                                    edm4hep::ClusterCollection *outClusters)
{

  // this vector will contain the input features for the calibration
  // i.e. the fraction of energy in each layer and the total energy
  std::vector<float> energiesInLayers(m_numLayersTotal + 1);

  // loop over the input clusters and perform the calibration
  for (size_t j = 0; j < inClusters->size(); ++j)
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
    for (size_t k = 0; k < m_numLayersTotal; ++k)
    {
      energiesInLayers[k] /= ecl;
    }
    energiesInLayers[m_numLayersTotal] = ecl;
    verbose() << "Calibration inputs:" << endmsg;
    for (size_t k = 0; k < energiesInLayers.size(); ++k)
    {
      verbose() << "    f" << k << " : " << energiesInLayers[k] << endmsg;
    }

    // run the MVA calibration and correct the cluster energy
    float corr = 1.0;
    // Create a single Ort tensor
    std::vector<Ort::Value> input_tensors;
    input_tensors.emplace_back(vec_to_tensor<float>(energiesInLayers, input_shapes));

    // double-check the dimensions of the input tensor
    // assert(input_tensors[0].IsTensor() && input_tensors[0].GetTensorTypeAndShapeInfo().GetShape() == input_shapes);

    // pass data through model
    try
    {
      std::vector<Ort::Value> output_tensors = ortSession->Run(input_names,
                                                               input_tensors,
                                                               output_names,
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
    verbose() << "Corrected cluster energy: " << ecl * corr << endmsg;
  }

  return StatusCode::SUCCESS;
}

void CalibrateCaloClusters::calcEnergiesInLayers(edm4hep::Cluster cluster,
                                                 std::vector<float> &energiesInLayer)
{
  // reset vector with energies per layer
  std::fill(energiesInLayer.begin(), energiesInLayer.end(), 0.0);

  // loop over the various readouts/subsystems/...
  int startPositionToFill = 0;
  for (size_t k = 0; k < m_readoutNames.size(); k++)
  {
    if (k > 0)
    {
      startPositionToFill += m_numLayers[k - 1];
    }
    int systemID = m_systemIDs[k];
    int firstLayer = m_firstLayerIDs[k];
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
      energiesInLayer[startPositionToFill + layer - firstLayer] += cell->getEnergy();
    }
  }
}
