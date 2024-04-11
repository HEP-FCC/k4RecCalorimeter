#include "AugmentClustersFCCee.h"

// k4FWCore
#include "k4Interface/IGeoSvc.h"

// edm4hep
#include "edm4hep/ClusterCollection.h"

// DD4hep
#include "DD4hep/Detector.h"

// ROOT
#include "TLorentzVector.h"
#include "TString.h"
#include "TVector3.h"
#include "TMath.h"

DECLARE_COMPONENT(AugmentClustersFCCee)

AugmentClustersFCCee::AugmentClustersFCCee(const std::string &name, ISvcLocator *svcLoc)
    : Gaudi::Algorithm(name, svcLoc),
      m_geoSvc("GeoSvc", "AugmentClustersFCCee")
{
  declareProperty("inClusters", m_inClusters, "Input clusters");
  declareProperty("outClusters", m_outClusters, "Output clusters");
}

StatusCode AugmentClustersFCCee::initialize()
{
  {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure())
      return sc;
  }

  // check if readouts exist
  for (size_t k = 0; k < m_readoutNames.size(); k++)
  {
    std::string readoutName = m_readoutNames[k];
    if (m_geoSvc->getDetector()->readouts().find(readoutName) == m_geoSvc->getDetector()->readouts().end())
    {
      error() << "Readout <<" << readoutName << ">> does not exist." << endmsg;
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
  if (m_systemIDs.size() != m_detectorNames.size())
  {
    error() << "Sizes of systemIDs vector and detectorNames vector do not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }

  // initialise the list of metadata for the clusters
  // append to the metadata of the input clusters
  std::vector<std::string> showerShapeDecorations = m_inShapeParameterHandle.get();
  for (size_t k = 0; k < m_detectorNames.size(); k++)
  {
    const char *detector = m_detectorNames[k].c_str();
    for (unsigned layer = 0; layer < m_numLayers[k]; layer++)
    {
      showerShapeDecorations.push_back(Form("f%s%d", detector, layer));
      showerShapeDecorations.push_back(Form("theta%s%d", detector, layer));
      showerShapeDecorations.push_back(Form("phi%s%d", detector, layer));
    }
  }
  m_showerShapeHandle.put(showerShapeDecorations);

  return StatusCode::SUCCESS;
}

StatusCode AugmentClustersFCCee::finalize()
{
  return Gaudi::Algorithm::finalize();
}

StatusCode AugmentClustersFCCee::execute(const EventContext &evtCtx) const
{
  (void)evtCtx; // event context not used

  // get the input collection with clusters
  const edm4hep::ClusterCollection *inClusters = m_inClusters.get();

  // create the new output collection
  edm4hep::ClusterCollection *outClusters = m_outClusters.createAndPut();

  // total number of layers
  size_t numLayersTotal = 0;
  for (size_t i = 0; i < m_numLayers.size(); ++i)
  {
    numLayersTotal += m_numLayers[i];
  }
  
  // loop over the clusters, clone them, and calculate the shape parameters to store with them
  for (const auto &cluster : *inClusters)
  {
    // clone original cluster
    auto newCluster = cluster.clone();
    outClusters->push_back(newCluster);

    // calculate the energy deposited in each layer
    // also, find out if cluster is around -pi..pi transition
    double E(0.0);
    std::vector<double> sumEnLayer;
    sumEnLayer.assign(numLayersTotal, 0.);
    double phiMin = 9999.;
    double phiMax = -9999.;

    // loop over each system/readout
    int startPositionToFill = 0;
    for (size_t k = 0; k < m_readoutNames.size(); k++)
    {
      if (k > 0)
      {
        startPositionToFill += m_numLayers[k - 1];
      }
      int systemID = m_systemIDs[k];
      std::string layerField = m_layerFieldNames[k];
      std::string readoutName = m_readoutNames[k];
      dd4hep::DDSegmentation::BitFieldCoder *decoder = m_geoSvc->getDetector()->readout(readoutName).idSpec().decoder();

      // loop over the cells
      for (auto cell = newCluster.hits_begin(); cell != newCluster.hits_end(); cell++)
      {
        dd4hep::DDSegmentation::CellID cID = cell->getCellID();
        int sysId = decoder->get(cID, "system");
        if (sysId != systemID)
          continue;
        uint layer = decoder->get(cID, layerField);
        sumEnLayer[layer+startPositionToFill] += cell->getEnergy();
        E += cell->getEnergy();
        TVector3 v = TVector3(cell->getPosition().x, cell->getPosition().y, cell->getPosition().z);
        double phi = v.Phi();
        if (phi < phiMin)
          phiMin = phi;
        if (phi > phiMax)
          phiMax = phi;
      }
    }
  
    // any number close to two pi should do, because if a cluster contains
    // the -pi<->pi transition, phiMin should be close to -pi and phiMax close to pi
    bool isClusterPhiNearPi = false;
    if (phiMax - phiMin > 6.)
      isClusterPhiNearPi = true;
    debug() << "phiMin, phiMax : " << phiMin << " " << phiMax << endmsg;
    debug() << "Cluster is near phi=pi : " << isClusterPhiNearPi << endmsg;

    // calculate the theta positions with log(E) weighting in each layer
    // for phi use standard E weighting
    std::vector<double> sumThetaLayer;
    std::vector<double> sumPhiLayer;
    std::vector<double> sumWeightLayer;
    sumThetaLayer.assign(numLayersTotal, 0);
    sumPhiLayer.assign(numLayersTotal, 0);
    sumWeightLayer.assign(numLayersTotal, 0);
    // loop over each system/readout
    startPositionToFill = 0;
    // rather than a loop over the systems, could first determine systemID
    // from cellID, match it against the systemIDs and skip cell if match not found
    // can do it in a separate PR
    for (size_t k = 0; k < m_readoutNames.size(); k++)
    {
      if (k > 0)
      {
        startPositionToFill += m_numLayers[k - 1];
      }
      int systemID = m_systemIDs[k];
      std::string layerField = m_layerFieldNames[k];
      std::string readoutName = m_readoutNames[k];
      dd4hep::DDSegmentation::BitFieldCoder *decoder = m_geoSvc->getDetector()->readout(readoutName).idSpec().decoder();

      for (auto cell = newCluster.hits_begin(); cell != newCluster.hits_end(); cell++)
      {
        dd4hep::DDSegmentation::CellID cID = cell->getCellID();
        int sysId = decoder->get(cID, "system");
        if (sysId != systemID)
          continue;
        uint layer = decoder->get(cID, layerField);
        double eCell = cell->getEnergy();
        double weightLog = std::max(0., m_thetaRecalcLayerWeights[k][layer] + log(eCell / sumEnLayer[layer]));
        TVector3 v = TVector3(cell->getPosition().x, cell->getPosition().y, cell->getPosition().z);
        double theta = v.Theta();
        double phi = v.Phi();
        // for clusters that are around the -pi<->pi transition, we want to avoid averaging
        // over phi values that might differ by 2pi. in that case, for cells with negative
        // phi we add two pi, so that we average phi values all close to pi
        if (isClusterPhiNearPi && phi < 0.)
          phi += TMath::TwoPi();
        if (m_thetaRecalcLayerWeights[k][layer]<0)
          sumThetaLayer[layer+startPositionToFill] += (eCell * theta);
        else
          sumThetaLayer[layer+startPositionToFill] += (weightLog * theta);
        sumWeightLayer[layer+startPositionToFill] += weightLog;
        sumPhiLayer[layer+startPositionToFill] += (eCell * phi);
      }
    }

    // save energy and theta/phi positions per layer in shape parameters
    startPositionToFill = 0;
    for (size_t k = 0; k < m_readoutNames.size(); k++)
    {
      if (k > 0)
      {
        startPositionToFill += m_numLayers[k - 1];
      }
      for (unsigned layer = 0; layer < m_numLayers[k]; layer++)
      {
        if (m_thetaRecalcLayerWeights[k][layer]<0)
        {
          if (sumEnLayer[layer+startPositionToFill] != 0.0)
          {
            sumEnLayer[layer+startPositionToFill] /= sumWeightLayer[layer+startPositionToFill];
          }
        }
        else
        {
          if (sumWeightLayer[layer+startPositionToFill] != 0.0)
          {
            sumThetaLayer[layer+startPositionToFill] /= sumWeightLayer[layer+startPositionToFill];
          }
        }
        if (sumEnLayer[layer+startPositionToFill] != 0.0)
        {
          sumPhiLayer[layer+startPositionToFill] /= sumEnLayer[layer+startPositionToFill];
        }
        // make sure phi is in range -pi..pi
        if (sumPhiLayer[layer+startPositionToFill] > TMath::Pi())
          sumPhiLayer[layer+startPositionToFill] -= TMath::TwoPi();
        newCluster.addToShapeParameters(sumEnLayer[layer+startPositionToFill] / E);
        newCluster.addToShapeParameters(sumThetaLayer[layer+startPositionToFill]);
        newCluster.addToShapeParameters(sumPhiLayer[layer+startPositionToFill]);
      }
    }
  }
  return StatusCode::SUCCESS;
}
