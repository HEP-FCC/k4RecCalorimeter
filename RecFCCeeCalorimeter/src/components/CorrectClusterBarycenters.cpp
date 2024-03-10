#include "CorrectClusterBarycenters.h"

// k4FWCore
#include "k4Interface/IGeoSvc.h"

// edm4hep
#include "edm4hep/ClusterCollection.h"

// DD4hep
#include "DD4hep/Detector.h"
// #include "DD4hep/Readout.h"
#include "DDSegmentation/MultiSegmentation.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"
#include "detectorSegmentations/FCCSWGridModuleThetaMerged_k4geo.h"

// ROOT
#include "TLorentzVector.h"

DECLARE_COMPONENT(CorrectClusterBarycenters)

CorrectClusterBarycenters::CorrectClusterBarycenters(const std::string& name, ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc),
      m_geoSvc("GeoSvc", "CorrectClusterBarycenters")
{
  declareProperty("clusters", m_inClusters, "Input clusters (input)");
  declareProperty("correctedClusters", m_correctedClusters, "Corrected clusters (output)");
}

StatusCode CorrectClusterBarycenters::initialize()
{
  {
    StatusCode sc = GaudiAlgorithm::initialize();
    if (sc.isFailure())
      return sc;
  }

  // check if readouts exist
  if (m_geoSvc->getDetector()->readouts().find(m_readoutName) == m_geoSvc->getDetector()->readouts().end()) {
    error() << "Readout <<" << m_readoutName << ">> does not exist." << endmsg;
    return StatusCode::FAILURE;
  }
  // retrieve PhiTheta segmentation
  m_segmentationPhiTheta = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo*>(
      m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());
  // m_segmentationMulti = dynamic_cast<dd4hep::DDSegmentation::MultiSegmentation*>(
  //     m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());
  if (m_segmentationPhiTheta == nullptr
      // && m_segmentationMulti == nullptr
      ) {
    error() << "There is no phi-theta segmentation." << endmsg;
    return StatusCode::FAILURE;
  }

  m_decoder = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder();
  return StatusCode::SUCCESS;

}

StatusCode CorrectClusterBarycenters::finalize()
{
  return GaudiAlgorithm::finalize();
}

StatusCode CorrectClusterBarycenters::execute()
{
  // Get the input collection with clusters
  const edm4hep::ClusterCollection *inClusters = m_inClusters.get();
  edm4hep::ClusterCollection *correctedClusters = m_correctedClusters.createAndPut();

  uint systemId = m_systemId;
  const dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo *segmentation = nullptr;
  if (m_segmentationPhiTheta != nullptr)
  {
    segmentation = m_segmentationPhiTheta;
  }
  for (const auto &cluster : *inClusters)
  {
    double oldEnergy = 0;
    TVector3 pos(cluster.getPosition().x, cluster.getPosition().y, cluster.getPosition().z);
    double oldTheta = pos.Theta();
    double oldPhi = pos.Phi();
    for (auto cell = cluster.hits_begin(); cell != cluster.hits_end(); cell++)
    {
      oldEnergy += cell->getEnergy();
    }
    verbose() << " OLD ENERGY = " << oldEnergy << " from " << cluster.hits_size() << " cells" << endmsg;
    verbose() << " OLD CLUSTER ENERGY = " << cluster.getEnergy() << endmsg;

    // Do everything only using the first defined calorimeter (default: Ecal barrel)
    double oldThetaId = -1;
    double oldPhiId = -1;
    if (m_segmentationPhiTheta != nullptr)
    {
      oldThetaId = int(floor((oldTheta + 0.5 * segmentation->gridSizeTheta() - segmentation->offsetTheta()) / segmentation->gridSizeTheta()));
      oldPhiId = int(floor((oldPhi + 0.5 * segmentation->gridSizePhi() - segmentation->offsetPhi()) / segmentation->gridSizePhi()));
    }
    // 0. Create new cluster, copy information from input
    auto newCluster = correctedClusters->create();
    double energy = 0;
    newCluster.setPosition(cluster.getPosition());
    for (auto cell = cluster.hits_begin(); cell != cluster.hits_end(); cell++)
    {
      // if (m_segmentationMulti != nullptr)
      // {
      //   segmentation = dynamic_cast<const dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo *>(&m_segmentationMulti->subsegmentation(cell->getCellID()));
      //   oldThetaId = int(floor((oldTheta + 0.5 * segmentation->gridSizeTheta() - segmentation->offsetTheta()) / segmentation->gridSizeTheta()));
      //   oldPhiId = int(floor((oldPhi + 0.5 * segmentation->gridSizePhi() - segmentation->offsetPhi()) / segmentation->gridSizePhi()));
      // }
      if (m_decoder->get(cell->getCellID(), "system") == systemId)
      {
        uint layerId = m_decoder->get(cell->getCellID(), "layer");
        if (m_nPhiFinal[layerId] > 0 && m_nThetaFinal[layerId] > 0)
        {
          uint thetaId = m_decoder->get(cell->getCellID(), "theta");
          uint phiId = m_decoder->get(cell->getCellID(), "phi");
          if (thetaId >= (oldThetaId - m_halfThetaFin[layerId]) && thetaId <= (oldThetaId + m_halfThetaFin[layerId]) &&
              phiId >= phiNeighbour((oldPhiId - m_halfPhiFin[layerId]), segmentation->phiBins()) && phiId <= phiNeighbour((oldPhiId + m_halfPhiFin[layerId]), segmentation->phiBins()))
          {

            newCluster.addToHits(*cell);
            energy += cell->getEnergy();

          }
        }
      }
    }
    newCluster.setEnergy(energy);

    // 1. Correct theta position with log-weighting

    // get current pseudorapidity
    std::vector<double> sumEnLayer;
    std::vector<double> sumThetaLayer;
    std::vector<double> sumWeightLayer;
    sumEnLayer.assign(m_numLayers, 0);
    sumThetaLayer.assign(m_numLayers, 0);
    sumWeightLayer.assign(m_numLayers, 0);
    // first check the energy deposited in each layer
    for (auto cell = newCluster.hits_begin(); cell != newCluster.hits_end(); cell++)
    {
      dd4hep::DDSegmentation::CellID cID = cell->getCellID();
      uint layer = m_decoder->get(cID, m_layerFieldName) + m_firstLayerId;
      sumEnLayer[layer] += cell->getEnergy();
    }
    // repeat but calculating theta barycentre in each layer
    for (auto cell = newCluster.hits_begin(); cell != newCluster.hits_end(); cell++)
    {
      // if (m_segmentationMulti != nullptr)
      // {
      //   segmentation = dynamic_cast<const dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo *>(&m_segmentationMulti->subsegmentation(cell->getCellID()));
      // }
      dd4hep::DDSegmentation::CellID cID = cell->getCellID();
      uint layer = m_decoder->get(cID, m_layerFieldName) + m_firstLayerId;
      double weightLog = std::max(0., m_thetaRecalcLayerWeights[layer] + log(cell->getEnergy() / sumEnLayer[layer]));
      double theta = segmentation->theta(cell->getCellID());
      sumThetaLayer[layer] += (weightLog * theta);
      sumWeightLayer[layer] += weightLog;
    }
    // calculate theta position weighting with energy deposited in layer
    // this energy is a good estimator of 1/sigma^2 of (theta_barycentre-theta_MC) distribution
    // double layerWeight = 0;
    // double sumLayerWeight = 0;
    // double sumLayerWeight2point = 0;
    // double newTheta = 0;
    // double newThetaErrorRes = 0;
    // double newThetaErrorRes2point = 0;
    // for (uint iLayer = 0; iLayer < m_numLayers; iLayer++)
    // {
    //   if (sumWeightLayer[iLayer] > 1e-10)
    //   {
    //     sumThetaLayer[iLayer] /= sumWeightLayer[iLayer];
    //     newTheta += sumThetaLayer[iLayer] * sumEnLayer[iLayer];
    //     // layerWeight = 1. / (pow(m_thetaLayerResolutionSampling[iLayer], 2) / energy +  pow(m_thetaLayerResolutionConst[iLayer], 2));
    //     // sumLayerWeight += layerWeight;
    //     // newThetaErrorRes += sumThetaLayer[iLayer] * layerWeight;
    //     // if (iLayer == 1 || iLayer == 2) {
    //     // newThetaErrorRes2point += sumThetaLayer[iLayer] * layerWeight;
    //     // sumLayerWeight2point += layerWeight;
    //     // }
    //   }
    // }
    // newTheta /= energy;
    // // newThetaErrorRes /= sumLayerWeight;
    // // newThetaErrorRes2point /= sumLayerWeight2point;
    // // alter Cartesian position of a cluster using new theta position
    // double radius = pos.Perp();
    // double phi = pos.Phi();
    // auto newClusterPosition = edm4hep::Vector3f(radius * cos(phi), radius * sin(phi), radius * sinh(newTheta));
    // // newCluster.setPosition(newClusterPosition);
    // // NOT SAVING NEW POSITION FOR NOW
    newCluster.setPosition(cluster.getPosition());
  }
  return StatusCode::SUCCESS;
}

unsigned int CorrectClusterBarycenters::phiNeighbour(int aIPhi, int aMaxPhi) const
{
  if (aIPhi < 0)
  {
    return aMaxPhi + aIPhi;
  }
  else if (aIPhi >= aMaxPhi)
  {
    return aIPhi % aMaxPhi;
  }
  return aIPhi;
}
