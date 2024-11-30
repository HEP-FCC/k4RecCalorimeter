#include "CaloTowerToolFCCee.h"

// k4FWCore
#include "k4Interface/IGeoSvc.h"

// edm4hep
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/Cluster.h"
#include "edm4hep/MutableCluster.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DD4hep/Readout.h"
#include "DDSegmentation/MultiSegmentation.h"

// k4geo
#include "detectorSegmentations/FCCSWGridModuleThetaMerged_k4geo.h"
#include "detectorSegmentations/FCCSWGridPhiTheta_k4geo.h"
#include "detectorSegmentations/FCCSWEndcapTurbine_k4geo.h"
#include "detectorSegmentations/FCCSWHCalPhiTheta_k4geo.h"
#include "detectorSegmentations/FCCSWHCalPhiRow_k4geo.h"


DECLARE_COMPONENT(CaloTowerToolFCCee)

CaloTowerToolFCCee::CaloTowerToolFCCee(const std::string& type, const std::string& name, const IInterface* parent)
    : AlgTool(type, name, parent), m_geoSvc("GeoSvc", name) {
  declareProperty("ecalBarrelCells", m_ecalBarrelCells, "");
  declareProperty("ecalEndcapCells", m_ecalEndcapCells, "");
  declareProperty("ecalFwdCells", m_ecalFwdCells, "");
  declareProperty("hcalBarrelCells", m_hcalBarrelCells, "");
  declareProperty("hcalExtBarrelCells", m_hcalExtBarrelCells, "");
  declareProperty("hcalEndcapCells", m_hcalEndcapCells, "");
  declareProperty("hcalFwdCells", m_hcalFwdCells, "");
  declareInterface<ITowerToolThetaModule>(this);
}

StatusCode CaloTowerToolFCCee::initialize() {
  if (AlgTool::initialize().isFailure()) {
    return StatusCode::FAILURE;
  }
 
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service. "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }

  std::vector<double> listPhiMax;
  std::vector<double> listThetaMax;
  listPhiMax.reserve(7);
  listThetaMax.reserve(7);

  std::pair<double, double> tmpPair;
  if (retrievePhiThetaExtrema(m_ecalBarrelReadoutName, tmpPair) == StatusCode::SUCCESS) {
    listPhiMax.push_back(tmpPair.first);
    listThetaMax.push_back(tmpPair.second);
  }
  if (retrievePhiThetaExtrema(m_ecalEndcapReadoutName, tmpPair) == StatusCode::SUCCESS) {
    listPhiMax.push_back(tmpPair.first);
    listThetaMax.push_back(tmpPair.second);
  }
  if (retrievePhiThetaExtrema(m_ecalFwdReadoutName, tmpPair) == StatusCode::SUCCESS) {
    listPhiMax.push_back(tmpPair.first);
    listThetaMax.push_back(tmpPair.second);
  }
  if (retrievePhiThetaExtrema(m_hcalBarrelReadoutName, tmpPair) == StatusCode::SUCCESS) {
    listPhiMax.push_back(tmpPair.first);
    listThetaMax.push_back(tmpPair.second);
  }
  if (retrievePhiThetaExtrema(m_hcalExtBarrelReadoutName, tmpPair) == StatusCode::SUCCESS) {
    listPhiMax.push_back(tmpPair.first);
    listThetaMax.push_back(tmpPair.second);
  }
  if (retrievePhiThetaExtrema(m_hcalEndcapReadoutName, tmpPair) == StatusCode::SUCCESS) {
    listPhiMax.push_back(tmpPair.first);
    listThetaMax.push_back(tmpPair.second);
  }
  if (retrievePhiThetaExtrema(m_hcalFwdReadoutName, tmpPair) == StatusCode::SUCCESS) {
    listPhiMax.push_back(tmpPair.first);
    listThetaMax.push_back(tmpPair.second);
  }
  // Maximum theta & phi of the calorimeter system
  m_phiMax = *std::max_element(listPhiMax.begin(), listPhiMax.end());
  m_thetaMax = *std::max_element(listThetaMax.begin(), listThetaMax.end());
  debug() << "Detector limits: phiMax " << m_phiMax << " thetaMax " << m_thetaMax << endmsg;

  // very small number (epsilon) substructed from the edges to ensure correct division
  float epsilon = 0.0001;
  // number of phi bins
  m_nPhiTower = ceil(2 * (m_phiMax - epsilon) / m_deltaPhiTower);
  // number of theta bins
  m_nThetaTower = ceil(2 * (fabs(m_thetaMax - M_PI/2.) - epsilon) / m_deltaThetaTower);
  debug() << "Towers: thetaMax " << m_thetaMax << ", deltaThetaTower " << m_deltaThetaTower << ", nThetaTower " << m_nThetaTower
          << endmsg;
  debug() << "Towers: phiMax " << m_phiMax << ", deltaPhiTower " << m_deltaPhiTower << ", nPhiTower " << m_nPhiTower
          << endmsg;

  return StatusCode::SUCCESS;
}

StatusCode CaloTowerToolFCCee::finalize() { 
  for (auto& towerInMap : m_cellsInTowers) {
    towerInMap.second.clear();
  }
  return AlgTool::finalize();
}

StatusCode  CaloTowerToolFCCee::retrievePhiThetaExtrema(std::string aReadoutName, std::pair<double, double> &phiThetaPair) {

  double phiMax = -1;
  double thetaMax = -1;

  
  // check if readout exists & retrieve Module-Theta segmentation
  // if readout does not exist, reconstruction without this calorimeter part will be performed

  std::pair<dd4hep::DDSegmentation::Segmentation*, SegmentationType> tmpPair;
  tmpPair = retrieveSegmentation(aReadoutName);
  dd4hep::DDSegmentation::Segmentation* aSegmentation = tmpPair.first;
  SegmentationType aType= tmpPair.second;
  if (aSegmentation != nullptr && aType == SegmentationType::kWrong) {
    error() << "Wrong type of segmentation" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_useHalfTower) {
    m_decoder = m_geoSvc->getDetector()->readout(aReadoutName).idSpec().decoder();
  }
  
  if (aSegmentation != nullptr) {
    
    switch (aType) {
    case SegmentationType::kModuleTheta: {
      
      dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo* segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo*>(aSegmentation);
      phiMax = M_PI - M_PI/segmentation->nModules();
      thetaMax = M_PI - fabs(segmentation->offsetTheta()) + segmentation->gridSizeTheta() * .5;
      break;
    }
    case SegmentationType::kPhiTheta: {
      info() << "== Retrieving segmentation " << aSegmentation->name() << endmsg;
      dd4hep::DDSegmentation::FCCSWHCalPhiTheta_k4geo* segmentationHCal = dynamic_cast<dd4hep::DDSegmentation::FCCSWHCalPhiTheta_k4geo*>(aSegmentation);
      phiMax = M_PI - M_PI/segmentationHCal->phiBins();
      thetaMax = M_PI - fabs(segmentationHCal->offsetTheta()) + segmentationHCal->gridSizeTheta() * .5;
      break;
    }
    case SegmentationType::kPhiRow: {
      info() << "== Retrieving segmentation " << aSegmentation->name() << endmsg;
      dd4hep::DDSegmentation::FCCSWHCalPhiRow_k4geo* segmentationHCal = dynamic_cast<dd4hep::DDSegmentation::FCCSWHCalPhiRow_k4geo*>(aSegmentation);
      phiMax = M_PI - M_PI/segmentationHCal->phiBins();
      thetaMax = segmentationHCal->thetaMax();
      break;
    }
    case SegmentationType::kMulti: {
      double phi = -1;
      double theta = -1;
      info() << "== Retrieving multi segmentation " << aSegmentation->name() << endmsg;
      for (const auto& subSegm: dynamic_cast<dd4hep::DDSegmentation::MultiSegmentation*>(aSegmentation)->subSegmentations()) {
	dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo* segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo*>(subSegm.segmentation);
        phi = M_PI - M_PI/segmentation->nModules();
	theta = M_PI - fabs(segmentation->offsetTheta()) + segmentation->gridSizeTheta() * .5;
        if (theta > thetaMax) { thetaMax = theta;}
        if (phi > phiMax) { phiMax = phi;}
      }
      break;
    }
    case SegmentationType::kEndcapTurbine: {
      dd4hep::DDSegmentation::FCCSWEndcapTurbine_k4geo* ECTsegmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWEndcapTurbine_k4geo*>(aSegmentation);
      thetaMax = M_PI - fabs(ECTsegmentation->offsetTheta());
      phiMax = M_PI;
      break;
    } 
    case SegmentationType::kWrong: {
      info() << "== Retrieving WRONG segmentation" << endmsg;
      phiMax = -1;
      thetaMax = -1;
      return StatusCode::FAILURE;
    }
    default: {
      error() << " Unsupported segmentation" << endmsg;
      phiMax = -1;
      thetaMax = -1;
      return StatusCode::FAILURE;
    }
    }
  }
  phiThetaPair = std::make_pair(phiMax, thetaMax);
  return StatusCode::SUCCESS;
}

void CaloTowerToolFCCee::towersNumber(int& nTheta, int& nPhi) {

  nTheta = m_nThetaTower;
  nPhi = m_nPhiTower;
}

uint CaloTowerToolFCCee::buildTowers(std::vector<std::vector<float>>& aTowers, bool fillTowersCells) {
  uint totalNumberOfCells = 0;
  for (auto& towerInMap : m_cellsInTowers) {
    towerInMap.second.clear();
  }
  // 1. ECAL barrel
  // Get the input collection with calorimeter cells
  const edm4hep::CalorimeterHitCollection* ecalBarrelCells = m_ecalBarrelCells.get();
  debug() << "Input Ecal barrel cell collection size: " << ecalBarrelCells->size() << endmsg;
  // Loop over a collection of calorimeter cells and build calo towers
  if (ecalBarrelCells->size() >0) {
    CellsIntoTowers(aTowers, ecalBarrelCells, fillTowersCells);
    totalNumberOfCells += ecalBarrelCells->size();
  }

  // 2. ECAL endcap calorimeter
  const edm4hep::CalorimeterHitCollection* ecalEndcapCells = m_ecalEndcapCells.get();
  debug() << "Input Ecal endcap cell collection size: " << ecalEndcapCells->size() << endmsg;
  // Loop over a collection of calorimeter cells and build calo towers
  if (ecalEndcapCells->size() > 0) {
    CellsIntoTowers(aTowers, ecalEndcapCells, fillTowersCells);
    totalNumberOfCells += ecalEndcapCells->size();
  }
  /*
  // 3. ECAL forward calorimeter
  const edm4hep::CalorimeterHitCollection* ecalFwdCells = m_ecalFwdCells.get();
  debug() << "Input Ecal forward cell collection size: " << ecalFwdCells->size() << endmsg;
  // Loop over a collection of calorimeter cells and build calo towers
  if (ecalFwdCells->size() > 0) {
    CellsIntoTowers(aTowers, ecalFwdCells, fillTowersCells);
    totalNumberOfCells += ecalFwdCells->size();
  }
  */
  // 4. HCAL barrel
  const edm4hep::CalorimeterHitCollection* hcalBarrelCells = m_hcalBarrelCells.get();
  debug() << "Input hadronic barrel cell collection size: " << hcalBarrelCells->size() << endmsg;
  // Loop over a collection of calorimeter cells and build calo towers
  if (hcalBarrelCells->size()>0) {
    CellsIntoTowers(aTowers, hcalBarrelCells, fillTowersCells);
    totalNumberOfCells += hcalBarrelCells->size();
  }
  /*
  // 5. HCAL extended barrel
  const edm4hep::CalorimeterHitCollection* hcalExtBarrelCells = m_hcalExtBarrelCells.get();
  debug() << "Input hadronic extended barrel cell collection size: " << hcalExtBarrelCells->size() << endmsg;
  // Loop over a collection of calorimeter cells and build calo towers
  if (hcalExtBarrelCells->size() >0) {
    CellsIntoTowers(aTowers, hcalExtBarrelCells, fillTowersCells);
    totalNumberOfCells += hcalExtBarrelCells->size();
  }
  */
  // 6. HCAL endcap calorimeter
  const edm4hep::CalorimeterHitCollection* hcalEndcapCells = m_hcalEndcapCells.get();
  debug() << "Input Hcal endcap cell collection size: " << hcalEndcapCells->size() << endmsg;
  // Loop over a collection of calorimeter cells and build calo towers
  if (hcalEndcapCells->size() > 0) {
    CellsIntoTowers(aTowers, hcalEndcapCells,  fillTowersCells);
    totalNumberOfCells += hcalEndcapCells->size();
  }
  /*
  // 7. HCAL forward calorimeter
  const edm4hep::CalorimeterHitCollection* hcalFwdCells = m_hcalFwdCells.get();
  debug() << "Input Hcal forward cell collection size: " << hcalFwdCells->size() << endmsg;
  // Loop over a collection of calorimeter cells and build calo towers
  if (hcalFwdCells->size()>0) {
    CellsIntoTowers(aTowers, hcalFwdCells, fillTowersCells);
    totalNumberOfCells += hcalFwdCells->size();
  }
*/
  return totalNumberOfCells;
}

// Get the tower IDs in theta
// aTheta Position of the calorimeter cell in theta
uint CaloTowerToolFCCee::idTheta(float aTheta) const {
  uint id = floor(( m_thetaMax - aTheta ) / m_deltaThetaTower);
  return id;
}

// Get the tower IDs in phi
// aPhi Position of the calorimeter cell in phi
uint CaloTowerToolFCCee::idPhi(float aPhi) const {
  uint id = floor((aPhi + m_phiMax) / m_deltaPhiTower);
  return id;
}

// Get the theta position of the centre of the tower
float CaloTowerToolFCCee::theta(int aIdTheta) const {
  // middle of the tower
  return (m_thetaMax - (aIdTheta + 0.5) * m_deltaThetaTower);
}

// Get the phi position of the centre of the tower
float CaloTowerToolFCCee::phi(int aIdPhi) const {
  // middle of the tower
  return ((aIdPhi + 0.5) * m_deltaPhiTower - m_phiMax);
}

uint CaloTowerToolFCCee::phiNeighbour(int aIPhi) const {
  if (aIPhi < 0) {
    return m_nPhiTower + aIPhi;
  } else if (aIPhi >= m_nPhiTower) {
    return aIPhi % m_nPhiTower;
  }
  return aIPhi;
}

std::map<std::pair<uint, uint>, std::vector<edm4hep::CalorimeterHit>> CaloTowerToolFCCee::cellsInTowers() const { return m_cellsInTowers; }

// to fill the cell infomation into towers
void CaloTowerToolFCCee::CellsIntoTowers(std::vector<std::vector<float>>& aTowers,
                                    const edm4hep::CalorimeterHitCollection* aCells,
                                    bool fillTowersCells) {
  // Loop over a collection of calorimeter cells and build calo towers
  // tower index of the borders of the cell
  int iTheta = 0;
  int iPhi = 0;

  bool pass = true;
  for (const auto& cell : *aCells) {
    pass = true;
    if (m_useHalfTower) {
        uint layerId = m_decoder->get(cell.getCellID(), "layer");
        if ( layerId > m_max_layer ) {
          pass = false;
        }
    }
    if (pass) {
      float cellX = cell.getPosition().x;
      float cellY = cell.getPosition().y;
      float cellZ = cell.getPosition().z;
      float cellTheta = atan2(sqrt(cellX * cellX + cellY * cellY), cellZ);
      float cellPhi   = atan2(cellY, cellX);
      iTheta = idTheta(cellTheta);
      iPhi = idPhi(cellPhi);
      aTowers[iTheta][phiNeighbour(iPhi)] +=
	cell.getEnergy() * sin(cellTheta);

      if (fillTowersCells) {
        m_cellsInTowers[std::make_pair(iTheta, phiNeighbour(iPhi))].push_back(cell.clone());
	int ncells = m_cellsInTowers[std::make_pair(iTheta, phiNeighbour(iPhi))].size();

	if ( ncells > 5)
	  verbose() << "NUM CELLs IN TOWER : " << ncells << endmsg;
      }
    }
  }
}

std::pair<dd4hep::DDSegmentation::Segmentation*, CaloTowerToolFCCee::SegmentationType> CaloTowerToolFCCee::retrieveSegmentation(std::string aReadoutName) {
  dd4hep::DDSegmentation::Segmentation* segmentation = nullptr;
  if (m_geoSvc->getDetector()->readouts().find(aReadoutName) == m_geoSvc->getDetector()->readouts().end()) {
    info() << "Readout does not exist! Please check if it is correct. Processing without it." << endmsg;
  }
  else
  {
    info() << "Readout " << aReadoutName << " found." << endmsg;
    segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo*>(
      m_geoSvc->getDetector()->readout(aReadoutName).segmentation().segmentation());
    if (segmentation == nullptr)
    {
      segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWHCalPhiTheta_k4geo*>(m_geoSvc->getDetector()->readout(aReadoutName).segmentation().segmentation());
      if (segmentation == nullptr)
      {
        segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWHCalPhiRow_k4geo*>(m_geoSvc->getDetector()->readout(aReadoutName).segmentation().segmentation());
        if (segmentation == nullptr)
        {
          segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWEndcapTurbine_k4geo*>(m_geoSvc->getDetector()->readout(aReadoutName).segmentation().segmentation());
          if (segmentation == nullptr)
          {
	    segmentation = dynamic_cast<dd4hep::DDSegmentation::MultiSegmentation*>( m_geoSvc->getDetector()->readout(aReadoutName).segmentation().segmentation());
	    if (segmentation == nullptr)
            {
	      warning() << "There is no module-theta, phi-theta, phi-row, endcap turbine, or multi- segmentation for the readout " << aReadoutName << " defined." << endmsg;
	      return std::make_pair(nullptr, SegmentationType::kWrong);
	    }
            else
            {
	      // check if multisegmentation contains only module-theta sub-segmentations
	      dd4hep::DDSegmentation::Segmentation* subsegmentation = nullptr;
	      for (const auto& subSegm: dynamic_cast<dd4hep::DDSegmentation::MultiSegmentation*>(segmentation)->subSegmentations())
              {
	        subsegmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo*>(subSegm.segmentation);
	        if (subsegmentation == nullptr)
                {
                  warning() << "At least one of the sub-segmentations in MultiSegmentation named " << aReadoutName << " is not a module-theta grid." << endmsg;
	          return std::make_pair(nullptr, SegmentationType::kWrong);
	        }
	      }
	      return std::make_pair(segmentation, SegmentationType::kMulti);
	    }
          } else return std::make_pair(segmentation, SegmentationType::kEndcapTurbine);
        } else return std::make_pair(segmentation, SegmentationType::kPhiRow);
      } else return std::make_pair(segmentation, SegmentationType::kPhiTheta);
    } else return std::make_pair(segmentation, SegmentationType::kModuleTheta);
  }
  return std::make_pair(segmentation, SegmentationType::kWrong);
}

void CaloTowerToolFCCee::attachCells(float theta, float phi, uint halfThetaFin, uint halfPhiFin, edm4hep::MutableCluster& aEdmCluster, edm4hep::CalorimeterHitCollection* aEdmClusterCells, bool aEllipse) {
  int thetaId = idTheta(theta);
  int phiId = idPhi(phi);
  std::vector<dd4hep::DDSegmentation::CellID> seen_cellIDs;
  if (aEllipse) {
    for (int iTheta = thetaId - halfThetaFin; iTheta <= int(thetaId + halfThetaFin); iTheta++) {
      for (int iPhi = phiId - halfPhiFin; iPhi <= int(phiId + halfPhiFin); iPhi++) {
        if (pow( (thetaId - iTheta) / (halfThetaFin + 0.5), 2) + pow( (phiId - iPhi) / (halfPhiFin + 0.5), 2) < 1) {
          for (auto cell : m_cellsInTowers[std::make_pair(iTheta, phiNeighbour(iPhi))]) {
            if (std::find(seen_cellIDs.begin(), seen_cellIDs.end(), cell.getCellID()) != seen_cellIDs.end()) { // towers can be smaller than cells in which case a cell belongs to several towers
               continue;
            }
            seen_cellIDs.push_back(cell.getCellID());
            auto cellclone = cell.clone();
            aEdmClusterCells->push_back(cellclone);
            aEdmCluster.addToHits(cellclone);
          }
        }
      }
    }
  } else {
    for (int iTheta = thetaId - halfThetaFin; iTheta <= int(thetaId + halfThetaFin); iTheta++) {
      for (int iPhi = phiId - halfPhiFin; iPhi <= int(phiId + halfPhiFin); iPhi++) {
        for (auto cell : m_cellsInTowers[std::make_pair(iTheta, phiNeighbour(iPhi))]) {
          if (std::find(seen_cellIDs.begin(), seen_cellIDs.end(), cell.getCellID()) != seen_cellIDs.end()) { // towers can be smaller than cells in which case a cell belongs to several towers
            continue;
          }
          seen_cellIDs.push_back(cell.getCellID());
          auto cellclone = cell.clone();
          aEdmClusterCells->push_back(cellclone);
          aEdmCluster.addToHits(cellclone);
        }
      }
    }
  }
  return;
}
