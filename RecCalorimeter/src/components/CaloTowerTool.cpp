#include "CaloTowerTool.h"

// FCCSW
#include "k4Interface/IGeoSvc.h"

// datamodel
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/Cluster.h"
#include "edm4hep/MutableCluster.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DD4hep/Readout.h"

DECLARE_COMPONENT(CaloTowerTool)

CaloTowerTool::CaloTowerTool(const std::string& type, const std::string& name, const IInterface* parent)
    : GaudiTool(type, name, parent), m_geoSvc("GeoSvc", name) {
  declareProperty("ecalBarrelCells", m_ecalBarrelCells, "");
  declareProperty("ecalEndcapCells", m_ecalEndcapCells, "");
  declareProperty("ecalFwdCells", m_ecalFwdCells, "");
  declareProperty("hcalBarrelCells", m_hcalBarrelCells, "");
  declareProperty("hcalExtBarrelCells", m_hcalExtBarrelCells, "");
  declareProperty("hcalEndcapCells", m_hcalEndcapCells, "");
  declareProperty("hcalFwdCells", m_hcalFwdCells, "");
  declareProperty("positionsECalBarrelTool", m_cellPositionsECalBarrelTool,
                  "Handle for tool to retrieve cell positions in ECal Barrel");

  declareInterface<ITowerTool>(this);
}

StatusCode CaloTowerTool::initialize() {
  if (GaudiTool::initialize().isFailure()) {
    return StatusCode::FAILURE;
  }
  
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service. "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }

  // check if readouts exist & retrieve Module-Theta segmentations
  // if readout does not exist, reconstruction without this calorimeter part will be performed
  std::pair<dd4hep::DDSegmentation::Segmentation*, SegmentationType> tmpPair;
  info() << "Retrieving Ecal barrel segmentation" << endmsg;
  tmpPair = retrieveSegmentation(m_ecalBarrelReadoutName);
  m_ecalBarrelSegmentation = tmpPair.first;
  m_ecalBarrelSegmentationType = tmpPair.second;
  if (tmpPair.first != nullptr && tmpPair.second == SegmentationType::kWrong) {
    error() << "Wrong type of segmentation" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_useHalfTower) {
    m_decoder = m_geoSvc->lcdd()->readout(m_ecalBarrelReadoutName).idSpec().decoder();
  }
  //info() << "Retrieving Ecal endcap segmentation" << endmsg;
  //tmpPair = retrieveSegmentation(m_ecalEndcapReadoutName);
  //m_ecalEndcapSegmentation = tmpPair.first;
  //m_ecalEndcapSegmentationType = tmpPair.second;
  //if (tmpPair.first != nullptr && tmpPair.second == SegmentationType::kWrong) {
  //  error() << "Wrong type of segmentation" << endmsg;
  //  return StatusCode::FAILURE;
  //}

  //info() << "Retrieving Ecal forward segmentation" << endmsg;
  //tmpPair = retrieveSegmentation(m_ecalFwdReadoutName);
  //m_ecalFwdSegmentation = tmpPair.first;
  //m_ecalFwdSegmentationType = tmpPair.second;
  //if (tmpPair.first != nullptr && tmpPair.second == SegmentationType::kWrong) {
  //  error() << "Wrong type of segmentation" << endmsg;
  //  return StatusCode::FAILURE;
  //}
  //info() << "Retrieving Hcal barrel segmentation" << endmsg;
  //tmpPair = retrieveSegmentation(m_hcalBarrelReadoutName);
  //m_hcalBarrelSegmentation = tmpPair.first;
  //m_hcalBarrelSegmentationType = tmpPair.second;
  //if (tmpPair.first != nullptr && tmpPair.second == SegmentationType::kWrong) {
  //  error() << "Wrong type of segmentation" << endmsg;
  //  return StatusCode::FAILURE;
  //}
  //info() << "Retrieving Hcal extended barrel segmentation" << endmsg;
  //tmpPair = retrieveSegmentation(m_hcalExtBarrelReadoutName);
  //m_hcalExtBarrelSegmentation = tmpPair.first;
  //m_hcalExtBarrelSegmentationType = tmpPair.second;
  //if (tmpPair.first != nullptr && tmpPair.second == SegmentationType::kWrong) {
  //  error() << "Wrong type of segmentation" << endmsg;
  //  return StatusCode::FAILURE;
  //}
  //info() << "Retrieving Hcal endcap segmentation" << endmsg;
  //tmpPair = retrieveSegmentation(m_hcalEndcapReadoutName);
  //m_hcalEndcapSegmentation = tmpPair.first;
  //m_hcalEndcapSegmentationType = tmpPair.second;
  //if (tmpPair.first != nullptr && tmpPair.second == SegmentationType::kWrong) {
  //  error() << "Wrong type of segmentation" << endmsg;
  //  return StatusCode::FAILURE;
  //}
  //info() << "Retrieving Hcal forward segmentation" << endmsg;
  //tmpPair = retrieveSegmentation(m_hcalFwdReadoutName);
  //m_hcalFwdSegmentation = tmpPair.first;
  //m_hcalFwdSegmentationType = tmpPair.second;
  //if (tmpPair.first != nullptr && tmpPair.second == SegmentationType::kWrong) {
  //  error() << "Wrong type of segmentation" << endmsg;
  //  return StatusCode::FAILURE;
  //}
  return StatusCode::SUCCESS;
}

StatusCode CaloTowerTool::finalize() { 
  for (auto& towerInMap : m_cellsInTowers) {
    towerInMap.second.clear();
    }
  return GaudiTool::finalize();
  }

std::pair<double, double> CaloTowerTool::retrievePhiThetaExtrema(dd4hep::DDSegmentation::Segmentation* aSegmentation, SegmentationType aType) {
  double phiMax = -1;
  double thetaMax = -1;
  dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged* segmentation = nullptr;
  if (aSegmentation != nullptr) {
    switch (aType) {
    case SegmentationType::kModuleTheta: {
      info() << "== Retrieving segmentation " << aSegmentation->name() << endmsg;
      segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged*>(aSegmentation);
      phiMax = M_PI - M_PI/segmentation->nModules();
      thetaMax = M_PI - fabs(segmentation->offsetTheta()) + segmentation->gridSizeTheta() * .5;
      break;
    }
    case SegmentationType::kMulti: {
      double phi = -1;
      double theta = -1;
      info() << "== Retrieving multi segmentation " << aSegmentation->name() << endmsg;
      for (const auto& subSegm: dynamic_cast<dd4hep::DDSegmentation::MultiSegmentation*>(aSegmentation)->subSegmentations()) {
        segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged*>(subSegm.segmentation);
        phi = M_PI - M_PI/segmentation->nModules();
        theta = M_PI - fabs(segmentation->offsetTheta()) + segmentation->gridSizeTheta() * .5;
        if (theta > thetaMax) { thetaMax = theta;}
        if (phi > phiMax) { phiMax = phi;}
      }
      break;
    }
    case SegmentationType::kWrong: {
      info() << "== Retrieving WRONG segmentation" << endmsg;
      phiMax = -1;
      thetaMax = -1;
      break;
    }
    }
  }
  return std::make_pair(phiMax, thetaMax);
}

tower CaloTowerTool::towersNumber() {

  std::vector<double> listPhiMax;
  std::vector<double> listThetaMax;
  //listPhiMax.reserve(7);
  //listThetaMax.reserve(7);
  listPhiMax.reserve(1);
  listThetaMax.reserve(1);

  std::pair<double, double> tmpPair;
  tmpPair = retrievePhiThetaExtrema(m_ecalBarrelSegmentation, m_ecalBarrelSegmentationType);
  listPhiMax.push_back(tmpPair.first);
  listThetaMax.push_back(tmpPair.second);
  //tmpPair = retrievePhiThetaExtrema(m_ecalEndcapSegmentation, m_ecalEndcapSegmentationType);
  //listPhiMax.push_back(tmpPair.first);
  //listThetaMax.push_back(tmpPair.second);
  //tmpPair = retrievePhiThetaExtrema(m_ecalFwdSegmentation, m_ecalFwdSegmentationType);
  //listPhiMax.push_back(tmpPair.first);
  //listThetaMax.push_back(tmpPair.second);
  //tmpPair = retrievePhiThetaExtrema(m_hcalBarrelSegmentation, m_hcalBarrelSegmentationType);
  //listPhiMax.push_back(tmpPair.first);
  //listThetaMax.push_back(tmpPair.second);
  //tmpPair = retrievePhiThetaExtrema(m_hcalExtBarrelSegmentation, m_hcalExtBarrelSegmentationType);
  //listPhiMax.push_back(tmpPair.first);
  //listThetaMax.push_back(tmpPair.second);
  //tmpPair = retrievePhiThetaExtrema(m_hcalEndcapSegmentation, m_hcalEndcapSegmentationType);
  //listPhiMax.push_back(tmpPair.first);
  //listThetaMax.push_back(tmpPair.second);
  //tmpPair = retrievePhiThetaExtrema(m_hcalFwdSegmentation, m_hcalFwdSegmentationType);
  //listPhiMax.push_back(tmpPair.first);
  //listThetaMax.push_back(tmpPair.second);

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

  tower total;
  total.theta = m_nThetaTower;
  total.phi = m_nPhiTower;
  return total;
}

uint CaloTowerTool::buildTowers(std::vector<std::vector<float>>& aTowers, bool fillTowersCells) {
  uint totalNumberOfCells = 0;
  for (auto& towerInMap : m_cellsInTowers) {
    towerInMap.second.clear();
  }
  // 1. ECAL barrel
  // Get the input collection with calorimeter cells
  const edm4hep::CalorimeterHitCollection* ecalBarrelCells = m_ecalBarrelCells.get();
  debug() << "Input Ecal barrel cell collection size: " << ecalBarrelCells->size() << endmsg;
  // Loop over a collection of calorimeter cells and build calo towers
  if (m_ecalBarrelSegmentation != nullptr) {
    CellsIntoTowers(aTowers, ecalBarrelCells, m_ecalBarrelSegmentation, m_ecalBarrelSegmentationType, fillTowersCells);
    totalNumberOfCells += ecalBarrelCells->size();
  }
/*
  // 2. ECAL endcap calorimeter
  const edm4hep::CalorimeterHitCollection* ecalEndcapCells = m_ecalEndcapCells.get();
  debug() << "Input Ecal endcap cell collection size: " << ecalEndcapCells->size() << endmsg;
  // Loop over a collection of calorimeter cells and build calo towers
  if (m_ecalEndcapSegmentation != nullptr) {
    CellsIntoTowers(aTowers, ecalEndcapCells, m_ecalEndcapSegmentation, m_ecalEndcapSegmentationType, fillTowersCells);
    totalNumberOfCells += ecalEndcapCells->size();
  }

  // 3. ECAL forward calorimeter
  const edm4hep::CalorimeterHitCollection* ecalFwdCells = m_ecalFwdCells.get();
  debug() << "Input Ecal forward cell collection size: " << ecalFwdCells->size() << endmsg;
  // Loop over a collection of calorimeter cells and build calo towers
  if (m_ecalFwdSegmentation != nullptr) {
    CellsIntoTowers(aTowers, ecalFwdCells, m_ecalFwdSegmentation, m_ecalFwdSegmentationType, fillTowersCells);
    totalNumberOfCells += ecalFwdCells->size();
  }

  // 4. HCAL barrel
  const edm4hep::CalorimeterHitCollection* hcalBarrelCells = m_hcalBarrelCells.get();
  debug() << "Input hadronic barrel cell collection size: " << hcalBarrelCells->size() << endmsg;
  // Loop over a collection of calorimeter cells and build calo towers
  if (m_hcalBarrelSegmentation != nullptr) {
    CellsIntoTowers(aTowers, hcalBarrelCells, m_hcalBarrelSegmentation, m_hcalBarrelSegmentationType, fillTowersCells);
    totalNumberOfCells += hcalBarrelCells->size();
  }

  // 5. HCAL extended barrel
  const edm4hep::CalorimeterHitCollection* hcalExtBarrelCells = m_hcalExtBarrelCells.get();
  debug() << "Input hadronic extended barrel cell collection size: " << hcalExtBarrelCells->size() << endmsg;
  // Loop over a collection of calorimeter cells and build calo towers
  if (m_hcalExtBarrelSegmentation != nullptr) {
    CellsIntoTowers(aTowers, hcalExtBarrelCells, m_hcalExtBarrelSegmentation, m_hcalExtBarrelSegmentationType, fillTowersCells);
    totalNumberOfCells += hcalExtBarrelCells->size();
  }

  // 6. HCAL endcap calorimeter
  const edm4hep::CalorimeterHitCollection* hcalEndcapCells = m_hcalEndcapCells.get();
  debug() << "Input Hcal endcap cell collection size: " << hcalEndcapCells->size() << endmsg;
  // Loop over a collection of calorimeter cells and build calo towers
  if (m_hcalEndcapSegmentation != nullptr) {
    CellsIntoTowers(aTowers, hcalEndcapCells, m_hcalEndcapSegmentation, m_hcalEndcapSegmentationType, fillTowersCells);
    totalNumberOfCells += hcalEndcapCells->size();
  }

  // 7. HCAL forward calorimeter
  const edm4hep::CalorimeterHitCollection* hcalFwdCells = m_hcalFwdCells.get();
  debug() << "Input Hcal forward cell collection size: " << hcalFwdCells->size() << endmsg;
  // Loop over a collection of calorimeter cells and build calo towers
  if (m_hcalFwdSegmentation != nullptr) {
    CellsIntoTowers(aTowers, hcalFwdCells, m_hcalFwdSegmentation, m_hcalFwdSegmentationType, fillTowersCells);
    totalNumberOfCells += hcalFwdCells->size();
  }
*/
  return totalNumberOfCells;
}

// Get the tower IDs in theta
// aTheta Position of the calorimeter cell in theta
uint CaloTowerTool::idTheta(float aTheta) const {
  uint id = floor(( m_thetaMax - aTheta ) / m_deltaThetaTower);
  return id;
}

// Get the tower IDs in phi
// aPhi Position of the calorimeter cell in phi
uint CaloTowerTool::idPhi(float aPhi) const {
  uint id = floor((aPhi + m_phiMax) / m_deltaPhiTower);
  return id;
}

// Get the theta position of the centre of the tower
float CaloTowerTool::theta(int aIdTheta) const {
  // middle of the tower
  return (m_thetaMax - (aIdTheta + 0.5) * m_deltaThetaTower);
}

// Get the phi position of the centre of the tower
float CaloTowerTool::phi(int aIdPhi) const {
  // middle of the tower
  return ((aIdPhi + 0.5) * m_deltaPhiTower - m_phiMax);
}

uint CaloTowerTool::phiNeighbour(int aIPhi) const {
  if (aIPhi < 0) {
    return m_nPhiTower + aIPhi;
  } else if (aIPhi >= m_nPhiTower) {
    return aIPhi % m_nPhiTower;
  }
  return aIPhi;
}

float CaloTowerTool::radiusForPosition() const { return m_radius; }

std::map<std::pair<uint, uint>, std::vector<edm4hep::CalorimeterHit>> CaloTowerTool::cellsInTowers() const { return m_cellsInTowers; }

// to fill the cell infomation into towers
void CaloTowerTool::CellsIntoTowers(std::vector<std::vector<float>>& aTowers,
                                    const edm4hep::CalorimeterHitCollection* aCells,
                                    dd4hep::DDSegmentation::Segmentation* aSegmentation, SegmentationType aType,
                                    bool fillTowersCells) {
  // Loop over a collection of calorimeter cells and build calo towers
  // tower index of the borders of the cell
  int iTheta = 0;
  int iPhi = 0;
  const dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged* segmentation = nullptr;
  const dd4hep::DDSegmentation::MultiSegmentation* multisegmentation = nullptr;
  if( aType == SegmentationType::kModuleTheta) {
      segmentation = dynamic_cast<const dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged*>(aSegmentation);
  } else if( aType == SegmentationType::kMulti) {
    multisegmentation = dynamic_cast<const dd4hep::DDSegmentation::MultiSegmentation*>(aSegmentation);
  }
  bool pass = true;
  for (const auto& cell : *aCells) {
    pass = true;
    // if multisegmentation is used - first find out which segmentation to use
    if( aType == SegmentationType::kMulti) {
      segmentation = dynamic_cast<const dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged*>(&multisegmentation->subsegmentation(cell.getCellID()));
    }
    if (m_useHalfTower) {
        uint layerId = m_decoder->get(cell.getCellID(), "layer");
        if ( layerId > m_max_layer ) {
          pass = false;
        }
    }
    if (pass) {
      dd4hep::Position posCell;
      posCell = m_cellPositionsECalBarrelTool->xyzPosition(cell.getCellID());
      // find to which tower(s) the cell belongs
      float cellTheta_old = segmentation->theta(cell.getCellID());
      float cellTheta = posCell.Theta();
      float cellPhi   = posCell.Phi();
      iTheta = idTheta(cellTheta);
      iPhi = idPhi(cellPhi);
      aTowers[iTheta][phiNeighbour(iPhi)] +=
        cell.getEnergy() * sin(segmentation->theta(cell.getCellID()));
      if (fillTowersCells) {
        m_cellsInTowers[std::make_pair(iTheta, phiNeighbour(iPhi))].push_back(cell.clone());
        if ( m_cellsInTowers[std::make_pair(iTheta, phiNeighbour(iPhi))].size() > 5 )
          verbose() << "NUM CELLs IN TOWER : " << m_cellsInTowers[std::make_pair(iTheta, phiNeighbour(iPhi))].size() << endmsg;
      }
    }
  }
}

std::pair<dd4hep::DDSegmentation::Segmentation*, CaloTowerTool::SegmentationType> CaloTowerTool::retrieveSegmentation(std::string aReadoutName) {
  dd4hep::DDSegmentation::Segmentation* segmentation = nullptr;
  if (m_geoSvc->lcdd()->readouts().find(aReadoutName) == m_geoSvc->lcdd()->readouts().end()) {
    info() << "Readout does not exist! Please check if it is correct. Processing without it." << endmsg;
  } else {
    info() << "Readout " << aReadoutName << " found." << endmsg;
    segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged*>(
      m_geoSvc->lcdd()->readout(aReadoutName).segmentation().segmentation());
    if (segmentation == nullptr) {
      segmentation = dynamic_cast<dd4hep::DDSegmentation::MultiSegmentation*>(
        m_geoSvc->lcdd()->readout(aReadoutName).segmentation().segmentation());
      if (segmentation == nullptr) {
        error() << "There is no module-theta or multi- segmentation for the readout " << aReadoutName << " defined." << endmsg;
      } else {
        // check if multisegmentation contains only module-theta sub-segmentations
        dd4hep::DDSegmentation::Segmentation* subsegmentation = nullptr;
        for (const auto& subSegm: dynamic_cast<dd4hep::DDSegmentation::MultiSegmentation*>(segmentation)->subSegmentations()) {
          subsegmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged*>(subSegm.segmentation);
          if (subsegmentation == nullptr) {
            error() << "At least one of the sub-segmentations in MultiSegmentation named " << aReadoutName << " is not a module-theta grid." << endmsg;
            return std::make_pair(nullptr, SegmentationType::kWrong);
          }
        }
        return std::make_pair(segmentation, SegmentationType::kMulti);
      }
    } else {
      return std::make_pair(segmentation, SegmentationType::kModuleTheta);
    }
  }
  return std::make_pair(segmentation, SegmentationType::kWrong);
}

void CaloTowerTool::attachCells(float theta, float phi, uint halfThetaFin, uint halfPhiFin, edm4hep::MutableCluster& aEdmCluster, edm4hep::CalorimeterHitCollection* aEdmClusterCells, bool aEllipse) {
  int thetaId = idTheta(theta);
  int phiId = idPhi(phi);
  int num1 = 0;
  int num2 = 0;
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
            num1++;
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
          num2++;
        }
      }
    }
  }
  return;
}
