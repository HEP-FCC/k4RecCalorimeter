#include "CaloTowerToolFCCee.h"

// edm4hep
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/Cluster.h"
#include "edm4hep/MutableCluster.h"

// DD4hep
#include "DD4hep/Detector.h"

DECLARE_COMPONENT(CaloTowerToolFCCee)

CaloTowerToolFCCee::CaloTowerToolFCCee(const std::string& type, const std::string& name, const IInterface* parent)
    : AlgTool(type, name, parent) {
  declareInterface<ITowerToolThetaModule>(this);
}

StatusCode CaloTowerToolFCCee::initialize() {
  if (AlgTool::initialize().isFailure()) {
    return StatusCode::FAILURE;
  }

  // create handles for input cell collections
  for (const auto& col : m_cellCollections) {
    debug() << "Creating handle for input cell (CalorimeterHit) collection : " << col << endmsg;
    try {
      m_cellCollectionHandles.push_back(
          new DataHandle<edm4hep::CalorimeterHitCollection>(col, Gaudi::DataHandle::Reader, this));
    } catch (...) {
      error() << "Error creating handle for input collection: " << col << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // very small number (epsilon) substructed from the edges to ensure correct division in case of rounding erros
  // (GM: issue not seen in my tests but kept to be on the safe side)
  double epsilon = 0.0001;
  // number of phi bins
  m_nPhiTower = ceil((m_phiMax - m_phiMin - epsilon) / m_deltaPhiTower);
  // number of theta bins
  m_nThetaTower = ceil((m_thetaMax - m_thetaMin - epsilon) / m_deltaThetaTower);
  debug() << "Towers: thetaMin " << m_thetaMin.value() <<  ", thetaMax " << m_thetaMax.value()
          << ", deltaThetaTower " << m_deltaThetaTower.value() << ", nThetaTower " << m_nThetaTower << endmsg;
  debug() << "Towers: phiMin " << m_phiMin.value() << ", phiMax " << m_phiMax.value()
          << ", deltaPhiTower " << m_deltaPhiTower.value() << ", nPhiTower " << m_nPhiTower << endmsg;

  return StatusCode::SUCCESS;
}

StatusCode CaloTowerToolFCCee::finalize() {
  for (auto& towerInMap : m_cellsInTowers) {
    towerInMap.second.clear();
  }

  for (size_t ih = 0; ih < m_cellCollectionHandles.size(); ih++)
    delete m_cellCollectionHandles[ih];
  return AlgTool::finalize();
}

void CaloTowerToolFCCee::towersNumber(int& nTheta, int& nPhi) {
  nTheta = m_nThetaTower;
  nPhi = m_nPhiTower;
}

uint CaloTowerToolFCCee::buildTowers(std::vector<std::vector<float>>& aTowers, bool fillTowersCells) {
  uint totalNumberOfCells = 0;
  uint totalNumberOfClusteredCells = 0;
  for (auto& towerInMap : m_cellsInTowers) {
    towerInMap.second.clear();
  }

  // Loop over input cell collections to build towers
  for (size_t ih = 0; ih < m_cellCollectionHandles.size(); ih++) {
    verbose() << "Processing collection " << ih << endmsg;
    const edm4hep::CalorimeterHitCollection* coll = m_cellCollectionHandles[ih]->get();
    debug() << "Input cell collection size: " << coll->size() << endmsg;
    // Loop over collection of calorimeter cells
    if (coll->size() > 0) {
      totalNumberOfClusteredCells += CellsIntoTowers(aTowers, coll, fillTowersCells);
      totalNumberOfCells += coll->size();
    }
  }

  debug() << "Total number of input cells: " << totalNumberOfCells << endmsg;
  debug() << "Total number of clustered input cells: " << totalNumberOfClusteredCells << endmsg;

  return totalNumberOfClusteredCells;
}

// Get the tower IDs in theta
// aTheta Position of the calorimeter cell in theta
// Note that the function returns an unsigned int so it
// assumes that aTheta is <= thetaMax.
// This is checked in CellsIntoTowers
uint CaloTowerToolFCCee::idTheta(float aTheta) const {
  uint id = floor((m_thetaMax - aTheta) / m_deltaThetaTower);
  return id;
}

// Get the tower IDs in phi
// aPhi Position of the calorimeter cell in phi
// Note hat the function returns an unsigned int so it
// assumes that aPhi is >= phiMin.
// This is checked in CellsIntoTowers
uint CaloTowerToolFCCee::idPhi(float aPhi) const {
  uint id = floor((aPhi - m_phiMin) / m_deltaPhiTower);
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
  return (m_phiMin + (aIdPhi + 0.5) * m_deltaPhiTower);
}

// Return phi index of tower, taking into account
// that phi is cyclic (if phimax - phimin=2pi, can have towers with
// aIPhi < 0 or >= m_nPhiTower)
uint CaloTowerToolFCCee::phiIndexTower(int aIPhi) const {
  if (aIPhi < 0) {
    return phiIndexTower(aIPhi + ceil(2*M_PI/m_deltaPhiTower));
  } else if (aIPhi >= m_nPhiTower) {
    return phiIndexTower(aIPhi - ceil(2*M_PI/m_deltaPhiTower));
  } else
    return aIPhi;
}

std::map<std::pair<uint, uint>, std::vector<edm4hep::CalorimeterHit>> CaloTowerToolFCCee::cellsInTowers() const {
  return m_cellsInTowers;
}

// to fill the cell infomation into towers
uint CaloTowerToolFCCee::CellsIntoTowers(std::vector<std::vector<float>>& aTowers,
                                         const edm4hep::CalorimeterHitCollection* aCells, bool fillTowersCells) {
  // Loop over a collection of calorimeter cells and build calo towers
  // tower index of the borders of the cell
  int iTheta = 0;
  int iPhi = 0;
  uint clusteredCells = 0;

  for (const auto& cell : *aCells) {
    float cellX = cell.getPosition().x;
    float cellY = cell.getPosition().y;
    float cellZ = cell.getPosition().z;
    float cellTheta = atan2(sqrt(cellX * cellX + cellY * cellY), cellZ);
    float cellPhi = atan2(cellY, cellX);
    // skip cells outside of specified ranges
    if (cellTheta<m_thetaMin || cellTheta>m_thetaMax) {
      warning() << "Cell theta " << cellTheta << " outside of theta range of towers, will not be clustered" << endmsg;
      continue;
    }
    if (cellPhi<m_phiMin || cellPhi>m_phiMax) {
      warning() << "Cell phi " << cellPhi << " outside of phi range of towers, will not be clustered" << endmsg;
      continue;
    }
    iTheta = idTheta(cellTheta);
    iPhi = idPhi(cellPhi);
    //debug() << "Cell: x = " << cellX << " y = " << cellY << " z = " << cellZ << endmsg;
    //debug() << "Cell: theta = " << cellTheta << " phi = " << cellPhi << endmsg;
    //debug() << "Cell: iTheta = " << iTheta << " iPhi = " << iPhi << " iPhi(cyclic) = " << phiIndexTower(iPhi) << endmsg;
    aTowers[iTheta][phiIndexTower(iPhi)] += cell.getEnergy() * sin(cellTheta);
    if (fillTowersCells) {
      clusteredCells++;
      m_cellsInTowers[std::make_pair(iTheta, phiIndexTower(iPhi))].push_back(cell);
      int ncells = m_cellsInTowers[std::make_pair(iTheta, phiIndexTower(iPhi))].size();

      if (ncells > 5)
        verbose() << "NUM CELLs IN TOWER : " << ncells << endmsg;
    }
  }

  return clusteredCells;
}

void CaloTowerToolFCCee::attachCells(float theta, float phi, uint halfThetaFin, uint halfPhiFin,
                                     edm4hep::MutableCluster& aEdmCluster,
                                     edm4hep::CalorimeterHitCollection* aEdmClusterCells, bool aEllipse) {
  int thetaId = idTheta(theta);
  int phiId = idPhi(phi);
  std::vector<dd4hep::DDSegmentation::CellID> seen_cellIDs;

  std::vector<float> subDetectorEnergies(m_nSubDetectors);

  if (aEllipse) {
    for (int iTheta = thetaId - halfThetaFin; iTheta <= int(thetaId + halfThetaFin); iTheta++) {
      for (int iPhi = phiId - halfPhiFin; iPhi <= int(phiId + halfPhiFin); iPhi++) {
        if (pow((thetaId - iTheta) / (halfThetaFin + 0.5), 2) + pow((phiId - iPhi) / (halfPhiFin + 0.5), 2) < 1) {
          for (auto cell : m_cellsInTowers[std::make_pair(iTheta, phiIndexTower(iPhi))]) {
            if (std::find(seen_cellIDs.begin(), seen_cellIDs.end(), cell.getCellID()) !=
                seen_cellIDs.end()) { // towers can be smaller than cells in which case a cell belongs to several towers
              continue;
            }
            seen_cellIDs.push_back(cell.getCellID());
            // if aEdmClusterCells it not nullptr, the user wants the clustered cells to be put into a new collection
            // otherwise we just set links to the existing cells
            if (aEdmClusterCells) {
              auto cellclone = cell.clone();
              aEdmClusterCells->push_back(cellclone);
              aEdmCluster.addToHits(cellclone);
            }
            else {
              aEdmCluster.addToHits(cell);
            }

            if (m_nSubDetectors>0) {
              // caloID: 1 = ecal, 2 = hcal, 3 = yoke - see how m_caloid is computed and encoded in cell type in
              // https://github.com/HEP-FCC/k4RecCalorimeter/blob/main/RecCalorimeter/src/components/CreatePositionedCaloCells.cpp
              int caloID = ((cell.getType() / 10) % 10) -1 ;
              if (caloID < 0 or caloID > (int) m_nSubDetectors) {
                warning() << "Wrong caloID " << caloID << endmsg;
              }
              else {
                subDetectorEnergies[caloID] += cell.getEnergy();
              }
            }
          }
        }
      }
    }
  } else {
    for (int iTheta = thetaId - halfThetaFin; iTheta <= int(thetaId + halfThetaFin); iTheta++) {
      for (int iPhi = phiId - halfPhiFin; iPhi <= int(phiId + halfPhiFin); iPhi++) {
        for (auto cell : m_cellsInTowers[std::make_pair(iTheta, phiIndexTower(iPhi))]) {
          if (std::find(seen_cellIDs.begin(), seen_cellIDs.end(), cell.getCellID()) !=
              seen_cellIDs.end()) { // towers can be smaller than cells in which case a cell belongs to several towers
            continue;
          }
          seen_cellIDs.push_back(cell.getCellID());
          if (aEdmClusterCells) {
            auto cellclone = cell.clone();
            aEdmClusterCells->push_back(cellclone);
            aEdmCluster.addToHits(cellclone);
          }
          else {
            aEdmCluster.addToHits(cell);
          }
          if (m_nSubDetectors>0) {
            int caloID = ((cell.getType() / 10) % 10 - 1);
            if (caloID < 0 or caloID > (int) m_nSubDetectors) {
              warning() << "Wrong caloID " << caloID << endmsg;
            }
            else {
              subDetectorEnergies[caloID] += cell.getEnergy();
            }
          }
        }
      }
    }
  }
  for (auto energy: subDetectorEnergies) {
    aEdmCluster.addToSubdetectorEnergies(energy);
  }
  return;
}
