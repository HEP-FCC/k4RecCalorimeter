#include "NoiseCaloCellsFromFileBaseTool.h"
#include "k4FWCore/GaudiChecks.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"

// k4FWCore
#include "k4Interface/IGeoSvc.h"

// DD4hep
#include "DD4hep/Detector.h"

// ROOT
#include "TFile.h"
#include "TSystem.h"

StatusCode NoiseCaloCellsFromFileBaseTool::initialize() {
  m_elecNoiseHistoName = m_elecNoiseRMSHistoName;

  K4_GAUDI_CHECK(base_class::initialize());
  K4_GAUDI_CHECK(m_randSvc = service<IRndmGenSvc>("RndmGenSvc", false));
  K4_GAUDI_CHECK(m_gauss.initialize(m_randSvc, Rndm::Gauss(0., 1.)));

  debug() << "Filter noise threshold: " << m_filterThreshold << "*sigma" << endmsg;

  return StatusCode::SUCCESS;
}

template <class C>
void NoiseCaloCellsFromFileBaseTool::addRandomCellNoiseT(C& aCells) const {
  for (auto& p : aCells) {
    auto [rms, offset] = getNoisePerCell(p.first);
    p.second += offset + rms * m_gauss.shoot();
  }
}

void NoiseCaloCellsFromFileBaseTool::addRandomCellNoise(std::unordered_map<CellID, double>& aCells) const {
  addRandomCellNoiseT(aCells);
}

void NoiseCaloCellsFromFileBaseTool::addRandomCellNoise(std::vector<std::pair<CellID, double>>& aCells) const {
  addRandomCellNoiseT(aCells);
}

template <typename C>
void NoiseCaloCellsFromFileBaseTool::filterCellNoiseT(C& aCells) const {
  // Erase a cell if it has energy below a threshold from the vector
  std::erase_if(aCells, [&](auto& p) {
    auto [rms, offset] = getNoisePerCell(p.first);
    double e = p.second - offset;
    if (m_useAbsInFilter)
      e = std::abs(e);
    return e < m_filterThreshold * rms;
  });
}

void NoiseCaloCellsFromFileBaseTool::filterCellNoise(std::unordered_map<CellID, double>& aCells) const {
  filterCellNoiseT(aCells);
}

void NoiseCaloCellsFromFileBaseTool::filterCellNoise(std::vector<std::pair<CellID, double>>& aCells) const {
  filterCellNoiseT(aCells);
}
