#include "NoiseCaloCellsFlatTool.h"
#include <GaudiKernel/StatusCode.h>

DECLARE_COMPONENT(NoiseCaloCellsFlatTool)

NoiseCaloCellsFlatTool::NoiseCaloCellsFlatTool(const std::string& type, const std::string& name,
                                               const IInterface* parent)
    : AlgTool(type, name, parent) {
  declareInterface<INoiseCaloCellsTool>(this);
}

StatusCode NoiseCaloCellsFlatTool::initialize() {
  {
    StatusCode sc = AlgTool::initialize();
    if (sc.isFailure()) return sc;
  }

  // Initialize random service
  m_randSvc = service("RndmGenSvc", false);
  if (!m_randSvc) {
    error() << "Couldn't get RndmGenSvc" << endmsg;
    return StatusCode::FAILURE;
  }
  {
    StatusCode sc = m_gauss.initialize(m_randSvc, Rndm::Gauss(0., 1.));
    if (sc.isFailure()) {
      error() << "Failed to initialize Gaussian random number generator!" << endmsg;
    }
  }

  info() << "RMS of the cell noise: " << m_cellNoiseRMS * 1.e3 << " MeV" << endmsg;
  info() << "Offset of the cell noise: " << m_cellNoiseOffset * 1.e3 << " MeV" << endmsg;
  info() << "Filter noise threshold: " << m_filterThreshold << "*sigma" << endmsg;
  return StatusCode::SUCCESS;
}

void NoiseCaloCellsFlatTool::addRandomCellNoise(std::unordered_map<uint64_t, double>& aCells) {
  std::for_each(aCells.begin(), aCells.end(),
                [this](std::pair<const uint64_t, double>& p) { p.second += (m_cellNoiseOffset + (m_gauss.shoot() * m_cellNoiseRMS)); });
}

void NoiseCaloCellsFlatTool::filterCellNoise(std::unordered_map<uint64_t, double>& aCells) {
  // Erase a cell if it has energy below a threshold
  double threshold = m_cellNoiseOffset + m_filterThreshold * m_cellNoiseRMS;
  auto it = aCells.begin();
  while ((it = std::find_if(it, aCells.end(), [&threshold](std::pair<const uint64_t, double>& p) {
            return bool(p.second < threshold);
          })) != aCells.end()) {
    aCells.erase(it++);
  }
}

StatusCode NoiseCaloCellsFlatTool::finalize() { return AlgTool::finalize(); }
