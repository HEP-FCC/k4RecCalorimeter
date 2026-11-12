#include "NoiseCaloCellsFromFileBarrelTool.h"
#include "k4FWCore/GaudiChecks.h"

DECLARE_COMPONENT(NoiseCaloCellsFromFileBarrelTool)

StatusCode NoiseCaloCellsFromFileBarrelTool::initialize() {
  K4_GAUDI_CHECK(NoiseCaloCellsFromFileBaseTool::initialize());
  K4_GAUDI_CHECK(m_cellPositionsTool.retrieve());
  return StatusCode::SUCCESS;
}

StatusCode NoiseCaloCellsFromFileBarrelTool::initBinning(NoiseData& data,
                                                         const k4::recCalo::ICaloIndexer& indexer) const {
  data.m_bins.resize(indexer.cellIDs().size());

  // All histograms have same binning, all bins with same size
  // Using the histograms in the first layer to define binning.
  const TH1* h_rms = data.m_histoElecNoiseRMS.at(0).get();
  const TH1* h_offset = m_setNoiseOffset ? data.m_histoElecNoiseOffset.at(0).get() : nullptr;
  for (uint64_t id : indexer.cellIDs()) {
    unsigned ndx = indexer.index(id);
    double cellTheta = m_cellPositionsTool->xyzPosition(id).Theta();
    int ibinRMS = h_rms->FindFixBin(cellTheta);
    int ibinOffset = h_offset ? h_offset->FindFixBin(cellTheta) : 0;
    data.m_bins.at(ndx) = std::make_pair(ibinRMS, ibinOffset);
  }

  return StatusCode::SUCCESS;
}
