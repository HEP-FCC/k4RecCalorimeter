#include "NoiseCaloCellsFromFileTurbineEndcapTool.h"

DECLARE_COMPONENT(NoiseCaloCellsFromFileTurbineEndcapTool)

// clang-format off
unsigned NoiseCaloCellsFromFileTurbineEndcapTool::getBin(const char* what,
                                                         const TH1& h,
                                                         unsigned iRho,
                                                         unsigned iZ) const {
  // clang-format on
  unsigned NbinsZ = h.GetNbinsX();
  unsigned NbinsRho = h.GetNbinsY();

  unsigned ibin = 0;
  if (iRho > NbinsRho || iZ > NbinsZ) {
    error() << "bins outside range of the " << what << " histograms! Bins: " << iRho << "," << iZ
            << " , Nbins in histogram: " << NbinsRho << "," << NbinsZ << endmsg;
  } else {
    ibin = h.GetBin(iZ, iRho);
  }
  return ibin;
}

StatusCode NoiseCaloCellsFromFileTurbineEndcapTool::initBinning(NoiseData& data,
                                                                const k4::recCalo::ICaloIndexer& indexer) const {
  data.m_bins.resize(indexer.cellIDs().size());

  for (uint64_t id : indexer.cellIDs()) {
    unsigned ndx = indexer.index(id);
    unsigned iHist = m_decoder->get(id, m_index_activeField);
    unsigned iRho = m_decoder->get(id, "rho") + 1;
    unsigned iZ = m_decoder->get(id, "z") + 1;

    unsigned ibinRMS = getBin("RMS", *data.m_histoElecNoiseRMS.at(iHist), iRho, iZ);
    unsigned ibinOffset = 0;
    if (m_setNoiseOffset)
      ibinOffset = getBin("Offset", *data.m_histoElecNoiseOffset.at(iHist), iRho, iZ);

    data.m_bins.at(ndx) = std::make_pair(ibinRMS, ibinOffset);
  }

  return StatusCode::SUCCESS;
}
