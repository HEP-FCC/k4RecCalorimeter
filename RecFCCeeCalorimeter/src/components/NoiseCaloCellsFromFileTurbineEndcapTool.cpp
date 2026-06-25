#include "NoiseCaloCellsFromFileTurbineEndcapTool.h"

DECLARE_COMPONENT(NoiseCaloCellsFromFileTurbineEndcapTool)

unsigned NoiseCaloCellsFromFileTurbineEndcapTool::getBin(CellID aCellId) const {
  unsigned iHist = m_decoder->get(aCellId, m_index_activeField);
  unsigned iRho = m_decoder->get(aCellId, "rho") + 1;
  unsigned iZ = m_decoder->get(aCellId, "z") + 1;

  unsigned NbinsZ = m_histoElecNoiseRMS.at(iHist)->GetNbinsX();
  unsigned NbinsRho = m_histoElecNoiseRMS.at(iHist)->GetNbinsY();

  unsigned ibin = 0;
  if (iRho > NbinsRho || iZ > NbinsZ) {
    error() << "bins outside range of the histograms! Bins: " << iRho << "," << iZ
            << " , Nbins in histogram: " << NbinsRho << "," << NbinsZ << endmsg;
  } else {
    ibin = m_histoElecNoiseRMS.at(iHist)->GetBin(iZ, iRho);
  }

  return ibin;
}
