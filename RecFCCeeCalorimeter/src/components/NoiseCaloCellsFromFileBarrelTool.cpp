#include "NoiseCaloCellsFromFileBarrelTool.h"

DECLARE_COMPONENT(NoiseCaloCellsFromFileBarrelTool)

unsigned NoiseCaloCellsFromFileBarrelTool::getBin(CellID aCellId) const {
  double cellTheta = m_cellPositionsTool->xyzPosition(aCellId).Theta();

  // All histograms have same binning, all bins with same size
  // Using the histogram in the first layer to get the bin size
  int Nbins = m_histoElecNoiseRMS.at(0)->GetNbinsX();
  int ibin = m_histoElecNoiseRMS.at(0)->FindFixBin(cellTheta);

  if (ibin > Nbins) {
    error() << "theta outside range of the histograms! Cell theta: " << cellTheta << " Nbins in histogram: " << Nbins
            << endmsg;
    ibin = Nbins;
  }

  return ibin;
}
