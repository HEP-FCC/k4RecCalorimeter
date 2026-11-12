#include "ReadNoiseFromFileTool.h"
#include "DD4hep/Detector.h"

#include "detectorSegmentations/FCCSWGridPhiEta_k4geo.h"
#include "k4Interface/IGeoSvc.h"

#include "TH1F.h"

DECLARE_COMPONENT(ReadNoiseFromFileTool)

StatusCode ReadNoiseFromFileTool::initBinning(NoiseData& data, const k4::recCalo::ICaloIndexer& indexer) const {
  /// PhiEta segmentation
  const auto* segmentation = dynamic_cast<const dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo*>(
      m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());
  if (segmentation == nullptr) {
    error() << "There is no phi-eta segmentation!!!!" << endmsg;
    return StatusCode::FAILURE;
  }

  data.m_bins.resize(indexer.cellIDs().size());

  auto deltaEta = [](const TH1* h) -> std::pair<int, double> {
    if (!h)
      return std::make_pair(0, 0.);
    int Nbins = h->GetNbinsX();
    double delta = (h->GetBinLowEdge(Nbins) + h->GetBinWidth(Nbins) - h->GetBinLowEdge(1)) / Nbins;
    return std::make_pair(Nbins, delta);
  };
  auto [NbinsRMS, deltaEtaRMS] = deltaEta(data.m_histoElecNoiseRMS.at(0).get());
  auto [NbinsOffset, deltaEtaOffset] = deltaEta(m_setNoiseOffset ? data.m_histoElecNoiseOffset.at(0).get() : nullptr);

  for (uint64_t id : indexer.cellIDs()) {
    unsigned ndx = indexer.index(id);
    double cellEta = segmentation->eta(id);

    {
      int ibin = floor(fabs(cellEta) / deltaEtaRMS) + 1;
      if (ibin > NbinsRMS) {
        error() << "eta outside range of the RMS histograms! Cell eta: " << cellEta
                << " Nbins in histogram: " << NbinsRMS << endmsg;
        ibin = NbinsRMS;
      }
      data.m_bins.at(ndx).first = ibin;
    }
    {
      int ibin = floor(fabs(cellEta) / deltaEtaOffset) + 1;
      if (ibin > NbinsOffset) {
        error() << "eta outside range of the offset histograms! Cell eta: " << cellEta
                << " Nbins in histogram: " << NbinsOffset << endmsg;
        ibin = NbinsOffset;
      }
      data.m_bins.at(ndx).second = ibin;
    }
  }

  return StatusCode::SUCCESS;
}
