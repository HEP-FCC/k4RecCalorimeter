#include "ReadNoiseFromFileTool.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"

// k4FWCore
#include "k4Interface/IGeoSvc.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DD4hep/Readout.h"
#include "DDSegmentation/Segmentation.h"

// ROOT
#include "TSystem.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"

DECLARE_COMPONENT(ReadNoiseFromFileTool)

ReadNoiseFromFileTool::ReadNoiseFromFileTool(const std::string& type, const std::string& name, const IInterface* parent)
    : GaudiTool(type, name, parent) {
  declareInterface<INoiseConstTool>(this);
}

StatusCode ReadNoiseFromFileTool::initialize() {
  // Get GeoSvc
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service. "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }
  m_segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo*>(m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());
  if (m_segmentation == nullptr) {
    error() << "There is no phi-eta segmentation!!!!" << endmsg;
    return StatusCode::FAILURE;
  }

  // open and check file, read the histograms with noise constants
  if (ReadNoiseFromFileTool::initNoiseFromFile().isFailure()) {
    error() << "Couldn't open file with noise constants!!!" << endmsg;
    return StatusCode::FAILURE;
  }

  // Take readout bitfield decoder from GeoSvc
  m_decoder = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder();

  StatusCode sc = GaudiTool::initialize();
  if (sc.isFailure()) return sc;

  return sc;
}

StatusCode ReadNoiseFromFileTool::finalize() {
  StatusCode sc = GaudiTool::finalize();
  return sc;
}

StatusCode ReadNoiseFromFileTool::initNoiseFromFile() {
  // Check if file exists
  if (m_noiseFileName.empty()) {
    error() << "Name of the file with the noise values not provided!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (gSystem->AccessPathName(m_noiseFileName.value().c_str())) {
    error() << "Provided file with the noise values not found!" << endmsg;
    error() << "File path: " << m_noiseFileName.value() << endmsg;
    return StatusCode::FAILURE;
  }
  std::unique_ptr<TFile> noiseFile(TFile::Open(m_noiseFileName.value().c_str(), "READ"));
  if (noiseFile->IsZombie()) {
    error() << "Unable to open the file with the noise values!" << endmsg;
    error() << "File path: " << m_noiseFileName.value() << endmsg;
    return StatusCode::FAILURE;
  } else {
    info() << "Using the following file with the noise values: "
           << m_noiseFileName.value() << endmsg;
  }

  std::string elecNoiseLayerHistoName, pileupLayerHistoName;
  std::string elecNoiseOffsetLayerHistoName, pileupOffsetLayerHistoName;
  // Read the histograms with electronics noise and pileup from the file
  for (unsigned i = 0; i < m_numRadialLayers; i++) {
    elecNoiseLayerHistoName = m_elecNoiseHistoName + std::to_string(i + 1);
    debug() << "Getting histogram with a name " << elecNoiseLayerHistoName << endmsg;
    m_histoElecNoiseRMS.push_back(*dynamic_cast<TH1F*>(noiseFile->Get(elecNoiseLayerHistoName.c_str())));
    if (m_histoElecNoiseRMS.at(i).GetNbinsX() < 1) {
      error() << "Histogram  " << elecNoiseLayerHistoName
              << " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_setNoiseOffset){
      elecNoiseOffsetLayerHistoName = m_elecNoiseOffsetHistoName + std::to_string(i + 1);
      debug() << "Getting histogram with a name " << elecNoiseOffsetLayerHistoName << endmsg;
      m_histoElecNoiseOffset.push_back(*dynamic_cast<TH1F*>(noiseFile->Get(elecNoiseOffsetLayerHistoName.c_str())));
      if (m_histoElecNoiseOffset.at(i).GetNbinsX() < 1) {
	error() << "Histogram  " << elecNoiseOffsetLayerHistoName
		<< " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
	return StatusCode::FAILURE;
      }
    }
    if (m_addPileup) {
      pileupLayerHistoName = m_pileupHistoName + std::to_string(i + 1);
      debug() << "Getting histogram with a name " << pileupLayerHistoName << endmsg;
      m_histoPileupNoiseRMS.push_back(*dynamic_cast<TH1F*>(noiseFile->Get(pileupLayerHistoName.c_str())));
      if (m_histoPileupNoiseRMS.at(i).GetNbinsX() < 1) {
        error() << "Histogram  " << pileupLayerHistoName
                << " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
        return StatusCode::FAILURE;
      }
      if (m_setNoiseOffset == true){
	pileupOffsetLayerHistoName = m_pileupOffsetHistoName + std::to_string(i + 1);
	debug() << "Getting histogram with a name " << pileupOffsetLayerHistoName << endmsg;
	m_histoPileupOffset.push_back(*dynamic_cast<TH1F*>(noiseFile->Get(pileupOffsetLayerHistoName.c_str())));
	if (m_histoElecNoiseOffset.at(i).GetNbinsX() < 1) {
	  error() << "Histogram  " << pileupOffsetLayerHistoName
		  << " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
	  return StatusCode::FAILURE;
	}
      }
    }
  }
  // Check if we have same number of histograms (all layers) for pileup and electronics noise
  if (m_histoElecNoiseRMS.size() == 0) {
    error() << "No histograms with noise found!!!!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_addPileup) {
    if (m_histoElecNoiseRMS.size() != m_histoPileupNoiseRMS.size()) {
      error() << "Missing histograms! Different number of histograms for electronics noise and pileup!!!!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  return StatusCode::SUCCESS;
}

double ReadNoiseFromFileTool::getNoiseRMSPerCell(uint64_t aCellId) {

  double elecNoiseRMS = 0.;
  double pileupNoiseRMS = 0.;

  // Get cell coordinates: eta and radial layer
  dd4hep::DDSegmentation::CellID cID = aCellId;
  double cellEta = m_segmentation->eta(cID);

  unsigned cellLayer = m_decoder->get(cID, m_activeFieldName);

  // All histograms have same binning, all bins with same size
  // Using the histogram in the first layer to get the bin size
  unsigned index = 0;
  if (m_histoElecNoiseRMS.size() != 0) {
    int Nbins = m_histoElecNoiseRMS.at(index).GetNbinsX();
    double deltaEtaBin =
        (m_histoElecNoiseRMS.at(index).GetBinLowEdge(Nbins) + m_histoElecNoiseRMS.at(index).GetBinWidth(Nbins) -
         m_histoElecNoiseRMS.at(index).GetBinLowEdge(1)) /
        Nbins;
    // find the eta bin for the cell
    int ibin = floor(fabs(cellEta) / deltaEtaBin) + 1;
    if (ibin > Nbins) {
      error() << "eta outside range of the histograms! Cell eta: " << cellEta << " Nbins in histogram: " << Nbins
              << endmsg;
      ibin = Nbins;
    }
    // Check that there are not more layers than the constants are provided for
    if (cellLayer < m_histoElecNoiseRMS.size()) {
      elecNoiseRMS = m_histoElecNoiseRMS.at(cellLayer).GetBinContent(ibin);
      if (m_addPileup) {
        pileupNoiseRMS = m_histoPileupNoiseRMS.at(cellLayer).GetBinContent(ibin);
      }
    } else {
      error()
          << "More radial layers than we have noise for!!!! Using the last layer for all histograms outside the range."
          << endmsg;
    }
  } else {
    error() << "No histograms with noise constants!!!!! " << endmsg;
  }

  // Total noise: electronics noise + pileup
  double totalNoiseRMS = sqrt(elecNoiseRMS*elecNoiseRMS + pileupNoiseRMS*pileupNoiseRMS) * m_scaleFactor;

  if (totalNoiseRMS < 1e-6) {
    warning() << "Zero noise: cell eta " << cellEta << " layer " << cellLayer << " noise " << totalNoiseRMS << endmsg;
  }

  return totalNoiseRMS;
}

double ReadNoiseFromFileTool::getNoiseOffsetPerCell(uint64_t aCellId) {

  if (!m_setNoiseOffset) return 0.;

  double elecNoiseOffset = 0.;
  double pileupNoiseOffset = 0.;

  // Get cell coordinates: eta and radial layer
  dd4hep::DDSegmentation::CellID cID = aCellId;
  double cellEta = m_segmentation->eta(cID);
  unsigned cellLayer = m_decoder->get(cID, m_activeFieldName);

  // All histograms have same binning, all bins with same size
  // Using the histogram in the first layer to get the bin size
  unsigned index = 0;
  if (m_histoElecNoiseOffset.size() != 0) {
    int Nbins = m_histoElecNoiseOffset.at(index).GetNbinsX();
    double deltaEtaBin =
        (m_histoElecNoiseOffset.at(index).GetBinLowEdge(Nbins) + m_histoElecNoiseOffset.at(index).GetBinWidth(Nbins) -
         m_histoElecNoiseOffset.at(index).GetBinLowEdge(1)) /
        Nbins;
    // find the eta bin for the cell
    int ibin = floor(fabs(cellEta) / deltaEtaBin) + 1;
    if (ibin > Nbins) {
      error() << "eta outside range of the histograms! Cell eta: " << cellEta << " Nbins in histogram: " << Nbins
              << endmsg;
      ibin = Nbins;
    }
    // Check that there are not more layers than the constants are provided for
    if (cellLayer < m_histoElecNoiseOffset.size()) {
      elecNoiseOffset = m_histoElecNoiseOffset.at(cellLayer).GetBinContent(ibin);
      if (m_addPileup) {
        pileupNoiseOffset = m_histoPileupOffset.at(cellLayer).GetBinContent(ibin);
      }
    } else {
      error()
          << "More radial layers than we have noise for!!!! Using the last layer for all histograms outside the range."
          << endmsg;
    }
  } else {
    error() << "No histograms with noise offset!!!!! " << endmsg;
  }

  // Total noise: electronics noise + pileup
  double totalNoiseOffset = sqrt(elecNoiseOffset*elecNoiseOffset + pileupNoiseOffset*pileupNoiseOffset) * m_scaleFactor; // shouldnt the offset be summed linearly?

  // No warning: offset is usually zero or close
  //if (totalNoiseOffset < 1e-6) {
  //  warning() << "Zero noise offset: cell eta " << cellEta << " layer " << cellLayer << " noise " << totalNoiseOffset << endmsg;
  //}

  return totalNoiseOffset;
}
