#include "ReadNoiseVsThetaFromFileTool.h"

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

DECLARE_COMPONENT(ReadNoiseVsThetaFromFileTool)

ReadNoiseVsThetaFromFileTool::ReadNoiseVsThetaFromFileTool(const std::string& type, const std::string& name, const IInterface* parent)
    : GaudiTool(type, name, parent) {
  declareInterface<INoiseConstTool>(this);
  declareProperty("cellPositionsTool", m_cellPositionsTool, "Handle for tool to retrieve cell positions");
}

StatusCode ReadNoiseVsThetaFromFileTool::initialize() {
  // Get GeoSvc
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service. "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }

  // Check if cell position tool available if m_useSeg==false; if tool not
  // available, try using segmentation instead
  if (!m_useSeg){
    if (!m_cellPositionsTool.retrieve()) {
      error() << "Unable to retrieve cell positions tool, and useSegmentation is false." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // Get segmentation
  if (m_useSeg) {
    m_segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo*>(
      m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());
    if (m_segmentation == nullptr) {
      error() << "There is no module-theta segmentation." << endmsg;
      return StatusCode::FAILURE;
    }
    else
      info() << "Found module-theta segmentation." << endmsg;
  }

  // open and check file, read the histograms with noise constants
  if (ReadNoiseVsThetaFromFileTool::initNoiseFromFile().isFailure()) {
    error() << "Couldn't open file with noise constants!!!" << endmsg;
    return StatusCode::FAILURE;
  }

  // Take readout bitfield decoder from GeoSvc
  m_decoder = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder();
  if (m_decoder == nullptr) {
    error() << "Cannot create decore for readout " << m_readoutName << endmsg;
    return StatusCode::FAILURE;
  }

  StatusCode sc = GaudiTool::initialize();
  if (sc.isFailure()) return sc;

  return sc;
}

StatusCode ReadNoiseVsThetaFromFileTool::finalize() {
  StatusCode sc = GaudiTool::finalize();
  return sc;
}

StatusCode ReadNoiseVsThetaFromFileTool::initNoiseFromFile() {
  // check if file exists
  if (m_noiseFileName.empty()) {
    error() << "Name of the file with the noise values not provided!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (gSystem->AccessPathName(m_noiseFileName.value().c_str())) {
    error() << "Provided file with the noise values not found!" << endmsg;
    error() << "File path: " << m_noiseFileName.value() << endmsg;
    return StatusCode::FAILURE;
  }
  std::unique_ptr<TFile> inFile(TFile::Open(m_noiseFileName.value().c_str(), "READ"));
  if (inFile->IsZombie()) {
    error() << "Couldn't open the file with the noise values!" << endmsg;
    error() << "File path: " << m_noiseFileName.value() << endmsg;
    return StatusCode::FAILURE;
  } else {
    info() << "Opening the following file with the noise values: "
           << m_noiseFileName << endmsg;
  }

  std::string elecNoiseLayerHistoName, pileupLayerHistoName;
  std::string elecNoiseOffsetLayerHistoName, pileupOffsetLayerHistoName;
  // Read the histograms with electronics noise and pileup from the file
  for (unsigned i = 0; i < m_numRadialLayers; i++) {
    elecNoiseLayerHistoName = m_elecNoiseHistoName + std::to_string(i + 1);
    debug() << "Getting histogram with a name " << elecNoiseLayerHistoName << endmsg;
    m_histoElecNoiseConst.push_back(*dynamic_cast<TH1F*>(inFile->Get(elecNoiseLayerHistoName.c_str())));
    if (m_histoElecNoiseConst.at(i).GetNbinsX() < 1) {
      error() << "Histogram  " << elecNoiseLayerHistoName
              << " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_setNoiseOffset){
      elecNoiseOffsetLayerHistoName = m_elecNoiseOffsetHistoName + std::to_string(i + 1);
      debug() << "Getting histogram with a name " << elecNoiseOffsetLayerHistoName << endmsg;
      m_histoElecNoiseOffset.push_back(*dynamic_cast<TH1F*>(inFile->Get(elecNoiseOffsetLayerHistoName.c_str())));
      if (m_histoElecNoiseOffset.at(i).GetNbinsX() < 1) {
        error() << "Histogram  " << elecNoiseOffsetLayerHistoName
                << " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
        return StatusCode::FAILURE;
      }
    }
    if (m_addPileup) {
      pileupLayerHistoName = m_pileupHistoName + std::to_string(i + 1);
      debug() << "Getting histogram with a name " << pileupLayerHistoName << endmsg;
      m_histoPileupConst.push_back(*dynamic_cast<TH1F*>(inFile->Get(pileupLayerHistoName.c_str())));
      if (m_histoPileupConst.at(i).GetNbinsX() < 1) {
        error() << "Histogram  " << pileupLayerHistoName
                << " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
        return StatusCode::FAILURE;
      }
      if (m_setNoiseOffset == true){
        pileupOffsetLayerHistoName = m_pileupOffsetHistoName + std::to_string(i + 1);
        debug() << "Getting histogram with a name " << pileupOffsetLayerHistoName << endmsg;
        m_histoPileupOffset.push_back(*dynamic_cast<TH1F*>(inFile->Get(pileupOffsetLayerHistoName.c_str())));
        if (m_histoElecNoiseOffset.at(i).GetNbinsX() < 1) {
          error() << "Histogram  " << pileupOffsetLayerHistoName
                  << " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
          return StatusCode::FAILURE;
        }
      }
    }
  }
  // Check if we have same number of histograms (all layers) for pileup and electronics noise
  if (m_histoElecNoiseConst.size() == 0 ) {
    error() << "No histograms with noise found!!!!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_addPileup) {
    if (m_histoElecNoiseConst.size() != m_histoPileupConst.size()) {
      error() << "Missing histograms! Different number of histograms for electronics noise and pileup!!!!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  return StatusCode::SUCCESS;
}

double ReadNoiseVsThetaFromFileTool::getNoiseConstantPerCell(uint64_t aCellId) {

  double elecNoise = 0.;
  double pileupNoise = 0.;

  // Get cell coordinates: eta/theta and radial layer
  dd4hep::DDSegmentation::CellID cID = aCellId;
  double cellTheta;
  // checked that for baseline theta-module merged segmentation
  // the two approaches give identical theta.
  // however, code based on positioning tool is more general
  // since it can be run for any segmentation class without the need
  // to change the interface of this tool..
  if (m_useSeg)
    cellTheta = m_segmentation->theta(aCellId);
  else
    cellTheta = m_cellPositionsTool->xyzPosition(aCellId).Theta();

  unsigned cellLayer = m_decoder->get(cID, m_activeFieldName);

  // All histograms have same binning, all bins with same size
  // Using the histogram in the first layer to get the bin size
  unsigned index = 0;
  if (m_histoElecNoiseConst.size() != 0) {
    int Nbins = m_histoElecNoiseConst.at(index).GetNbinsX();
    int ibin = m_histoElecNoiseConst.at(index).FindBin(cellTheta);
    if (ibin > Nbins) {
      error() << "theta outside range of the histograms! Cell theta: " << cellTheta << " Nbins in histogram: " << Nbins << endmsg;
      ibin = Nbins;
    }
    // Check that there are not more layers than the constants are provided for
    if (cellLayer < m_histoElecNoiseConst.size()) {
      elecNoise = m_histoElecNoiseConst.at(cellLayer).GetBinContent(ibin);
      if (m_addPileup) {
        pileupNoise = m_histoPileupConst.at(cellLayer).GetBinContent(ibin);
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
  double totalNoise = sqrt(pow(elecNoise, 2) + pow(pileupNoise, 2)) * m_scaleFactor;

  if (totalNoise < 1e-6) {
    warning() << "Zero noise: cell theta " << cellTheta << " layer "
              << cellLayer << " noise " << totalNoise << endmsg;
  }

  return totalNoise;
}

double ReadNoiseVsThetaFromFileTool::getNoiseOffsetPerCell(uint64_t aCellId) {

  if (!m_setNoiseOffset)
    return 0.;
  else {
    double elecNoise = 0.;
    double pileupNoise = 0.;

  // Get cell coordinates: eta and radial layer
  dd4hep::DDSegmentation::CellID cID = aCellId;
  double cellTheta;
  if (m_useSeg)
    cellTheta = m_segmentation->theta(aCellId);
  else
    cellTheta = m_cellPositionsTool->xyzPosition(aCellId).Theta();
  unsigned cellLayer = m_decoder->get(cID, m_activeFieldName);

  // All histograms have same binning, all bins with same size
  // Using the histogram in the first layer to get the bin size
  unsigned index = 0;
  if (m_histoElecNoiseOffset.size() != 0) {
    int Nbins = m_histoElecNoiseOffset.at(index).GetNbinsX();
    int ibin = m_histoElecNoiseOffset.at(index).FindBin(cellTheta);
    if (ibin > Nbins) {
      error() << "theta outside range of the histograms! Cell theta: "
              << cellTheta << " Nbins in histogram: " << Nbins << endmsg;
      ibin = Nbins;
    }

    // Check that there are not more layers than the constants are provided for
    if (cellLayer < m_histoElecNoiseOffset.size()) {
      elecNoise = m_histoElecNoiseOffset.at(cellLayer).GetBinContent(ibin);
      if (m_addPileup) {
        pileupNoise = m_histoPileupOffset.at(cellLayer).GetBinContent(ibin);
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
  double totalNoise = sqrt(pow(elecNoise, 2) + pow(pileupNoise, 2)) * m_scaleFactor;

  if (totalNoise < 1e-6) {
    warning() << "Zero noise: cell theta " << cellTheta << " layer " << cellLayer << " noise " << totalNoise << endmsg;
  }

  return totalNoise;
  }
}
