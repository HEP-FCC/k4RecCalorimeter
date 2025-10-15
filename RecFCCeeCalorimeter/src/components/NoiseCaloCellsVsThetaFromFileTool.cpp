#include "NoiseCaloCellsVsThetaFromFileTool.h"
#include "RecCaloCommon/k4RecCalorimeter_check.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"

// k4FWCore
#include "k4Interface/IGeoSvc.h"

// DD4hep
#include "DD4hep/Detector.h"

// ROOT
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TSystem.h"

DECLARE_COMPONENT(NoiseCaloCellsVsThetaFromFileTool)

StatusCode NoiseCaloCellsVsThetaFromFileTool::initialize() {
  K4RECCALORIMETER_CHECK( m_geoSvc.retrieve() );
  K4RECCALORIMETER_CHECK( m_cellPositionsTool.retrieve() );
  K4RECCALORIMETER_CHECK( m_randSvc = service<IRndmGenSvc> ("RndmGenSvc", false) );
  K4RECCALORIMETER_CHECK( m_gauss.initialize(m_randSvc, Rndm::Gauss(0., 1.)) );

  // open and check file, read the histograms with noise constants
  K4RECCALORIMETER_CHECK( initNoiseFromFile() );

  // Get decoder and index of layer field
  // put in try... block
  m_decoder = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder();
  m_index_activeField = m_decoder->index(m_activeFieldName);

  debug() << "Filter noise threshold: " << m_filterThreshold << "*sigma" << endmsg;

  K4RECCALORIMETER_CHECK( AlgTool::initialize() );

  return StatusCode::SUCCESS;
}

template <class C>
void NoiseCaloCellsVsThetaFromFileTool::addRandomCellNoiseT(C& aCells) const {
  for (auto& p : aCells) {
    p.second += getNoiseOffsetPerCell(p.first);
    p.second += (getNoiseRMSPerCell(p.first) * m_gauss.shoot());
  }
}

void NoiseCaloCellsVsThetaFromFileTool::addRandomCellNoise(std::unordered_map<uint64_t, double>& aCells) const {
  addRandomCellNoiseT(aCells);
}

void NoiseCaloCellsVsThetaFromFileTool::addRandomCellNoise(std::vector<std::pair<uint64_t, double> >& aCells) const {
  addRandomCellNoiseT (aCells);
}

template <typename C>
void NoiseCaloCellsVsThetaFromFileTool::filterCellNoiseT(C& aCells) const {
  // Erase a cell if it has energy below a threshold from the vector
  if (m_useAbsInFilter) {
    std::erase_if(aCells,
                  [&](auto& p) { return std::abs(p.second-getNoiseOffsetPerCell(p.first)) < m_filterThreshold * getNoiseRMSPerCell(p.first); });
  }
  else {
    std::erase_if(aCells,
                  [&](auto& p) { return p.second < getNoiseOffsetPerCell(p.first) + m_filterThreshold * getNoiseRMSPerCell(p.first); });
  }
}

void NoiseCaloCellsVsThetaFromFileTool::filterCellNoise(std::unordered_map<uint64_t, double>& aCells) const {
  filterCellNoiseT (aCells);
}

void NoiseCaloCellsVsThetaFromFileTool::filterCellNoise(std::vector<std::pair<uint64_t, double> >& aCells) const {
  filterCellNoiseT (aCells);
}

StatusCode NoiseCaloCellsVsThetaFromFileTool::initNoiseFromFile() {
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
    info() << "Using the following file with noise values: " << m_noiseFileName.value() << endmsg;
  }

  std::string elecNoiseRMSLayerHistoName, pileupNoiseRMSLayerHistoName;
  std::string elecNoiseOffsetLayerHistoName, pileupNoiseOffsetLayerHistoName;
  // Read the histograms with electronics noise and pileup from the file
  for (unsigned i = 0; i < m_numRadialLayers; i++) {
    elecNoiseRMSLayerHistoName = m_elecNoiseRMSHistoName + std::to_string(i + 1);
    debug() << "Getting histogram with a name " << elecNoiseRMSLayerHistoName << endmsg;
    m_histoElecNoiseRMS.push_back(*dynamic_cast<TH1F*>(noiseFile->Get(elecNoiseRMSLayerHistoName.c_str())));
    if (m_histoElecNoiseRMS.at(i).GetNbinsX() < 1) {
      error() << "Histogram  " << elecNoiseRMSLayerHistoName
              << " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_setNoiseOffset) {
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
      pileupNoiseRMSLayerHistoName = m_pileupNoiseRMSHistoName + std::to_string(i + 1);
      debug() << "Getting histogram with a name " << pileupNoiseRMSLayerHistoName << endmsg;
      m_histoPileupNoiseRMS.push_back(*dynamic_cast<TH1F*>(noiseFile->Get(pileupNoiseRMSLayerHistoName.c_str())));
      if (m_histoPileupNoiseRMS.at(i).GetNbinsX() < 1) {
        error() << "Histogram  " << pileupNoiseRMSLayerHistoName
                << " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
        return StatusCode::FAILURE;
      }
      if (m_setNoiseOffset) {
        pileupNoiseOffsetLayerHistoName = m_pileupNoiseOffsetHistoName + std::to_string(i + 1);
        debug() << "Getting histogram with a name " << pileupNoiseOffsetLayerHistoName << endmsg;
        m_histoPileupNoiseOffset.push_back(
            *dynamic_cast<TH1F*>(noiseFile->Get(pileupNoiseOffsetLayerHistoName.c_str())));
        if (m_histoElecNoiseOffset.at(i).GetNbinsX() < 1) {
          error() << "Histogram  " << pileupNoiseOffsetLayerHistoName
                  << " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
          return StatusCode::FAILURE;
        }
      }
    }
  }
  noiseFile->Close();

  // Check if we have same number of histograms (all layers) for pileup and electronics noise
  if (m_histoElecNoiseRMS.size() == 0) {
    error() << "No histograms with noise found!!!!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_addPileup) {
    if (m_histoElecNoiseRMS.size() != m_histoPileupNoiseRMS.size()) {
      error() << "Missing histograms! Different number of histograms for electronics noise RMS and pileup noise RMS!!!!"
              << endmsg;
      return StatusCode::FAILURE;
    }
  }
  if (m_setNoiseOffset) {
    if (m_histoElecNoiseOffset.size() != m_histoElecNoiseRMS.size()) {
      error() << "Missing histograms! Different number of histograms for electronics noise RMS and offset!!!!"
              << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_addPileup) {
      if (m_histoPileupNoiseOffset.size() != m_histoElecNoiseRMS.size()) {
        error() << "Missing histograms! Different number of histograms for electronics noise RMS and pileup noise "
                   "offset!!!!"
                << endmsg;
        return StatusCode::FAILURE;
      }
    }
  }

  return StatusCode::SUCCESS;
}

double NoiseCaloCellsVsThetaFromFileTool::getNoiseRMSPerCell(uint64_t aCellId) const {

  double elecNoiseRMS = 0.;
  double pileupNoiseRMS = 0.;

  double cellTheta = m_cellPositionsTool->xyzPosition(aCellId).Theta();
  unsigned cellLayer = m_decoder->get(aCellId, m_index_activeField);

  // All histograms have same binning, all bins with same size
  // Using the histogram in the first layer to get the bin size
  if (m_histoElecNoiseRMS.size() != 0) {
    unsigned index = 0;
    int ibin = m_histoElecNoiseRMS.at(index).FindFixBin(cellTheta);
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
  double totalNoiseRMS = 0;
  if (m_addPileup) {
    totalNoiseRMS = sqrt(elecNoiseRMS * elecNoiseRMS + pileupNoiseRMS * pileupNoiseRMS) * m_scaleFactor;
  } else { // avoid useless math operations if no pileup
    totalNoiseRMS = elecNoiseRMS * m_scaleFactor;
  }

  if (totalNoiseRMS < 1e-6) {
    warning() << "Zero noise RMS: cell theta " << cellTheta << " layer " << cellLayer << " noise RMS " << totalNoiseRMS
              << endmsg;
  }

  return totalNoiseRMS;
}

double NoiseCaloCellsVsThetaFromFileTool::getNoiseOffsetPerCell(uint64_t aCellId) const {

  if (!m_setNoiseOffset)
    return 0.;

  double elecNoiseOffset = 0.;
  double pileupNoiseOffset = 0.;

  // Get cell coordinates: theta and radial layer
  double cellTheta = m_cellPositionsTool->xyzPosition(aCellId).Theta();
  unsigned cellLayer = m_decoder->get(aCellId, m_activeFieldName);

  // All histograms have same binning, all bins with same size
  // Using the histogram in the first layer to get the bin size
  if (m_histoElecNoiseOffset.size() != 0) {
    unsigned index = 0;
    int Nbins = m_histoElecNoiseOffset.at(index).GetNbinsX();
    int ibin = m_histoElecNoiseOffset.at(index).FindFixBin(cellTheta);
    if (ibin > Nbins) {
      error() << "theta outside range of the histograms! Cell theta: " << cellTheta << " Nbins in histogram: " << Nbins
              << endmsg;
      ibin = Nbins;
    }

    // Check that there are not more layers than the constants are provided for
    if (cellLayer < m_histoElecNoiseOffset.size()) {
      elecNoiseOffset = m_histoElecNoiseOffset.at(cellLayer).GetBinContent(ibin);
      if (m_addPileup) {
        pileupNoiseOffset = m_histoPileupNoiseOffset.at(cellLayer).GetBinContent(ibin);
      }
    } else {
      error()
          << "More radial layers than we have noise for!!!! Using the last layer for all histograms outside the range."
          << endmsg;
    }
  } else {
    error() << "No histograms with noise offset!!!!! " << endmsg;
  }

  // Total noise offset: electronics noise + pileup (linear or quadratic sum??)
  // double totalNoiseOffset = sqrt(pow(elecNoiseOffset, 2) + pow(pileupNoiseOffset, 2)) * m_scaleFactor;
  double totalNoiseOffset = (elecNoiseOffset + pileupNoiseOffset) * m_scaleFactor;

  return totalNoiseOffset;
}
