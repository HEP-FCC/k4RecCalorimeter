#include "NoiseCaloCellsFromFileBaseTool.h"
#include "RecCaloCommon/k4RecCalorimeter_check.h"

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
  K4RECCALORIMETER_CHECK( m_geoSvc.retrieve() );
  K4RECCALORIMETER_CHECK( m_cellPositionsTool.retrieve() );
  K4RECCALORIMETER_CHECK( m_randSvc = service<IRndmGenSvc> ("RndmGenSvc", false) );
  K4RECCALORIMETER_CHECK( m_gauss.initialize(m_randSvc, Rndm::Gauss(0., 1.)) );

  // open and check file, read the histograms with noise constants
  K4RECCALORIMETER_CHECK( initNoiseFromFile() );

  // Get decoder and get index of layer/wheel field .. should put in try catch block
  m_decoder = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder();
  m_index_activeField = m_decoder->index(m_activeFieldName);

  debug() << "Filter noise threshold: " << m_filterThreshold << "*sigma" << endmsg;

  K4RECCALORIMETER_CHECK( AlgTool::initialize() );

  return StatusCode::SUCCESS;
}

template <class C>
void NoiseCaloCellsFromFileBaseTool::addRandomCellNoiseT(C& aCells) const {
  for (auto& p : aCells) {
    p.second += getNoiseOffsetPerCell(p.first);
    p.second += (getNoiseRMSPerCell(p.first) * m_gauss.shoot());
  }
}

void NoiseCaloCellsFromFileBaseTool::addRandomCellNoise(std::unordered_map<CellID, double>& aCells) const {
  addRandomCellNoiseT(aCells);
}

void NoiseCaloCellsFromFileBaseTool::addRandomCellNoise(std::vector<std::pair<CellID, double> >& aCells) const {
  addRandomCellNoiseT(aCells);
}

template <typename C>
void NoiseCaloCellsFromFileBaseTool::filterCellNoiseT(C& aCells) const {
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

void NoiseCaloCellsFromFileBaseTool::filterCellNoise(std::unordered_map<CellID, double>& aCells) const {
  filterCellNoiseT (aCells);
}

void NoiseCaloCellsFromFileBaseTool::filterCellNoise(std::vector<std::pair<CellID, double> >& aCells) const {
  filterCellNoiseT (aCells);
}

StatusCode NoiseCaloCellsFromFileBaseTool::readHistograms(TFile* noiseFile) {
  std::string elecNoiseRMSLayerHistoName, pileupNoiseRMSLayerHistoName;
  std::string elecNoiseOffsetLayerHistoName, pileupNoiseOffsetLayerHistoName;
  debug() << "Retrieving histograms" << endmsg;
  // Read the histograms with electronics noise and pileup from the file
  for (unsigned i = 0; i < m_numHistograms; i++) {
    elecNoiseRMSLayerHistoName = m_elecNoiseRMSHistoName + std::to_string(i + 1);
    debug() << "Getting histogram with a name " << elecNoiseRMSLayerHistoName << endmsg;
    m_histoElecNoiseRMS.push_back(dynamic_cast<TH1*>(noiseFile->Get(elecNoiseRMSLayerHistoName.c_str())));
    if (m_histoElecNoiseRMS.at(i)->GetNbinsX() < 1) {
      error() << "Histogram  " << elecNoiseRMSLayerHistoName
              << " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
      return StatusCode::FAILURE;
    }
    m_histoElecNoiseRMS.at(i)->SetDirectory(0); // prevent deletion when file is closed
    if (m_setNoiseOffset) {
      elecNoiseOffsetLayerHistoName = m_elecNoiseOffsetHistoName + std::to_string(i + 1);
      debug() << "Getting histogram with a name " << elecNoiseOffsetLayerHistoName << endmsg;
      m_histoElecNoiseOffset.push_back(dynamic_cast<TH1*>(noiseFile->Get(elecNoiseOffsetLayerHistoName.c_str())));
      if (m_histoElecNoiseOffset.at(i)->GetNbinsX() < 1) {
        error() << "Histogram  " << elecNoiseOffsetLayerHistoName
                << " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
        return StatusCode::FAILURE;
      }
      m_histoElecNoiseOffset.at(i)->SetDirectory(0); // prevent deletion when file is closed
    }
    if (m_addPileup) {
      pileupNoiseRMSLayerHistoName = m_pileupNoiseRMSHistoName + std::to_string(i + 1);
      debug() << "Getting histogram with a name " << pileupNoiseRMSLayerHistoName << endmsg;
      m_histoPileupNoiseRMS.push_back(dynamic_cast<TH1*>(noiseFile->Get(pileupNoiseRMSLayerHistoName.c_str())));
      if (m_histoPileupNoiseRMS.at(i)->GetNbinsX() < 1) {
        error() << "Histogram  " << pileupNoiseRMSLayerHistoName
                << " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
        return StatusCode::FAILURE;
      }
      m_histoPileupNoiseRMS.at(i)->SetDirectory(0); // prevent deletion when file is closed
      if (m_setNoiseOffset) {
        pileupNoiseOffsetLayerHistoName = m_pileupNoiseOffsetHistoName + std::to_string(i + 1);
        debug() << "Getting histogram with a name " << pileupNoiseOffsetLayerHistoName << endmsg;
        m_histoPileupNoiseOffset.push_back(
            dynamic_cast<TH1*>(noiseFile->Get(pileupNoiseOffsetLayerHistoName.c_str())));
        if (m_histoPileupNoiseOffset.at(i)->GetNbinsX() < 1) {
          error() << "Histogram  " << pileupNoiseOffsetLayerHistoName
                  << " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
          return StatusCode::FAILURE;
        }
        m_histoPileupNoiseOffset.at(i)->SetDirectory(0); // prevent deletion when file is closed
      }
    }
  }
  return StatusCode::SUCCESS;
}

StatusCode NoiseCaloCellsFromFileBaseTool::initNoiseFromFile() {
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
  if (!noiseFile || noiseFile->IsZombie()) {
    error() << "Unable to open the file with the noise values!" << endmsg;
    error() << "File path: " << m_noiseFileName.value() << endmsg;
    return StatusCode::FAILURE;
  } else {
    info() << "Using the following file with noise values: " << m_noiseFileName.value() << endmsg;
  }

  K4RECCALORIMETER_CHECK(readHistograms(noiseFile.get()));

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

unsigned NoiseCaloCellsFromFileBaseTool::getIndexHistogram(CellID aCellId) const {
  return m_decoder->get(aCellId, m_index_activeField);
}

double NoiseCaloCellsFromFileBaseTool::getNoiseRMSPerCell(CellID aCellId) const {

  double elecNoiseRMS = 0.;
  double pileupNoiseRMS = 0.;

  if (m_histoElecNoiseRMS.size() != 0) {
    // retrieve index of histogram
    unsigned indexHistogram = getIndexHistogram(aCellId);
    // retrieve bin in histogram
    int ibin = getBin(aCellId);

    if (indexHistogram < m_histoElecNoiseRMS.size()) {
      elecNoiseRMS = m_histoElecNoiseRMS.at(indexHistogram)->GetBinContent(ibin);
      if (m_addPileup) {
        pileupNoiseRMS = m_histoPileupNoiseRMS.at(indexHistogram)->GetBinContent(ibin);
      }
    } else {
      error()
          << "Index of histogram is larger than the number of noise histograms that are available !!!! Using the last one for all histograms outside the range."
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
      warning() << "Zero noise RMS: cell ID " << aCellId
              << endmsg;
  }

  return totalNoiseRMS;
}

double NoiseCaloCellsFromFileBaseTool::getNoiseOffsetPerCell(CellID aCellId) const {

  if (!m_setNoiseOffset)
    return 0.;

  double elecNoiseOffset = 0.;
  double pileupNoiseOffset = 0.;

  if (m_histoElecNoiseOffset.size() != 0) {
    // retrieve index of histogram
    unsigned indexHistogram = getIndexHistogram(aCellId);
    // retrieve bin in histogram
    int ibin = getBin(aCellId);

    if (indexHistogram < m_histoElecNoiseOffset.size()) {
      elecNoiseOffset = m_histoElecNoiseOffset.at(indexHistogram)->GetBinContent(ibin);
      if (m_addPileup) {
        pileupNoiseOffset = m_histoPileupNoiseOffset.at(indexHistogram)->GetBinContent(ibin);
      }
    } else {
      error()
          << "Index of histogram is larger than the number of noise histograms that are available !!!! Using the last one for all histograms outside the range."
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


std::pair<double, double>
NoiseCaloCellsFromFileBaseTool::getNoisePerCell(CellID aCellId) const
{
  return std::make_pair (getNoiseRMSPerCell(aCellId),
                         getNoiseOffsetPerCell(aCellId));
}
