#include "NoiseCaloCellsFromFileTool.h"
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

DECLARE_COMPONENT(NoiseCaloCellsFromFileTool)

StatusCode NoiseCaloCellsFromFileTool::initialize() {

  K4RECCALORIMETER_CHECK( m_geoSvc.retrieve() );
  K4RECCALORIMETER_CHECK( m_randSvc = service<IRndmGenSvc> ("RndmGenSvc", false) );
  K4RECCALORIMETER_CHECK( m_gauss.initialize(m_randSvc, Rndm::Gauss(0., 1.)) );

  // open and check file, read the histograms with noise constants
  K4RECCALORIMETER_CHECK( initNoiseFromFile() );

  // Check if cell position tool available
  if (!m_cellPositionsTool.retrieve() and !m_useSeg) {
    info() << "Unable to retrieve cell positions tool, try eta-phi segmentation." << endmsg;
    // Get PhiEta segmentation
    m_segmentationPhiEta = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo*>(
        m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());
    if (m_segmentationPhiEta == nullptr) {
      error() << "There is no phi-eta segmentation." << endmsg;
      return StatusCode::FAILURE;
    } else
      info() << "Found phi-eta segmentation." << endmsg;
  }
  // Get PhiEta segmentation
  m_segmentationPhiEta = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo*>(
      m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());
  if (m_segmentationPhiEta == nullptr) {
    m_segmentationMulti = dynamic_cast<dd4hep::DDSegmentation::MultiSegmentation*>(
        m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());
    if (m_segmentationMulti == nullptr) {
      error() << "There is no phi-eta or multi- segmentation for the readout " << m_readoutName << " defined."
              << endmsg;
      return StatusCode::FAILURE;
    } else {
      // check if multisegmentation contains only phi-eta sub-segmentations
      const dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo* subsegmentation = nullptr;
      for (const auto& subSegm : m_segmentationMulti->subSegmentations()) {
        subsegmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo*>(subSegm.segmentation);
        if (subsegmentation == nullptr) {
          error() << "At least one of the sub-segmentations in MultiSegmentation named " << m_readoutName
                  << " is not a phi-eta grid." << endmsg;
          return StatusCode::FAILURE;
        }
      }
    }
  }

  m_decoder = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder();
  m_index_activeField = m_decoder->index(m_activeFieldName);

  debug() << "Filter noise threshold: " << m_filterThreshold << "*sigma" << endmsg;

  K4RECCALORIMETER_CHECK( AlgTool::initialize() );

  return StatusCode::SUCCESS;
}

template <typename C>
void NoiseCaloCellsFromFileTool::addRandomCellNoiseT (C& aCells) const
{
  for (auto& p : aCells) {
    p.second += getNoiseRMSPerCell(p.first) * m_gauss.shoot();
  }
}

void NoiseCaloCellsFromFileTool::addRandomCellNoise(std::unordered_map<uint64_t, double>& aCells) const
{
  addRandomCellNoiseT(aCells);
}

void NoiseCaloCellsFromFileTool::addRandomCellNoise(std::vector<std::pair<uint64_t, double> >& aCells) const
{
  addRandomCellNoiseT (aCells);
}

template <typename C>
void NoiseCaloCellsFromFileTool::filterCellNoiseT (C& aCells) const
{
  // Erase a cell if it has energy bellow a threshold from the vector
  if (m_useAbsInFilter) {
    std::erase_if (aCells,
                   [&](auto& p) { return std::abs(p.second) < m_filterThreshold * getNoiseRMSPerCell(p.first); });
  }
  else {
    std::erase_if (aCells,
                   [&](auto& p) { return p.second < m_filterThreshold * getNoiseRMSPerCell(p.first); });
  }
}

void NoiseCaloCellsFromFileTool::filterCellNoise(std::unordered_map<uint64_t, double>& aCells) const
{
  filterCellNoiseT(aCells);
}

void NoiseCaloCellsFromFileTool::filterCellNoise(std::vector<std::pair<uint64_t, double> >& aCells) const
{
  filterCellNoiseT(aCells);
}

StatusCode NoiseCaloCellsFromFileTool::initNoiseFromFile() {
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

  std::string elecNoiseLayerHistoName, pileupLayerHistoName;
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
    if (m_addPileup) {
      pileupLayerHistoName = m_pileupHistoName + std::to_string(i + 1);
      debug() << "Getting histogram with a name " << pileupLayerHistoName << endmsg;
      m_histoPileupNoiseRMS.push_back(*dynamic_cast<TH1F*>(noiseFile->Get(pileupLayerHistoName.c_str())));
      if (m_histoPileupNoiseRMS.at(i).GetNbinsX() < 1) {
        error() << "Histogram  " << pileupLayerHistoName
                << " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
        return StatusCode::FAILURE;
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
      error() << "Missing histograms! Different number of histograms for electronics noise and pileup!!!!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  return StatusCode::SUCCESS;
}

double NoiseCaloCellsFromFileTool::getNoiseRMSPerCell(uint64_t aCellId) const {
  const dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo* segmentation = m_segmentationPhiEta;
  if (segmentation == nullptr) {
    segmentation = dynamic_cast<const dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo*>(
        &m_segmentationMulti->subsegmentation(aCellId));
  }

  double elecNoiseRMS = 0.;
  double pileupNoiseRMS = 0.;

  double cellEta;
  if (m_useSeg)
    cellEta = m_segmentationPhiEta->eta(aCellId);
  else
    cellEta = m_cellPositionsTool->xyzPosition(aCellId).Eta();
  unsigned cellLayer = m_decoder->get(aCellId, m_index_activeField);

  // All histograms have same binning, all bins with same size
  // Using the histogram in the first layer to get the bin size
  unsigned index = 0;
  if (m_histoElecNoiseRMS.size() != 0) {
    int ibin = m_histoElecNoiseRMS.at(index).FindFixBin(fabs(cellEta));
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
    warning() << "Zero noise RMS: cell eta " << cellEta << " layer " << cellLayer << " noise " << totalNoiseRMS
              << endmsg;
  }

  return totalNoiseRMS;
}
