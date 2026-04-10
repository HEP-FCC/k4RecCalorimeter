#include "NoiseCaloCellsTurbineEndcapFromFileTool.h"
#include "RecCaloCommon/k4RecCalorimeter_check.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"

// k4FWCore
#include "k4Interface/IGeoSvc.h"

// DD4hep
#include "DD4hep/Detector.h"

// ROOT
#include "TFile.h"
#include "TH2F.h"
#include "TMath.h"
#include "TSystem.h"

DECLARE_COMPONENT(NoiseCaloCellsTurbineEndcapFromFileTool)

StatusCode NoiseCaloCellsTurbineEndcapFromFileTool::initialize() {
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
void NoiseCaloCellsTurbineEndcapFromFileTool::addRandomCellNoiseT(C& aCells) const {
  for (auto& p : aCells) {
    p.second += getNoiseOffsetPerCell(p.first);
    p.second += (getNoiseRMSPerCell(p.first) * m_gauss.shoot());
  }
}

void NoiseCaloCellsTurbineEndcapFromFileTool::addRandomCellNoise(std::unordered_map<uint64_t, double>& aCells) const {
  addRandomCellNoiseT(aCells);
}

void NoiseCaloCellsTurbineEndcapFromFileTool::addRandomCellNoise(std::vector<std::pair<uint64_t, double> >& aCells) const {
  addRandomCellNoiseT (aCells);
}

template <typename C>
void NoiseCaloCellsTurbineEndcapFromFileTool::filterCellNoiseT(C& aCells) const {
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

void NoiseCaloCellsTurbineEndcapFromFileTool::filterCellNoise(std::unordered_map<uint64_t, double>& aCells) const {
  filterCellNoiseT (aCells);
}

void NoiseCaloCellsTurbineEndcapFromFileTool::filterCellNoise(std::vector<std::pair<uint64_t, double> >& aCells) const {
  filterCellNoiseT (aCells);
}

StatusCode NoiseCaloCellsTurbineEndcapFromFileTool::initNoiseFromFile() {
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
  
  std::string elecNoiseRMSWheelHistoName, pileupNoiseRMSWheelHistoName;
  std::string elecNoiseOffsetWheelHistoName, pileupNoiseOffsetWheelHistoName;
  // Read the histograms with electronics noise and pileup from the file
  for (unsigned iWheel = 0; iWheel < m_numWheels; iWheel++) {
    elecNoiseRMSWheelHistoName = m_elecNoiseRMSHistoName +std::to_string(iWheel);
    debug() << "Getting histogram with a name " <<  elecNoiseRMSWheelHistoName << endmsg;
    m_histoElecNoiseRMS.push_back(* dynamic_cast<TH2F*>(noiseFile->Get(elecNoiseRMSWheelHistoName.c_str())));
    if (m_histoElecNoiseRMS.at(iWheel).GetNbinsX() < 1) {
      error() << "Histogram  " << elecNoiseRMSWheelHistoName 
	      << " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_setNoiseOffset) {
      elecNoiseOffsetWheelHistoName = m_elecNoiseOffsetHistoName +std::to_string(iWheel);
      debug() << "Getting histogram with a name " << elecNoiseOffsetWheelHistoName  << endmsg;
      m_histoElecNoiseOffset.push_back(*dynamic_cast<TH2F*>(noiseFile->Get(elecNoiseOffsetWheelHistoName.c_str())));
      if (m_histoElecNoiseOffset.at(iWheel).GetNbinsX() < 1) {
	error() << "Histogram  " <<  elecNoiseOffsetWheelHistoName
		<< " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
	return StatusCode::FAILURE;
      }
    }
    if (m_addPileup) {
      pileupNoiseRMSWheelHistoName = m_pileupNoiseRMSHistoName +std::to_string(iWheel);
      debug() << "Getting histogram with a name " <<  pileupNoiseRMSWheelHistoName << endmsg;
      m_histoPileupNoiseRMS.push_back(* dynamic_cast<TH2F*>(noiseFile->Get(pileupNoiseRMSWheelHistoName.c_str())));
      if (m_histoPileupNoiseRMS.at(iWheel).GetNbinsX() < 1) {
	error() << "Histogram  " << pileupNoiseRMSWheelHistoName 
		<< " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
	return StatusCode::FAILURE;
      }
      if (m_setNoiseOffset) {
	pileupNoiseOffsetWheelHistoName = m_pileupNoiseOffsetHistoName +std::to_string(iWheel);
	debug() << "Getting histogram with a name " << pileupNoiseOffsetWheelHistoName  << endmsg;
	m_histoPileupNoiseOffset.push_back(*dynamic_cast<TH2F*>(noiseFile->Get(pileupNoiseOffsetWheelHistoName.c_str())));
	if (m_histoPileupNoiseOffset.at(iWheel).GetNbinsX() < 1) {
	  error() << "Histogram  " <<  pileupNoiseOffsetWheelHistoName
		  << " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
	  return StatusCode::FAILURE;
	}
      }
    }
  }
  noiseFile->Close();

  return StatusCode::SUCCESS;
}

double NoiseCaloCellsTurbineEndcapFromFileTool::getNoiseRMSPerCell(uint64_t aCellId) const {

  double elecNoiseRMS = 0.;
  double pileupNoiseRMS = 0.;

  unsigned iWheel = m_decoder->get(aCellId, "wheel");
  unsigned iRho = m_decoder->get(aCellId, "rho");
  unsigned iZ = m_decoder->get(aCellId, "z");

  // All histograms have same binning, all bins with same size
  // Using the histogram in the first layer to get the bin size
  if (m_histoElecNoiseRMS.size() >  0) {
     // Check that there are not more wheels than the constants are provided for
    if (iWheel < m_histoElecNoiseRMS.size()) {      
      elecNoiseRMS = m_histoElecNoiseRMS.at(iWheel).GetBinContent(iZ+1, iRho+1);
      debug() << "Here with noise RMS = " << elecNoiseRMS << endmsg;
      if (m_addPileup) {
	pileupNoiseRMS = m_histoPileupNoiseRMS.at(iWheel).GetBinContent(iZ+1, iRho+1);
      }
    } else {
      error() << "More wheels than we have noise for!!!!! " << endmsg;
    }
  
  } else {
    error() << "No histograms with RMS noise constants!!!!! " << endmsg;
  }

  // Total noise: electronics noise + pileup
  double totalNoiseRMS = 0;
  if (m_addPileup) {
    totalNoiseRMS = sqrt(elecNoiseRMS * elecNoiseRMS + pileupNoiseRMS * pileupNoiseRMS) * m_scaleFactor;
  } else { // avoid useless math operations if no pileup
    totalNoiseRMS = elecNoiseRMS * m_scaleFactor;
  }
  if (totalNoiseRMS < 1e-6) {
    warning() << "Zero noise RMS:  iZ, iRho " << iZ << "," << iRho  << " noise RMS " << totalNoiseRMS
              << endmsg;
  }
  
  return totalNoiseRMS;
}

double NoiseCaloCellsTurbineEndcapFromFileTool::getNoiseOffsetPerCell(uint64_t aCellId) const {

  if (!m_setNoiseOffset)
    return 0.;

  double elecNoiseOffset = 0.;
  double pileupNoiseOffset = 0.;
  
  unsigned iWheel = m_decoder->get(aCellId, "wheel");
  unsigned iRho = m_decoder->get(aCellId, "rho");
  unsigned iZ = m_decoder->get(aCellId, "z");
  
  // All histograms have same binning, all bins with same size
  // Using the histogram in the first layer to get the bin size
  if (m_histoElecNoiseOffset.size() >  0) {
    // Check that there are not more wheels than the constants are provided for
    if (iWheel < m_histoElecNoiseOffset.size()) {
      elecNoiseOffset = m_histoElecNoiseOffset.at(iWheel).GetBinContent(iZ+1, iRho+1);
      if (m_addPileup) {
	pileupNoiseOffset = m_histoPileupNoiseOffset.at(iWheel).GetBinContent(iZ+1, iRho+1);
      }
    } else {
      error() << "More wheels than we have noise for!!!!! " << endmsg;
    }
    
  } else {
    error() << "No histograms with Offset noise constants!!!!! " << endmsg;
  }
  
  // Total noise: electronics noise + pileup
  double totalNoiseOffset = 0;
  if (m_addPileup) {
    totalNoiseOffset = sqrt(elecNoiseOffset * elecNoiseOffset + pileupNoiseOffset * pileupNoiseOffset) * m_scaleFactor;
  } else { // avoid useless math operations if no pileup
    totalNoiseOffset = elecNoiseOffset * m_scaleFactor;
  }
  
  if (totalNoiseOffset < 1e-6) {
    warning() << "Zero noise Offset:  iZ, iRho " << iZ << "," << iRho << " noise Offset " << totalNoiseOffset
              << endmsg;
  }
  
  return totalNoiseOffset;
}

std::pair<double, double> NoiseCaloCellsTurbineEndcapFromFileTool::getNoisePerCell(uint64_t aCellId) const
{
  return std::make_pair (getNoiseRMSPerCell(aCellId),
                         getNoiseOffsetPerCell(aCellId));
}
