#include "CreateFCCeeEndcapTurbineNoiseLevelMap.h"

#include "DD4hep/Detector.h"
#include "detectorCommon/DetUtils_k4geo.h"
#include "detectorSegmentations/FCCSWEndcapTurbine_k4geo.h"
#include "k4Interface/IGeoSvc.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"

DECLARE_COMPONENT(CreateFCCeeEndcapTurbineNoiseLevelMap)

CreateFCCeeEndcapTurbineNoiseLevelMap::CreateFCCeeEndcapTurbineNoiseLevelMap(const std::string& aName, ISvcLocator* aSL)
    : base_class(aName, aSL) {
  declareProperty( "outputFileName", m_outputFileName, "Name of the output file");
}

CreateFCCeeEndcapTurbineNoiseLevelMap::~CreateFCCeeEndcapTurbineNoiseLevelMap() {}

StatusCode CreateFCCeeEndcapTurbineNoiseLevelMap::initialize() {
  // Initialize necessary Gaudi components
  if (Service::initialize().isFailure()) {
    error() << "Unable to initialize Service()" << endmsg;
    return StatusCode::FAILURE;
  }
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service. "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }
  std::unordered_map<uint64_t, std::pair<double,double>> map;

  //////////////////////////////////
  /// SEGMENTED VOLUME           ///
  //////////////////////////////////

  for (uint iSys = 0; iSys < m_readoutNamesSegmented.size(); iSys++) {
    // Check if readouts exist
    info() << "Readout: " << m_readoutNamesSegmented[iSys] << endmsg;
    if (m_geoSvc->getDetector()->readouts().find(m_readoutNamesSegmented[iSys]) == m_geoSvc->getDetector()->readouts().end()) {
      error() << "Readout <<" << m_readoutNamesSegmented[iSys] << ">> does not exist." << endmsg;
      return StatusCode::FAILURE;
    }
    // get segmentation
    dd4hep::DDSegmentation::FCCSWEndcapTurbine_k4geo* segmentation;
    segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWEndcapTurbine_k4geo*>(
										m_geoSvc->getDetector()->readout(m_readoutNamesSegmented[iSys]).segmentation().segmentation());
    if (segmentation == nullptr) {
      error() << "There is no turbine segmentation!!!!" << endmsg;
      return StatusCode::FAILURE;
    }
    
    unsigned layerOffset[3];
    layerOffset[0] = 0;
    layerOffset[1] = segmentation->numCellsRhoCalib(0)*segmentation->numCellsZCalib(0);
    layerOffset[2] = layerOffset[1] + segmentation->numCellsRhoCalib(1)*segmentation->numCellsZCalib(1);
    
    auto decoder = m_geoSvc->getDetector()->readout(m_readoutNamesSegmented[iSys]).idSpec().decoder();
    // Loop over all cells in the calorimeter and retrieve existing cellIDs
    // Loop over active layers
    std::vector<std::pair<int, int>> extrema;
    extrema.push_back(std::make_pair(0, 2));  // wheel
    extrema.push_back(std::make_pair(0, 0)); // modules (set per wheel)
    extrema.push_back(std::make_pair(0, 0)); // rho (set per wheel)
    extrema.push_back(std::make_pair(0, 0)); // z (set per wheel)

    for (int iSide = -1; iSide < 2; iSide+=2) {
      for (unsigned int iWheel = 0; iWheel < 3; iWheel++)
	{
	  dd4hep::DDSegmentation::CellID volumeId = 0;
	  // Get VolumeID
	  // for ECAL OK (volume extends along z)
	  (*decoder)[m_fieldNamesSegmented[iSys]].set(volumeId, m_fieldValuesSegmented[iSys]);
	  (*decoder)["side"].set(volumeId, iSide);
	  (*decoder)["wheel"].set(volumeId, iWheel);
	  (*decoder)["module"].set(volumeId, 0);
	  (*decoder)["rho"].set(volumeId, 0);
	  (*decoder)["z"].set(volumeId, 0);
	  
	  int numModules = segmentation->nModules(iWheel);
	  int numCellsRho = segmentation->numCellsRho(iWheel);
	  int numCellsZ = segmentation->numCellsZ(iWheel);
	  int numCellsRhoCalib = segmentation->numCellsRhoCalib(iWheel);
	  int numCellsZCalib = segmentation->numCellsZCalib(iWheel);
	  
	  // extrema 1: 0, ID of last module
	  extrema[1] = std::make_pair(0, numModules - 1);
	  // extrema[2]: 0, ID of last rho cell
	  extrema[2] = std::make_pair(0, numCellsRho - 1);
	  // extrema[3]: 0, ID of last z cell
	  extrema[3] = std::make_pair(0, numCellsZ - 1);
	  
	  // Loop over segmentation cells
	  for (int imodule = 0; imodule < numModules; imodule++) {
	    for (int irho = 0; irho < numCellsRho; irho++) {
	      for (int iz = 0; iz < numCellsZ; iz++) {
		dd4hep::DDSegmentation::CellID cellId = volumeId;
		decoder->set(cellId, "module", imodule);
		decoder->set(cellId, "rho", irho);
		decoder->set(cellId, "z", iz);
		unsigned iLayerZ = iz/(numCellsZ/numCellsZCalib);
		unsigned iLayerRho = irho/(numCellsRho/numCellsRhoCalib);
		unsigned iLayer = layerOffset[iWheel] + iLayerRho*numCellsZCalib + iLayerZ;
		decoder->set(cellId, "layer", iLayer);

		uint64_t id = cellId;
		if (iSide==1 && iWheel == 2 && imodule == 113 && irho == 19 && iz == 1) {
		  debug() << "in test cell, iLayer = " << iLayer << " and cell ID = " << id << endmsg;
		}
		double noiseRMS = 1e-12;  // dummy small non-zero value
		double noiseOffset = 0.;
		map.insert( std::pair<uint64_t, std::pair<double, double> >(id, std::make_pair(noiseRMS, noiseOffset) ) );
	      }
	    }
	  }
	}
    }
  }
  // Check if output directory exists
  std::string outDirPath = gSystem->DirName(m_outputFileName.c_str());
  if (!gSystem->OpenDirectory(outDirPath.c_str())) {
    error() << "Output directory \"" << outDirPath
            << "\" does not exists! Please create it." << endmsg;
    return StatusCode::FAILURE;
  }

  std::unique_ptr<TFile> outFile(TFile::Open(m_outputFileName.c_str(), "RECREATE"));
  outFile->cd();
  TTree tree("noisyCells", "Tree with map of noise per cell");
  uint64_t saveCellId;
  double saveNoiseLevel;
  double saveNoiseOffset;
  tree.Branch("cellId", &saveCellId, "cellId/l");
  tree.Branch("noiseLevel", &saveNoiseLevel);
  tree.Branch("noiseOffset", &saveNoiseOffset);
  for (const auto& item : map) {
    saveCellId = item.first;
    saveNoiseLevel = item.second.first;
    saveNoiseOffset = item.second.second;
    tree.Fill();
  }
  outFile->Write();
  outFile->Close();

  return StatusCode::SUCCESS;
}

StatusCode CreateFCCeeEndcapTurbineNoiseLevelMap::finalize() { return Service::finalize(); }
