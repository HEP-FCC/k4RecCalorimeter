#include "CreateFCCeeCaloNoiseLevelMap.h"

#include "DD4hep/Detector.h"
#include "detectorCommon/DetUtils_k4geo.h"
#include "k4Interface/IGeoSvc.h"
#include "detectorSegmentations/FCCSWGridModuleThetaMerged_k4geo.h"

#include "TFile.h"
#include "TTree.h"

DECLARE_COMPONENT(CreateFCCeeCaloNoiseLevelMap)

CreateFCCeeCaloNoiseLevelMap::CreateFCCeeCaloNoiseLevelMap(const std::string& aName, ISvcLocator* aSL)
    : base_class(aName, aSL) {
  declareProperty("ECalBarrelNoiseTool", m_ecalBarrelNoiseTool, "Handle for the cells noise tool of Barrel ECal");
  declareProperty("HCalBarrelNoiseTool", m_hcalBarrelNoiseTool, "Handle for the cells noise tool of Barrel HCal");
  declareProperty( "outputFileName", m_outputFileName, "Name of the output file");
}

CreateFCCeeCaloNoiseLevelMap::~CreateFCCeeCaloNoiseLevelMap() {}

StatusCode CreateFCCeeCaloNoiseLevelMap::initialize() {
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

  if (!m_ecalBarrelNoiseTool.retrieve()) {
    error() << "Unable to retrieve the ECAL noise tool!!!" << endmsg;
    return StatusCode::FAILURE;
  }

  if (!m_hcalBarrelNoiseTool.empty()){
    if (!m_hcalBarrelNoiseTool.retrieve()) {
      error() << "Unable to retrieve the HCAL noise tool!!!" << endmsg;
      return StatusCode::FAILURE;
    }
  }
  
  std::unordered_map<uint64_t, std::pair<double,double>> map;

  for (uint iSys = 0; iSys < m_readoutNamesSegmented.size(); iSys++) {
    // Check if readouts exist
    info() << "Readout: " << m_readoutNamesSegmented[iSys] << endmsg;
    if (m_geoSvc->getDetector()->readouts().find(m_readoutNamesSegmented[iSys]) == m_geoSvc->getDetector()->readouts().end()) {
      error() << "Readout <<" << m_readoutNamesSegmented[iSys] << ">> does not exist." << endmsg;
      return StatusCode::FAILURE;
    }
    // get segmentation
    dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo* segmentation;
    segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo*>(
        m_geoSvc->getDetector()->readout(m_readoutNamesSegmented[iSys]).segmentation().segmentation());
    if (segmentation == nullptr) {
      error() << "There is no Module-theta segmentation!!!!" << endmsg;
      return StatusCode::FAILURE;
    }

    info() << "FCCSWGridModuleThetaMerged: size in Theta " << segmentation->gridSizeTheta() << " , bins in Module " << segmentation->nModules()
           << endmsg;
    info() << "FCCSWGridModuleThetaMerged: offset in Theta " << segmentation->offsetTheta() << endmsg;

    auto decoder = m_geoSvc->getDetector()->readout(m_readoutNamesSegmented[iSys]).idSpec().decoder();
    // Loop over all cells in the calorimeter and retrieve existing cellIDs
    // Loop over active layers
    std::vector<std::pair<int, int>> extrema;
    extrema.push_back(std::make_pair(0, m_activeVolumesNumbersSegmented[iSys] - 1));
    extrema.push_back(std::make_pair(0, 0));
    extrema.push_back(std::make_pair(0, 0));
    for (unsigned int ilayer = 0; ilayer < m_activeVolumesNumbersSegmented[iSys]; ilayer++) {
      dd4hep::DDSegmentation::CellID volumeId = 0;
      // Get VolumeID
      (*decoder)[m_fieldNamesSegmented[iSys]].set(volumeId, m_fieldValuesSegmented[iSys]);
      (*decoder)[m_activeFieldNamesSegmented[iSys]].set(volumeId, ilayer);
      (*decoder)["theta"].set(volumeId, 0);
      (*decoder)["module"].set(volumeId, 0);
      // Get number of segmentation cells within the active volume, for given layer      
      auto numCells = det::utils::numberOfCells(volumeId, *segmentation);
      // Range of module ID
      extrema[1] = std::make_pair(0, (numCells[0] - 1)*segmentation->mergedModules(ilayer));
      // Range and min of theta ID
      // for HCAL use alternative calculation
      if (m_fieldNamesSegmented[iSys] == "system" &&
	  m_fieldValuesSegmented[iSys] == m_hcalBarrelSysId){
	uint cellsTheta = ceil(( 2*m_activeVolumesTheta[ilayer] - segmentation->gridSizeTheta() ) / 2 / segmentation->gridSizeTheta()) * 2 + 1; //ceil( 2*m_activeVolumesRadii[ilayer] / segmentation->gridSizeTheta());
	uint minThetaID = int(floor(( - m_activeVolumesTheta[ilayer] + 0.5 * segmentation->gridSizeTheta() - segmentation->offsetTheta()) / segmentation->gridSizeTheta()));
	numCells[1]=cellsTheta;
	numCells[2]=minThetaID;
      }
      extrema[2] = std::make_pair(numCells[2], numCells[2] + (numCells[1] - 1)*segmentation->mergedThetaCells(ilayer));
      debug() << "Layer: " << ilayer << endmsg;
      debug() << "Number of segmentation cells in (module, theta): " << numCells << endmsg;
      // Loop over segmentation cells
      for (unsigned int imodule = 0; imodule < numCells[0]; imodule++) {
        for (unsigned int itheta = 0; itheta < numCells[1]; itheta++) {
	  dd4hep::DDSegmentation::CellID cellId = volumeId;
	  decoder->set(cellId, "module", imodule*segmentation->mergedModules(ilayer));
	  decoder->set(cellId, "theta", numCells[2] + itheta*segmentation->mergedThetaCells(ilayer));  // start from the minimum existing theta cell in this layer
	  uint64_t id = cellId;
	  double noise = 0.;
	  double noiseOffset = 0.;
	  if (m_fieldValuesSegmented[iSys] == m_hcalBarrelSysId){
	    noise = m_hcalBarrelNoiseTool->getNoiseConstantPerCell(id);
	    noiseOffset = m_hcalBarrelNoiseTool->getNoiseOffsetPerCell(id);
	  } else if (m_fieldValuesSegmented[iSys] == m_ecalBarrelSysId){
	    noise = m_ecalBarrelNoiseTool->getNoiseConstantPerCell(id);
            noiseOffset = m_ecalBarrelNoiseTool->getNoiseOffsetPerCell(id);
	  }
          map.insert( std::pair<uint64_t, std::pair<double, double> >(id, std::make_pair(noise, noiseOffset) ) );
        }
      }
    }
  }

  std::unique_ptr<TFile> file(TFile::Open(m_outputFileName.c_str(), "RECREATE"));
  file->cd();
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
  file->Write();
  file->Close();

  return StatusCode::SUCCESS;
}

StatusCode CreateFCCeeCaloNoiseLevelMap::finalize() { return Service::finalize(); }
