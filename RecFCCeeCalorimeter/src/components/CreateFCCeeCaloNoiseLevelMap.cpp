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

  if (!m_hcalBarrelNoiseTool.retrieve()) {
    error() << "Unable to retrieve the HCAL noise tool!!!" << endmsg;
    return StatusCode::FAILURE;
  }
  
  std::unordered_map<uint64_t, std::pair<double,double>> map;

  ///////////////////////////////////////
  /// SEGMENTED MODULE-THETA VOLUMES  ///
  ///////////////////////////////////////
  
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
      error() << "There is no Module-Theta segmentation!!!!" << endmsg;
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
      if(m_fieldNamesSegmented[iSys] == "system" &&
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
  } // end of SEGMENTED MODULE-THETA VOLUMES 

  //////////////////////////////////
  ///      NESTED VOLUMES        ///
  //////////////////////////////////

  for (uint iSys = 0; iSys < m_readoutNamesNested.size(); iSys++) {
    // Sanity check
    if (m_activeFieldNamesNested.size() != 3) {
      error() << "Property activeFieldNamesNested requires 3 names." << endmsg;
      return StatusCode::FAILURE;
    }
    // Check if readouts exist
    info() << "Readout: " << m_readoutNamesNested[iSys] << endmsg;
    if (m_geoSvc->getDetector()->readouts().find(m_readoutNamesNested[iSys]) == m_geoSvc->getDetector()->readouts().end()) {
      error() << "Readout <<" << m_readoutNamesNested[iSys] << ">> does not exist." << endmsg;
      return StatusCode::FAILURE;
    }
    auto decoder = m_geoSvc->getDetector()->readout(m_readoutNamesNested[iSys]).idSpec().decoder();
    // Get VolumeID
    dd4hep::DDSegmentation::CellID volumeId = 0;
    decoder->set(volumeId, m_fieldNameNested, m_fieldValuesNested[iSys]);
    // Get the total number of given hierarchy of active volumes
    auto highestVol = gGeoManager->GetTopVolume();
    std::vector<unsigned int> numVolumes;
    numVolumes.reserve(m_activeVolumeNamesNested.size());
    for (const auto& volName : m_activeVolumeNamesNested) {
      numVolumes.push_back(det::utils::countPlacedVolumes(highestVol, volName));
      info() << "Number of active volumes named " << volName << " is " << numVolumes.back() << endmsg;
      if (numVolumes.back() == 0) {
        error() << "Volume name " << volName << " not found! Check naming in detector description." << endmsg;
        return StatusCode::FAILURE;
      }
    }
    // First sort to figure out which volume is inside which one
    std::vector<std::pair<std::string, uint>> numVolumesMap;
    for (unsigned int it = 0; it < m_activeVolumeNamesNested.size(); it++) {
      // names of volumes (m_activeVolumeNamesNested) not needed anymore, only corresponding bitfield names are used
      // (m_activeFieldNamesNested)
      numVolumesMap.push_back(std::pair<std::string, uint>(m_activeFieldNamesNested[it], numVolumes[it]));
    }
    std::sort(
        numVolumesMap.begin(), numVolumesMap.end(),
        [](std::pair<std::string, uint> vol1, std::pair<std::string, uint> vol2) { return vol1.second < vol2.second; });
    // now recompute how many volumes exist within the larger volume
    for (unsigned int it = numVolumesMap.size() - 1; it > 0; it--) {
      if (numVolumesMap[it].second % numVolumesMap[it - 1].second != 0) {
        error() << "Given volumes are not nested in each other!" << endmsg;
        return StatusCode::FAILURE;
      }
      numVolumesMap[it].second /= numVolumesMap[it - 1].second;
    }
    // Debug calculation of total number of cells
    if (msgLevel() <= MSG::DEBUG) {
      unsigned int checkTotal = 1;
      for (const auto& vol : numVolumesMap) {
        debug() << "Number of active volumes named " << vol.first << " is " << vol.second << endmsg;
        checkTotal *= vol.second;
      }
      debug() << "Total number of cells ( " << numVolumesMap.back().first << " ) is " << checkTotal << endmsg;
    }
    // make sure the ordering above is as in property activeFieldNamesNested
    std::map<std::string, uint> activeVolumesNumbersNested;
    for (const auto& name : m_activeFieldNamesNested) {
      for (const auto& number : numVolumesMap) {
        if (name == number.first) {
          activeVolumesNumbersNested.insert(std::make_pair(number.first, number.second));
        }
      }
    }

    // Loop over all cells in the calorimeter and retrieve existing cellIDs
    // Loop over active layers
    std::vector<std::pair<int, int>> extrema;
    extrema.push_back(std::make_pair(0, activeVolumesNumbersNested.find(m_activeFieldNamesNested[0])->second - 1));
    extrema.push_back(std::make_pair(0, activeVolumesNumbersNested.find(m_activeFieldNamesNested[1])->second - 1));
    extrema.push_back(std::make_pair(0, activeVolumesNumbersNested.find(m_activeFieldNamesNested[2])->second - 1));
    for (unsigned int ilayer = 0; ilayer < activeVolumesNumbersNested.find(m_activeFieldNamesNested[0])->second;
         ilayer++) {
      for (unsigned int iphi = 0; iphi < activeVolumesNumbersNested.find(m_activeFieldNamesNested[1])->second; iphi++) {
        for (unsigned int iz = 0; iz < activeVolumesNumbersNested.find(m_activeFieldNamesNested[2])->second; iz++) {

	  dd4hep::DDSegmentation::CellID cID = volumeId;
          decoder->set(cID, m_activeFieldNamesNested[0], ilayer);
	  decoder->set(cID, m_activeFieldNamesNested[1], iphi);
	  decoder->set(cID, m_activeFieldNamesNested[2], iz);
	  
	  double noise = m_hcalBarrelNoiseTool->getNoiseConstantPerCell(cID);
	  double noiseOffset = m_hcalBarrelNoiseTool->getNoiseOffsetPerCell(cID);
	  
	  map.insert( std::pair<uint64_t, std::pair<double, double> >(cID, std::make_pair(noise, noiseOffset) ) );
        }
      }
    }
  } // end of NESTED VOLUMES

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
