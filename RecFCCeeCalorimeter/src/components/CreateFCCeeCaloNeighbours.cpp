#include "CreateFCCeeCaloNeighbours.h"

#include "DD4hep/Detector.h"
#include "DetCommon/DetUtils.h"
#include "k4Interface/IGeoSvc.h"

#include "TFile.h"
#include "TTree.h"

DECLARE_COMPONENT(CreateFCCeeCaloNeighbours)

CreateFCCeeCaloNeighbours::CreateFCCeeCaloNeighbours(const std::string& aName, ISvcLocator* aSL)
    : base_class(aName, aSL) {
  declareProperty( "outputFileName", m_outputFileName, "Name of the output file");
}

CreateFCCeeCaloNeighbours::~CreateFCCeeCaloNeighbours() {}

StatusCode CreateFCCeeCaloNeighbours::initialize() {
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
  std::unordered_map<uint64_t, std::vector<uint64_t>> map;
 
  // will be used for volume connecting
  //int eCalLastLayer;
  std::pair<int, int> extremaECalLastLayerModule;
  std::pair<int, int> extremaECalLastLayerTheta;
  //double eCalThetaOffset = 0;
  //double eCalThetaSize = 0;
  //double eCalPhiOffset = 0;
  //double eCalModuleSize = 0;
  //double hCalThetaOffset = 0;
  //double hCalThetaSize = 0;
  //double hCalPhiOffset = 0;
  //dd4hep::DDSegmentation::BitFieldCoder* decoderECalBarrel = nullptr;
  // will be used for volume connecting
  std::pair<int, int> extremaHCalFirstLayerPhi;
  std::pair<int, int> extremaHCalFirstLayerTheta;
  std::pair<int, int> extremaHCalFirstLayerZ;
  //dd4hep::DDSegmentation::BitFieldCoder* decoderHCalBarrel = nullptr;

  ////////////////////////////////////
  /// SEGMENTED THETA-MODULE VOLUMES  ///
  ////////////////////////////////////
  
  for (uint iSys = 0; iSys < m_readoutNamesSegmented.size(); iSys++) {
    // Check if readouts exist
    info() << "Readout: " << m_readoutNamesSegmented[iSys] << endmsg;
    if (m_geoSvc->getDetector()->readouts().find(m_readoutNamesSegmented[iSys]) == m_geoSvc->getDetector()->readouts().end()) {
      error() << "Readout <<" << m_readoutNamesSegmented[iSys] << ">> does not exist." << endmsg;
      return StatusCode::FAILURE;
    }
    
    // get Theta-Module Merged segmentation
    dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged* segmentation;
    segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged*>(
        m_geoSvc->getDetector()->readout(m_readoutNamesSegmented[iSys]).segmentation().segmentation());
    if (segmentation == nullptr) {
      error() << "There is no Theta-Module Merged segmentation!!!!" << endmsg;
      return StatusCode::FAILURE;
    }
    info() << "FCCSWGridModuleThetaMerged: size in Theta " << segmentation->gridSizeTheta() << " , bins in Module " << segmentation->nModules()
           << endmsg;
    info() << "FCCSWGridModuleThetaMerged: offset in Theta " << segmentation->offsetTheta() << endmsg;

    // retrieve decoders and other info needed for volume (ECal-HCal) connection
    auto decoder = m_geoSvc->getDetector()->readout(m_readoutNamesSegmented[iSys]).idSpec().decoder();
    if (m_fieldNamesSegmented[iSys] == "system" && m_fieldValuesSegmented[iSys] == 5) {
      //decoderECalBarrel = decoder;
      //eCalThetaSize = segmentation->gridSizeTheta();
      //eCalModuleSize = 2 * M_PI / segmentation->nModules();
      //eCalModuleSize = 2 * M_PI / segmentation->phiBins();
      //eCalThetaOffset = segmentation->offsetTheta();
      //eCalPhiOffset = segmentation->offsetPhi();
    }
    if (m_fieldNamesSegmented[iSys] == "system" && m_fieldValuesSegmented[iSys] == 8) {
      //decoderHCalBarrel = decoder;
      //hCalThetaSize = segmentation->gridSizeTheta();
      //hCalThetaOffset = segmentation->offsetTheta();
      //hCalPhiOffset = segmentation->offsetPhi();
    }

    // Loop over all cells in the calorimeter and retrieve existing cellIDs
    // Loop over active layers
    std::vector<std::pair<int, int>> extrema;
    // extrema[0]: min layer, n layers
    extrema.push_back(std::make_pair(0, m_activeVolumesNumbersSegmented[iSys] - 1));
    extrema.push_back(std::make_pair(0, 0));
    extrema.push_back(std::make_pair(0, 0));
    for (unsigned int ilayer = 0; ilayer < m_activeVolumesNumbersSegmented[iSys]; ilayer++) {
      dd4hep::DDSegmentation::CellID volumeId = 0;
      // Get VolumeID
      // m_fieldValuesSegmented: in .py systemValuesModuleTheta = [4]
      // m_activeFieldNamesSegmented: in .py activeFieldNamesModuleTheta = ["layer"]
      (*decoder)[m_fieldNamesSegmented[iSys]].set(volumeId, m_fieldValuesSegmented[iSys]);
      (*decoder)[m_activeFieldNamesSegmented[iSys]].set(volumeId, ilayer);
      (*decoder)["theta"].set(volumeId, 0);
      (*decoder)["module"].set(volumeId, 0);
      // Get number of segmentation cells within the active volume
      // numberOfCells: return Array of the number of cells in (module, theta) and the minimum theta ID.
      auto numCells = det::utils::numberOfCells(volumeId, *segmentation);
      // extrema 1: min module number (0), max module number
      extrema[1] = std::make_pair(0, (numCells[0] - 1)*segmentation->mergedModules(ilayer));
      // extrema[2]: min theta ID, n (merged) theta cells
      extrema[2] = std::make_pair(numCells[2], numCells[2] + (numCells[1] - 1)*segmentation->mergedThetaCells(ilayer));
	
      // for layer N-1 of ECal barrel,  will be used for volume connecting
      // should 5 be systemValuesModuleTheta instead?
      if (ilayer == (m_activeVolumesNumbersSegmented[iSys] - 1) && m_fieldNamesSegmented[iSys] == "system" &&
          m_fieldValuesSegmented[iSys] == 5) {
        //eCalLastLayer = m_activeVolumesNumbersSegmented[iSys] - 1;
        extremaECalLastLayerModule = std::make_pair(0, numCells[0] - 1);
        extremaECalLastLayerTheta = std::make_pair(numCells[2], numCells[1] + numCells[2] - 1);
      }
      else if(m_fieldNamesSegmented[iSys] == "system" &&
	      m_fieldValuesSegmented[iSys] == 8 && m_readoutNamesSegmented[iSys]=="BarHCal_Readout_phitheta"){
////
	uint cellsTheta = ceil(( 2*m_activeVolumesTheta[ilayer] - segmentation->gridSizeTheta() ) / 2 / segmentation->gridSizeTheta()) * 2 + 1; //ceil( 2*m_activeVolumesRadii[ilayer] / segmentation->gridSizeTheta());
	uint minThetaID = int(floor(( - m_activeVolumesTheta[ilayer] + 0.5 * segmentation->gridSizeTheta() - segmentation->offsetTheta()) / segmentation->gridSizeTheta()));
	numCells[1]=cellsTheta;
	numCells[2]=minThetaID;
	// for layer 0 of HCal barrel,  will be used for volume connecting
	if (ilayer == 0){	
	  extremaHCalFirstLayerPhi = std::make_pair(0, numCells[0] - 1);
	  extremaHCalFirstLayerTheta = std::make_pair(numCells[2], numCells[1] + numCells[2] - 1);	
	} 
      }
      debug() << "Layer: " << ilayer << endmsg;
      debug() << "Extrema[0]: " << extrema[0].first << " , " << extrema[0].second << endmsg;
      debug() << "Extrema[1]: " << extrema[1].first << " , " << extrema[1].second << endmsg;
      debug() << "Extrema[2]: " << extrema[2].first << " , " << extrema[2].second << endmsg;
      debug() << "Number of segmentation cells in (module,theta): " << numCells << endmsg;
      // Loop over segmentation cells
      for (int imodule = extrema[1].first; imodule <= extrema[1].second; imodule += segmentation->mergedModules(ilayer)) {
        for (int itheta = extrema[2].first; itheta <= extrema[2].second; itheta += segmentation->mergedThetaCells(ilayer)) {
	  dd4hep::DDSegmentation::CellID cellId = volumeId;
	  decoder->set(cellId, "module", imodule);
	  decoder->set(cellId, "theta", itheta);  // start from the minimum existing theta cell in this layer
	  uint64_t id = cellId;
	  map.insert(std::pair<uint64_t, std::vector<uint64_t>>(
								id,
								det::utils::neighbours_ModuleThetaMerged(
													 *segmentation,
													 *decoder,
													 {m_activeFieldNamesSegmented[iSys],
													  "module", "theta"},
													 extrema,
													 id,
													 m_includeDiagonalCells
													 )));
        }
      }
    }
    if (msgLevel() <= MSG::DEBUG) {
      std::vector<int> counter;
      counter.assign(40, 0);
      for (const auto& item : map) {
        counter[item.second.size()]++;
      }
      for (uint iCount = 0; iCount < counter.size(); iCount++) {
        if (counter[iCount] != 0) {
          info() << counter[iCount] << " cells have " << iCount << " neighbours" << endmsg;
        }
      }
    }
    info() << "total number of cells:  " <<  map.size() << endmsg;
  }

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
    // will be used for volume connecting
    if (m_fieldNameNested == "system" && m_fieldValuesNested[iSys] == 8) {
      //decoderHCalBarrel = decoder;
    }
    //hCalPhiOffset = m_hCalPhiOffset;
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
    // for layer 0 of HCal barrel
    if (m_fieldNameNested == "system" && m_fieldValuesNested[iSys] == 8) {
      extremaHCalFirstLayerPhi =
          std::make_pair(0, activeVolumesNumbersNested.find(m_activeFieldNamesNested[1])->second - 1);
      extremaHCalFirstLayerZ =
          std::make_pair(0, activeVolumesNumbersNested.find(m_activeFieldNamesNested[2])->second - 1);
    }
    for (unsigned int ilayer = 0; ilayer < activeVolumesNumbersNested.find(m_activeFieldNamesNested[0])->second;
         ilayer++) {
      for (unsigned int iphi = 0; iphi < activeVolumesNumbersNested.find(m_activeFieldNamesNested[1])->second; iphi++) {
        for (unsigned int iz = 0; iz < activeVolumesNumbersNested.find(m_activeFieldNamesNested[2])->second; iz++) {

	  dd4hep::DDSegmentation::CellID cID = volumeId;
          decoder->set(cID, m_activeFieldNamesNested[0], ilayer);
	  decoder->set(cID, m_activeFieldNamesNested[1], iphi);
	  decoder->set(cID, m_activeFieldNamesNested[2], iz);
          
          map.insert(std::pair<uint64_t, std::vector<uint64_t>>(
              cID, det::utils::neighbours(*decoder, {m_activeFieldNamesNested[0], m_activeFieldNamesNested[1],
                                                          m_activeFieldNamesNested[2]},
                                               extrema, cID, {false, true, false}, true)));
        }
      }
    }
    if (msgLevel() <= MSG::DEBUG) {
      std::vector<int> counter;
      counter.assign(40, 0);
      for (const auto& item : map) {
        counter[item.second.size()]++;
      }
      for (uint iCount = 0; iCount < counter.size(); iCount++) {
        if (counter[iCount] != 0) {
          info() << counter[iCount] << " cells have " << iCount << " neighbours" << endmsg;
        }
      }
    }
  }

  //////////////////////////////////////////////////
  ///      BARREL: connection ECAL + HCAL        ///
  /////////////////////////////////////////////////
  int count=0;
/*
  if (m_connectBarrels) {
    // first check if ECAL barrel (system==5) and HCal barrel (system==8) are configured
    if (decoderECalBarrel == nullptr || decoderHCalBarrel == nullptr) {
      error() << "ECAL barrel and/or HCal barrel are not configured correctly!" << endmsg;
      return StatusCode::FAILURE;
    }
    // print how many cells in each dimensions will be matched
    if(m_readoutNamesNested.size()!=0){
      info() << "ECAL layer " << eCalLastLayer << " is a neighbour of HCAL layer 0." << endmsg;
      info() << "ECAL phi cells " << extremaECalLastLayerModule.first << " - " << extremaECalLastLayerModule.second
	     << " will be matched to HCAL " << m_activeFieldNamesNested[1] << "(s) " << extremaHCalFirstLayerPhi.first
	     << " - " << extremaHCalFirstLayerPhi.second << endmsg;
      info() << "ECAL theta cells " << extremaECalLastLayerTheta.first << " - " << extremaECalLastLayerTheta.second
	     << " will be matched to HCAL " << m_activeFieldNamesNested[2] << "(s) " << extremaHCalFirstLayerZ.first
	     << " - " << extremaHCalFirstLayerZ.second << endmsg;
    }
    else{
      info() << "ECAL layer " << eCalLastLayer << " is a neighbour of HCAL layer 0." << endmsg;
      info() << "ECAL phi cells " << extremaECalLastLayerModule.first << " - " << extremaECalLastLayerModule.second
	     << " will be matched to HCAL cells " << extremaHCalFirstLayerPhi.first
	     << " - " << extremaHCalFirstLayerPhi.second << endmsg;
      info() << "ECAL theta cells " << extremaECalLastLayerTheta.first << " - " << extremaECalLastLayerTheta.second
	     << " will be matched to HCAL " << extremaHCalFirstLayerTheta.first
	     << " - " << extremaHCalFirstLayerTheta.second << endmsg;
    }
    
    std::unordered_map<uint, std::vector<uint>> thetaNeighbours;
    std::unordered_map<uint, std::vector<uint>> phiNeighbours;
    double hCalPhiSize = 2 * M_PI / (extremaHCalFirstLayerPhi.second - extremaHCalFirstLayerPhi.first + 1);
    // loop over z and find which theta cells to add
    if (m_readoutNamesNested.size()!=0){
      for (int iZ = 0; iZ < extremaHCalFirstLayerZ.second + 1; iZ++) {
	double lowZ = m_hCalZOffset + iZ * m_hCalZSize;
	double highZ = m_hCalZOffset + (iZ + 1) * m_hCalZSize;
	double lowTheta = 0, highTheta = 0;
	if (fabs(lowZ) < 1e-3) {
	  lowTheta = 0;
	} else {
	  lowTheta =
////TODO
            lowZ / fabs(lowZ) * atan(m_hCalRinner / lowZ); 
            //lowZ / fabs(lowZ) * (-log(fabs(tan(atan(m_hCalRinner / lowZ) / 2.))));  // theta = atan(m_hCalRinner / lowZ)
	}
	if (fabs(highZ) < 1e-3) {
	  highTheta = 0;
	} else {
	  highTheta = highZ / fabs(highZ) * (-log(fabs(tan(atan(m_hCalRinner / highZ) / 2.))));
	}
	debug() << "HCal z id  : " << iZ << endmsg;
	debug() << "HCal theta range  : " << lowTheta << " -  " << highTheta << endmsg;
	int lowId = floor((lowTheta - 0.5 * eCalThetaSize - eCalThetaOffset) / eCalThetaSize);
	int highId = floor((highTheta + 0.5 * eCalThetaSize - eCalThetaOffset) / eCalThetaSize);
	debug() << "ECal theta range  : " << lowId * eCalThetaSize + eCalThetaOffset << " -  "
		<< highId * eCalThetaSize + eCalThetaOffset << endmsg;
	std::vector<uint> neighbours;
	for (int idThetaToAdd = lowId; idThetaToAdd <= highId; idThetaToAdd++) {
	  if (idThetaToAdd >= extremaECalLastLayerTheta.first && idThetaToAdd <= extremaECalLastLayerTheta.second) {
	    neighbours.push_back(idThetaToAdd);
	  }
	}
	debug() << "HCal z id  : " << iZ << endmsg;
	debug() << "Found ECal Neighbours in theta : " << neighbours.size() << endmsg;
	for (auto id : neighbours) {
	  debug() << "ECal Neighbours id : " << id << endmsg;
	}
	thetaNeighbours.insert(std::pair<uint, std::vector<uint>>(iZ, neighbours));
      }
    }
    else{ // loop over theta cells to match in theta
      for (int iTheta = extremaHCalFirstLayerTheta.first; iTheta < extremaHCalFirstLayerTheta.second + 1; iTheta++) {
	double lowTheta = hCalThetaOffset + iTheta * hCalThetaSize;
	double highTheta = hCalThetaOffset + (iTheta + 1) * hCalThetaSize;
	debug() << "HCal theta range  : " << lowTheta << " -  " << highTheta << endmsg;
	int lowId = floor((lowTheta - 0.5 * eCalThetaSize - eCalThetaOffset) / eCalThetaSize);
	int highId = floor((highTheta + 0.5 * eCalThetaSize - eCalThetaOffset) / eCalThetaSize);
	debug() << "ECal theta range  : " << lowId * eCalThetaSize + eCalThetaOffset << " -  "
		<< highId * eCalThetaSize + eCalThetaOffset << endmsg;
	std::vector<uint> neighbours;
	for (int idThetaToAdd = lowId; idThetaToAdd <= highId; idThetaToAdd++) {
	  neighbours.push_back(det::utils::cyclicNeighbour(idThetaToAdd, extremaECalLastLayerTheta));
	}
	debug() << "HCal theta id  : " << iTheta << endmsg;
	debug() << "Found ECal Neighbours in theta : " << neighbours.size() << endmsg;
	for (auto id : neighbours) {
	  debug() << "ECal Neighbours id : " << id << endmsg;
	}
	thetaNeighbours.insert(std::pair<uint, std::vector<uint>>(iTheta, neighbours));
      }
    }
    // loop over phi and find which phi cells to add
    for (int iPhi = 0; iPhi < extremaHCalFirstLayerPhi.second +1; iPhi++) {
      double lowPhi = hCalPhiOffset + iPhi * hCalPhiSize;
      double highPhi = hCalPhiOffset + (iPhi + 1) * hCalPhiSize;
      debug() << "HCal phi range  : " << lowPhi << " -  " << highPhi << endmsg;
      int lowId = floor((lowPhi - 0.5 * eCalModuleSize - eCalPhiOffset) / eCalModuleSize);
      int highId = floor((highPhi + 0.5 * eCalModuleSize - eCalPhiOffset) / eCalModuleSize);
      debug() << "ECal phi range  : " << lowId * eCalModuleSize + eCalPhiOffset << " -  "
             << highId * eCalModuleSize + eCalPhiOffset << endmsg;
      std::vector<uint> neighbours;
      for (int idPhiToAdd = lowId; idPhiToAdd <= highId; idPhiToAdd++) {
        neighbours.push_back(det::utils::cyclicNeighbour(idPhiToAdd, extremaECalLastLayerModule));
      }
      debug() << "HCal phi id  : " << iPhi << endmsg;
      debug() << "Found ECal Neighbours in phi : " << neighbours.size() << endmsg;
      for (auto id : neighbours) {
        debug() << "ECal Neighbours id : " << id << endmsg;
      }
      phiNeighbours.insert(std::pair<uint, std::vector<uint>>(iPhi, neighbours));
    }
    // add neighbours to both ecal cell and hcal cells
    dd4hep::DDSegmentation::CellID ecalCellId = 0;
    dd4hep::DDSegmentation::CellID hcalCellId = 0;
    (*decoderECalBarrel)["system"].set(ecalCellId, 5);
    (*decoderECalBarrel)[m_activeFieldNamesSegmented[0]].set(ecalCellId, eCalLastLayer);
    (*decoderHCalBarrel)["system"].set(hcalCellId, 8);
    // loop over nested hcal cells
    if (m_readoutNamesNested.size()!=0){
      (*decoderHCalBarrel)[m_activeFieldNamesNested[0]].set(hcalCellId, 0);
      for (auto iZ : thetaNeighbours) {
	(*decoderHCalBarrel)[m_activeFieldNamesNested[2]].set(hcalCellId, iZ.first);
	for (auto iMod : phiNeighbours) {
	  (*decoderHCalBarrel)[m_activeFieldNamesNested[1]].set(hcalCellId, iMod.first);
	  for (auto iTheta : iZ.second) {
	    (*decoderECalBarrel)["theta"].set(ecalCellId, iTheta);
	    for (auto iPhi : iMod.second) {
	      (*decoderECalBarrel)["phi"].set(ecalCellId, iPhi);
	      map.find(hcalCellId)->second.push_back(ecalCellId);
	      map.find(ecalCellId)->second.push_back(hcalCellId);
	      count++;
	    }
	  }
	}
      }
    }
    // loop over segmented hcal cells
    else {
      (*decoderHCalBarrel)[m_activeFieldNamesSegmented[1]].set(hcalCellId, 0);
      for (auto iThetaHCal : thetaNeighbours) {
	(*decoderHCalBarrel)["theta"].set(hcalCellId, iThetaHCal.first); 
	for (auto iPhiHCal : phiNeighbours) {
	  (*decoderHCalBarrel)["phi"].set(hcalCellId, iPhiHCal.first);
	  for (auto iTheta : iThetaHCal.second) {
	    (*decoderECalBarrel)["theta"].set(ecalCellId, iTheta);
	    for (auto iPhi : iPhiHCal.second) {
	      (*decoderECalBarrel)["phi"].set(ecalCellId, iPhi);
	      map.find(hcalCellId)->second.push_back(ecalCellId);
	      map.find(ecalCellId)->second.push_back(hcalCellId);
	      count ++;
	    }
	  }
	}
      }
    }
  }
*/
  if (msgLevel() <= MSG::DEBUG) {
    std::vector<int> counter;
    counter.assign(40, 0);
    for (const auto& item : map) {
      counter[item.second.size()]++;
    }
    for (uint iCount = 0; iCount < counter.size(); iCount++) {
      if (counter[iCount] != 0) {
        debug() << counter[iCount] << " cells have " << iCount << " neighbours" << endmsg;
      }
    }
  }
  debug() << "cells with neighbours across Calo boundaries: " << count << endmsg;

  std::unique_ptr<TFile> file(TFile::Open(m_outputFileName.c_str(), "RECREATE"));
  file->cd();
  TTree tree("neighbours", "Tree with map of neighbours");
  uint64_t saveCellId;
  std::vector<uint64_t> saveNeighbours;
  tree.Branch("cellId", &saveCellId, "cellId/l");
  tree.Branch("neighbours", &saveNeighbours);
  for (const auto& item : map) {
    saveCellId = item.first;
    saveNeighbours = item.second;
    tree.Fill();
  }
  file->Write();
  file->Close();

  return StatusCode::SUCCESS;
}

StatusCode CreateFCCeeCaloNeighbours::finalize() { return Service::finalize(); }
