#include "LayerPhiEtaCaloTool.h"

// segm
#include "DD4hep/Detector.h"
#include "detectorCommon/DetUtils_k4geo.h"
#include "k4Interface/IGeoSvc.h"

DECLARE_COMPONENT(LayerPhiEtaCaloTool)

LayerPhiEtaCaloTool::LayerPhiEtaCaloTool(const std::string& type, const std::string& name,
                                                 const IInterface* parent)
    : GaudiTool(type, name, parent) {
  declareInterface<ICalorimeterTool>(this);
}

StatusCode LayerPhiEtaCaloTool::initialize() {
  StatusCode sc = GaudiTool::initialize();
  if (sc.isFailure()) return sc;
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service. "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_readoutName != "") {
    // Check if readouts exist
    info() << "Readout: " << m_readoutName << endmsg;
    if (m_geoSvc->getDetector()->readouts().find(m_readoutName) == m_geoSvc->getDetector()->readouts().end()) {
      error() << "Readout <<" << m_readoutName << ">> does not exist." << endmsg;
      return StatusCode::FAILURE;
    }
  }
  return sc;
}

StatusCode LayerPhiEtaCaloTool::finalize() { return GaudiTool::finalize(); }

StatusCode LayerPhiEtaCaloTool::prepareEmptyCells(std::unordered_map<uint64_t, double>& aCells) {
  // Get the total number of active volumes in the geometry
  auto highestVol = gGeoManager->GetTopVolume();
  unsigned int numLayers;
  if (!m_activeVolumesNumber) {
    numLayers = det::utils::countPlacedVolumes(highestVol, m_activeVolumeName);
  } else {
    // used when MergeLayers tool is used. To be removed once MergeLayer gets replaced by RedoSegmentation.
    numLayers = m_activeVolumesNumber;
  }
  info() << "Number of active layers " << numLayers << endmsg;
  
  // get PhiEta segmentation
  dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo* segmentation;
  segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo*>(
      m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());
  if (segmentation == nullptr) {
    error() << "There is no phi-eta segmentation!!!!" << endmsg;
    return StatusCode::FAILURE;
  }
  info() << "FCCSWGridPhiEta_k4geo: size in eta " << segmentation->gridSizeEta() << " , bins in phi " << segmentation->phiBins()
         << endmsg;
  info() << "FCCSWGridPhiEta_k4geo: offset in eta " << segmentation->offsetEta() << " , offset in phi "
         << segmentation->offsetPhi() << endmsg;

  // Take readout bitfield decoder from GeoSvc
  auto decoder = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder();
  if (m_fieldNames.size() != m_fieldValues.size()) {
    error() << "Volume readout field descriptors (names and values) have different size." << endmsg;
    return StatusCode::FAILURE;
  }

  // Get VolumeID
  dd4hep::DDSegmentation::VolumeID volumeID = 0;
  for (unsigned int it = 0; it < m_fieldNames.size(); it++) {
    decoder->set(volumeID, m_fieldNames[it], m_fieldValues[it]);
  }  

  if (m_activeVolumesEta.size() != numLayers){
    error() << "The given number of min eta is not equal to the number of layer!!!!" << endmsg;
    return StatusCode::FAILURE;
  }

  // Loop over all cells in the calorimeter and retrieve existing cellIDs
  // Loop over active layers
  for (unsigned int ilayer = 0; ilayer < numLayers; ilayer++) {
    // Get VolumeID
    decoder->set(volumeID, m_activeFieldName, ilayer);
    decoder->set(volumeID, "eta", 0);
    decoder->set(volumeID, "phi", 0);

    // Calculate number of cells per layer
    auto numCells = det::utils::numberOfCells(volumeID, *segmentation);
    uint cellsEta = ceil(( 2*m_activeVolumesEta[ilayer] - segmentation->gridSizeEta() ) / 2 / segmentation->gridSizeEta()) * 2 + 1; // ceil( 2*m_activeVolumesRadii[ilayer] / segmentation->gridSizeEta()) ;
    uint minEtaID = int(floor(( - m_activeVolumesEta[ilayer] + 0.5 * segmentation->gridSizeEta() - segmentation->offsetEta()) / segmentation->gridSizeEta()));
  
    numCells[1] = cellsEta; 
    numCells[2] = minEtaID; 
    debug() << "Segmentation cells  (Nphi, Neta, minEta): " << numCells << endmsg;
    // Loop over segmenation cells
    for (unsigned int iphi = 0; iphi < numCells[0]; iphi++) {
      for (unsigned int ieta = 0; ieta < numCells[1]; ieta++) {
        decoder->set(volumeID, "phi", iphi);
        decoder->set(volumeID, "eta", ieta + numCells[2]); // start from the minimum existing eta cell in this layer
        dd4hep::DDSegmentation::CellID cellId = volumeID;
	//debug() << "CellID: " <<  cellId<< endmsg;
        aCells.insert(std::pair<dd4hep::DDSegmentation::CellID, double>(cellId, 0));
      }
    }
  }
  return StatusCode::SUCCESS;
}
