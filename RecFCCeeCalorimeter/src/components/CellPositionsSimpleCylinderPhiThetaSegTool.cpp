#include "CellPositionsSimpleCylinderPhiThetaSegTool.h"

// EDM4hep
#include "edm4hep/CalorimeterHitCollection.h"

// DD4hep
#include "DDRec/DetectorData.h"

DECLARE_COMPONENT(CellPositionsSimpleCylinderPhiThetaSegTool)

CellPositionsSimpleCylinderPhiThetaSegTool::CellPositionsSimpleCylinderPhiThetaSegTool(
  const std::string& type,
  const std::string& name,
  const IInterface* parent)
  : AlgTool(type, name, parent) {
  declareInterface<ICellPositionsTool>(this);
}

StatusCode CellPositionsSimpleCylinderPhiThetaSegTool::initialize() {

  // base class initialization
  StatusCode sc = AlgTool::initialize();
  if (sc.isFailure()) return sc;
  
  // get geometry service
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry service." << endmsg;
    return StatusCode::FAILURE;
  }

  // get the detector
  dd4hep::Detector* detector = m_geoSvc->getDetector();
  if (!detector) {
    error() << "Unable to retrieve the detector." << endmsg;
    return StatusCode::FAILURE;
  }
  
  // get the detector element
  dd4hep::DetElement detectorEl = detector->detector(m_detectorName);
  if (!detectorEl.isValid()) {
    error() << "Unable to retrieve the detector element: " << m_detectorName << endmsg;
    return StatusCode::FAILURE;
  }

  // get info about rmin, rmax
  dd4hep::rec::LayeredCalorimeterData* theExtension =  detectorEl.extension<dd4hep::rec::LayeredCalorimeterData>();
  if (!theExtension) {
    error() << "The detector element does not have the required LayeredCalorimeterData extension." << endmsg;
    return StatusCode::FAILURE;
  }
  m_detRadius = (theExtension->extent[0] + theExtension->extent[1]) / 2.0;
  if (theExtension->extent[2] > 0.)
    m_detZ = (theExtension->extent[2] + theExtension->extent[3]) / 2.0;
  else
    m_detZ = 0.0;
  debug() << "Radius of detector " << m_detectorName.value() << " [mm] : " << m_detRadius/dd4hep::mm << endmsg;
  debug() << "Z of detector " << m_detectorName.value() << " [mm] : " << m_detZ/dd4hep::mm << endmsg;

   // get layer positions
  const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>* layers = &(theExtension->layers) ;
  for (unsigned int idxLayer = 0; idxLayer < layers->size(); ++idxLayer) {
    const dd4hep::rec::LayeredCalorimeterStruct::Layer & theLayer = layers->at(idxLayer);
    // distance from inner face of layer to origin
    double layerInnerPosition = theLayer.distance;
    // layer thickness
    double layerThickness = theLayer.sensitive_thickness;
    // calculate the position of the center of current layer
    double position = layerInnerPosition + layerThickness / 2.;
    m_layerPositions.push_back(position);
  }

  // get phi-theta segmentation
  m_segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo*>(
      detector->readout(m_readoutName).segmentation().segmentation());
  if (m_segmentation == nullptr) {
    error() << "There is no phi-theta segmentation!!!!" << endmsg;
    return StatusCode::FAILURE;
  }

  // Take readout bitfield decoder
  m_decoder = detector->readout(m_readoutName).idSpec().decoder();
  // check if decoder contains "layer"
  std::vector<std::string> fields;
  for (uint itField = 0; itField < m_decoder->size(); itField++) {
    fields.push_back((*m_decoder)[itField].name());
  }
  auto iter = std::find(fields.begin(), fields.end(), "layer");
  if (iter == fields.end()) {
    error() << "Readout does not contain field: 'layer'" << endmsg;
  }
  
  return sc;
}

void CellPositionsSimpleCylinderPhiThetaSegTool::getPositions(const edm4hep::CalorimeterHitCollection& aCells,
                                                              edm4hep::CalorimeterHitCollection& outputColl) {

  debug() << "Input collection size : " << aCells.size() << endmsg;
  // Loop through cell collection
  for (const auto& cell : aCells) {
    auto outSeg = CellPositionsSimpleCylinderPhiThetaSegTool::xyzPosition(cell.getCellID());
    auto edmPos = edm4hep::Vector3f();
    edmPos.x = outSeg.x() / dd4hep::mm;
    edmPos.y = outSeg.y() / dd4hep::mm;
    edmPos.z = outSeg.z() / dd4hep::mm;

    auto positionedHit = cell.clone();
    positionedHit.setPosition(edmPos);
    outputColl.push_back(positionedHit);

    // Debug information about cell position
    debug() << "Cell energy (GeV) : " << positionedHit.getEnergy() << "\tcellID " << positionedHit.getCellID() << endmsg;
    debug() << "Position of cell (mm) : \t" << outSeg.x() / dd4hep::mm << "\t" << outSeg.y() / dd4hep::mm << "\t"
            << outSeg.z() / dd4hep::mm << "\n"
            << endmsg;
  }
  debug() << "Output positions collection size: " << outputColl.size() << endmsg;
}

dd4hep::Position CellPositionsSimpleCylinderPhiThetaSegTool::xyzPosition(const uint64_t& aCellId) const {
  debug() << "Cell ID: " << aCellId << endmsg;

  int layer = m_decoder->get(aCellId, "layer");

  // determine radius of cell
  double radius;
  if (m_detZ == 0.0) {
    // barrel
    // radius = m_detRadius;
    radius = m_layerPositions[layer];
  }
  else {
    // endcap
    auto theta = m_segmentation->theta(aCellId);
    // radius = fabs(m_detZ * tan(theta));
    radius = fabs(m_layerPositions[layer] * tan(theta));
  }
  // get position scaled to R=1
  auto inSeg = m_segmentation->position(aCellId);
  // rescale by radius
  dd4hep::Position outSeg(inSeg.x() * radius, inSeg.y() * radius, inSeg.z() * radius);

  return outSeg;
}

int CellPositionsSimpleCylinderPhiThetaSegTool::layerId(const uint64_t& aCellId) {
  int layer;
  dd4hep::DDSegmentation::CellID cID = aCellId;
  layer = m_decoder->get(cID, "layer");
  return layer;
}

StatusCode CellPositionsSimpleCylinderPhiThetaSegTool::finalize() { return AlgTool::finalize(); }
