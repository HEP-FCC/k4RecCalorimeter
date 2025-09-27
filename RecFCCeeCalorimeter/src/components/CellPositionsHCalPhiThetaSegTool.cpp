#include "CellPositionsHCalPhiThetaSegTool.h"

#include "edm4hep/CalorimeterHitCollection.h"
#include <DDRec/DetectorData.h>

using dd4hep::DetElement;

DECLARE_COMPONENT(CellPositionsHCalPhiThetaSegTool)

CellPositionsHCalPhiThetaSegTool::CellPositionsHCalPhiThetaSegTool(const std::string& type, const std::string& name,
                                                                   const IInterface* parent)
    : AlgTool(type, name, parent) {
  declareInterface<ICellPositionsTool>(this);
}

StatusCode CellPositionsHCalPhiThetaSegTool::initialize() {
  StatusCode sc = AlgTool::initialize();
  if (sc.isFailure())
    return sc;
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry service." << endmsg;
    return StatusCode::FAILURE;
  }

  // get the segmentation class type
  m_segmentationType = m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation()->type();

  if (m_segmentationType == "FCCSWGridPhiTheta_k4geo") {
    // get GridPhiTheta segmentation
    m_segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo*>(
        m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());
  }
  if (m_segmentationType == "FCCSWHCalPhiTheta_k4geo") {
    // get PhiTheta segmentation
    m_segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWHCalPhiTheta_k4geo*>(
        m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());
  }
  if (m_segmentationType == "FCCSWHCalPhiRow_k4geo") {
    // get PhiRow segmentation
    m_segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWHCalPhiRow_k4geo*>(
        m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());
  }

  if (m_segmentation == nullptr) {
    error() << "There is no phi-theta segmentation!!!!" << endmsg;
    return StatusCode::FAILURE;
  }
  // Take readout bitfield decoder from GeoSvc
  m_decoder = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder();

  // check if decoder contains "layer"
  std::vector<std::string> fields;
  for (uint itField = 0; itField < m_decoder->size(); itField++) {
    fields.push_back((*m_decoder)[itField].name());
  }
  auto iter = std::find(fields.begin(), fields.end(), "layer");
  if (iter == fields.end()) {
    error() << "Readout does not contain field: 'layer'" << endmsg;
    return StatusCode::FAILURE;
  }

  // retrieve radii from the LayeredCalorimeterData extension
  dd4hep::Detector* detector = m_geoSvc->getDetector();
  if (!detector) {
    error() << "Unable to retrieve the detector." << endmsg;
    return StatusCode::FAILURE;
  }

  DetElement caloDetElem = detector->detector(m_detectorName);
  if (!caloDetElem.isValid()) {
    error() << "Unable to retrieve the detector element: " << m_detectorName << endmsg;
    return StatusCode::FAILURE;
  }

  dd4hep::rec::LayeredCalorimeterData* theExtension = caloDetElem.extension<dd4hep::rec::LayeredCalorimeterData>();
  if (!theExtension) {
    error() << "The detector element does not have the required LayeredCalorimeterData extension." << endmsg;
    return StatusCode::FAILURE;
  }

  // Debug information to check the number of layers retrieved from the LayeredCalorimeterData extension
  m_layersRetrieved = &(theExtension->layers);
  debug() << "Number of layers retrieved: " << m_layersRetrieved->size() << endmsg;

  if (m_detectorName == "HCalBarrel") {
    m_radii = CellPositionsHCalPhiThetaSegTool::calculateLayerRadiiBarrel();
  } else if (m_detectorName == "HCalThreePartsEndcap") {
    m_radii = CellPositionsHCalPhiThetaSegTool::calculateLayerRadiiEndcap();
  } else {
    error() << "Provided detector name in m_detectorName " << m_detectorName
            << " is not matching, expected inputs are HCalBarrel or HCalThreePartsEndcap" << endmsg;
    return StatusCode::FAILURE;
  }

  unsigned int numLayersProvided = 0;

  if (m_detectorName == "HCalThreePartsEndcap") {
    // Check that the vector containing number of layers in each cylinder is provided
    if (m_numLayersHCalThreeParts.empty()) {
      error() << "The vector m_numLayersHCalThreeParts is empty." << endmsg;
      return StatusCode::FAILURE;
    }
    // Check that the total number of layers provided in the steering file
    // matches the total number of layers in the geometry xml file
    for (unsigned int i = 0; i < m_numLayersHCalThreeParts.size(); i++) {
      numLayersProvided += m_numLayersHCalThreeParts[i];
    }

    if (m_layersRetrieved->size() != numLayersProvided) {
      error() << "Total number of radial layers provided in m_numLayersHCalThreeParts " << numLayersProvided
              << " does not match the numbers retrieved from the xml file " << m_layersRetrieved->size() << endmsg;
      return StatusCode::FAILURE;
    }
  }
  return sc;
}

void CellPositionsHCalPhiThetaSegTool::getPositions(const edm4hep::CalorimeterHitCollection& aCells,
                                                    edm4hep::CalorimeterHitCollection& outputColl) const {
  debug() << "Input collection size : " << aCells.size() << endmsg;
  // Loop through cell collection
  for (const auto& cell : aCells) {
    auto outSeg = CellPositionsHCalPhiThetaSegTool::xyzPosition(cell.getCellID());

    auto edmPos = edm4hep::Vector3f();
    edmPos.x = outSeg.x() / dd4hep::mm;
    edmPos.y = outSeg.y() / dd4hep::mm;
    edmPos.z = outSeg.z() / dd4hep::mm;

    auto positionedHit = cell.clone();
    positionedHit.setPosition(edmPos);
    outputColl.push_back(positionedHit);

    // Debug information about cell position
    debug() << "Cell energy (GeV) : " << positionedHit.getEnergy() << "\tcellID " << positionedHit.getCellID()
            << endmsg;
    debug() << "Position of cell (mm) : \t" << outSeg.x() / dd4hep::mm << "\t" << outSeg.y() / dd4hep::mm << "\t"
            << outSeg.z() / dd4hep::mm << endmsg;
  }
  debug() << "Output positions collection size: " << outputColl.size() << endmsg;
}

dd4hep::Position CellPositionsHCalPhiThetaSegTool::xyzPosition(const uint64_t& aCellId) const {
  // retrieve position for FCCSWHCalPhiTheta_k4geo and FCCSWHCalPhiRow_k4geo segmentation types
  if (m_segmentationType == "FCCSWHCalPhiTheta_k4geo" || m_segmentationType == "FCCSWHCalPhiRow_k4geo") {
    // get global position
    auto posCell = m_segmentation->position(aCellId);
    dd4hep::Position outSeg(posCell.x(), posCell.y(), posCell.z());

    debug() << "Layer : " << m_decoder->get(aCellId, "layer") << endmsg;
    debug() << "Global position : x = " << outSeg.x() << " y = " << outSeg.y() << " z = " << outSeg.z() << endmsg;

    return outSeg;
  }

  // following code is valid for FCCSWGridPhiTheta_k4geo segmentation type.

  dd4hep::DDSegmentation::CellID volumeId = aCellId;
  m_decoder->set(volumeId, "phi", 0);
  m_decoder->set(volumeId, "theta", 0);

  int layer = m_decoder->get(volumeId, "layer");
  // get radius in cm
  double radius = m_radii[layer];
  // get local position (for r=1)
  auto inSeg = m_segmentation->position(aCellId);
  // scale by radius to get global position
  dd4hep::Position outSeg(inSeg.x() * radius, inSeg.y() * radius, inSeg.z() * radius);

  // MM: TBD the z-coordinate still needs to be carefully validated
  // at the first glance it seems to be off in some cases in the Endcap
  debug() << "Layer : " << layer << "\tradius : " << radius << " cm" << endmsg;
  debug() << "Local position : x = " << inSeg.x() << " y = " << inSeg.y() << " z = " << inSeg.z() << endmsg;
  debug() << "Global position : x = " << outSeg.x() << " y = " << outSeg.y() << " z = " << outSeg.z() << endmsg;

  return outSeg;
}

int CellPositionsHCalPhiThetaSegTool::layerId(const uint64_t& aCellId) const {
  int layer;
  layer = m_decoder->get(aCellId, "layer");
  return layer;
}

// calculate layer radii from LayeredCalorimeterData extension
// which is included in the geometry description
std::vector<double> CellPositionsHCalPhiThetaSegTool::calculateLayerRadii(unsigned int startIndex,
                                                                          unsigned int endIndex) {
  std::vector<double> radii;

  for (unsigned int idxLayer = startIndex; idxLayer < endIndex; ++idxLayer) {
    const dd4hep::rec::LayeredCalorimeterStruct::Layer& theLayer = m_layersRetrieved->at(idxLayer);
    // inner radius of a given layer
    double layerInnerRadius = theLayer.distance;
    // radial dimension of a given layer
    double layerThickness = theLayer.sensitive_thickness;

    // Calculate mid-radius for the current layer (layerInnerRadius+layerOuterRadius)/2
    double layerMidRadius = layerInnerRadius + layerThickness / 2;

    radii.push_back(layerMidRadius);
  }
  return radii;
}

// calculateLayerRadii should be used for HCalTileBarrel which is formed
// by a single cylinder with a constant number of radial layers
std::vector<double> CellPositionsHCalPhiThetaSegTool::calculateLayerRadiiBarrel() {
  return calculateLayerRadii(0, m_layersRetrieved->size());
}

// calculateLayerRadiiEndcap should be used for HCalThreePartsEndcap which is formed
// by three cylinders along the z-coordinate and each cylinder has a different number of radial layers
std::vector<double> CellPositionsHCalPhiThetaSegTool::calculateLayerRadiiEndcap() {
  std::vector<double> radii;

  // Calculate radii for each part (cylinder) and merge the results
  unsigned int upperIndexLayerRadiiPartOne = m_numLayersHCalThreeParts[0];
  unsigned int upperIndexLayerRadiiPartTwo = upperIndexLayerRadiiPartOne + m_numLayersHCalThreeParts[1];

  // Part 1
  auto partOneRadii = calculateLayerRadii(0, upperIndexLayerRadiiPartOne);
  radii.insert(radii.end(), partOneRadii.begin(), partOneRadii.end());

  // Part 2
  auto partTwoRadii = calculateLayerRadii(upperIndexLayerRadiiPartOne, upperIndexLayerRadiiPartTwo);
  radii.insert(radii.end(), partTwoRadii.begin(), partTwoRadii.end());

  // Part 3
  auto partThreeRadii = calculateLayerRadii(upperIndexLayerRadiiPartTwo, m_layersRetrieved->size());
  radii.insert(radii.end(), partThreeRadii.begin(), partThreeRadii.end());

  return radii;
}

StatusCode CellPositionsHCalPhiThetaSegTool::finalize() { return AlgTool::finalize(); }
