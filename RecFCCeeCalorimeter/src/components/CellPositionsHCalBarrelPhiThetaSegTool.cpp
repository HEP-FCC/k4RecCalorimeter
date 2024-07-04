#include "CellPositionsHCalBarrelPhiThetaSegTool.h"

#include "edm4hep/CalorimeterHitCollection.h"
#include <DDRec/DetectorData.h>

using dd4hep::DetElement;


DECLARE_COMPONENT(CellPositionsHCalBarrelPhiThetaSegTool)

CellPositionsHCalBarrelPhiThetaSegTool::CellPositionsHCalBarrelPhiThetaSegTool(
    const std::string &type,
    const std::string &name,
    const IInterface *parent)
    : GaudiTool(type, name, parent)
{
  declareInterface<ICellPositionsTool>(this);
}

StatusCode CellPositionsHCalBarrelPhiThetaSegTool::initialize()
{
  StatusCode sc = GaudiTool::initialize();
  if (sc.isFailure())
    return sc;
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc)
  {
    error() << "Unable to locate Geometry service." << endmsg;
    return StatusCode::FAILURE;
  }
  // get PhiTheta segmentation
  m_segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo *>(
      m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());
  if (m_segmentation == nullptr)
  {
    error() << "There is no phi-theta segmentation!!!!" << endmsg;
    // return StatusCode::FAILURE;
  }
  // Take readout bitfield decoder from GeoSvc
  m_decoder = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder();
  m_volman = m_geoSvc->getDetector()->volumeManager();

  // check if decoder contains "layer"
  std::vector<std::string> fields;
  for (uint itField = 0; itField < m_decoder->size(); itField++)
  {
    fields.push_back((*m_decoder)[itField].name());
  }
  auto iter = std::find(fields.begin(), fields.end(), "layer");
  if (iter == fields.end())
  {
    error() << "Readout does not contain field: 'layer'" << endmsg;
  }

  // needed to retrieve radii from the dd4hep extension
  dd4hep::Detector* detector = m_geoSvc->getDetector();
  if (!detector)
  {
    error() << "Unable to retrieve the detector." << endmsg;
    return StatusCode::FAILURE;
  }

  DetElement caloDetElem = detector->detector(m_detectorName);
  if (!caloDetElem.isValid())
  {
    error() << "Unable to retrieve the detector element: " << m_detectorName << endmsg;
    return StatusCode::FAILURE;
  }

  dd4hep::rec::LayeredCalorimeterData* theExtension = caloDetElem.extension<dd4hep::rec::LayeredCalorimeterData>();
  if (!theExtension)
  {
    error() << "The detector element does not have the required LayeredCalorimeterData extension." << endmsg;
    return StatusCode::FAILURE;
  }

  // Debug information to verify the layers retrieved
  // m_layersRetrieved is a pointer, then it needs to point to the address of the vector
  m_layersRetrieved = &(theExtension->layers);
  debug() << "Number of layers retrieved: " << m_layersRetrieved->size() << endmsg;

  if (m_detectorName=="HCalBarrel"){
    m_radii = CellPositionsHCalBarrelPhiThetaSegTool::calculateLayerRadiiBarrel();
  }
  else if (m_detectorName=="HCalThreePartsEndcap"){
    m_radii = CellPositionsHCalBarrelPhiThetaSegTool::calculateLayerRadiiEndcap();
  }
  else{
    error() << "Provided detector name in m_detectorName " << m_detectorName << " is not matching, expected inputs are HCalBarrel or HCalThreePartsEndcap" << endmsg;
    return StatusCode::FAILURE;
  }

  unsigned int numLayersProvided = 0;

  if (m_detectorName=="HCalThreePartsEndcap")
  {
    for (unsigned int i=0; i<m_numLayersHCalThreeParts.size(); i++) 
    {
      numLayersProvided += m_numLayersHCalThreeParts[i];
    }

    if (m_layersRetrieved->size() != numLayersProvided)
    {
      error() << "Total number of radial layers provided in m_numLayersHCalThreeParts " << numLayersProvided << " does not match the numbers retrieved from the xml file " << m_layersRetrieved->size() << endmsg;
      return StatusCode::FAILURE;
    }
  }
  return sc;
}

void CellPositionsHCalBarrelPhiThetaSegTool::getPositions(const edm4hep::CalorimeterHitCollection &aCells,
                                                          edm4hep::CalorimeterHitCollection &outputColl)
{
  debug() << "Input collection size : " << aCells.size() << endmsg;
  // Loop through cell collection
  for (const auto &cell : aCells)
  {
    auto outSeg = CellPositionsHCalBarrelPhiThetaSegTool::xyzPosition(cell.getCellID());

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
            << outSeg.z() / dd4hep::mm << endmsg;
  }
  debug() << "Output positions collection size: " << outputColl.size() << endmsg;
}


dd4hep::Position CellPositionsHCalBarrelPhiThetaSegTool::xyzPosition(const uint64_t &aCellId) const
{
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

  // MM: the z-coordinate is not correct for Endcap 
  debug() << "Layer : " << layer << "\tradius : " << radius << " cm" << endmsg;
  debug() << "Local position : x = " << inSeg.x() << " y = "  << inSeg.y() << " z = "  << inSeg.z() << endmsg;
  debug() << "Global position : x = " << outSeg.x() << " y = "  << outSeg.y() << " z = "  << outSeg.z() << endmsg;

  return outSeg;
}


int CellPositionsHCalBarrelPhiThetaSegTool::layerId(const uint64_t &aCellId)
{
  int layer;
  layer = m_decoder->get(aCellId, "layer");
  return layer;
}

std::vector<double> CellPositionsHCalBarrelPhiThetaSegTool::calculateLayerRadii(unsigned int startIndex, unsigned int endIndex) {
  double totalDepth = 0;
  std::vector<double> radii;

  for (unsigned int idxLayer = startIndex; idxLayer < endIndex; ++idxLayer) {
    const dd4hep::rec::LayeredCalorimeterStruct::Layer & theLayer = m_layersRetrieved->at(idxLayer);
    double distanceFirstLayer = theLayer.distance;
    double layerThickness = theLayer.sensitive_thickness;

    // Initialize totalDepth for the first layer in the range
    if (idxLayer == startIndex) {
      totalDepth = distanceFirstLayer + layerThickness;
    } else {
      totalDepth += layerThickness;
    }

    // Calculate radius for the current layer
    double radiusLayer = totalDepth - layerThickness / 2;
    radii.push_back(radiusLayer);
  }
  return radii;
}

// to be used for HCalTileBarrel 
std::vector<double> CellPositionsHCalBarrelPhiThetaSegTool::calculateLayerRadiiBarrel() {
  return calculateLayerRadii(0, m_layersRetrieved->size());
}

// to be used for HCalThreePartsEndcap
std::vector<double> CellPositionsHCalBarrelPhiThetaSegTool::calculateLayerRadiiEndcap() {
  std::vector<double> radii;

  // Calculate radii for each part and merge the results
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

StatusCode CellPositionsHCalBarrelPhiThetaSegTool::finalize() { return GaudiTool::finalize(); }
