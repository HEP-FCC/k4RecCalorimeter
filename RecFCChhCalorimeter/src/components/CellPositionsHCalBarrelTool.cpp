#include "CellPositionsHCalBarrelTool.h"

#include "edm4hep/CalorimeterHitCollection.h"

DECLARE_COMPONENT(CellPositionsHCalBarrelTool)

CellPositionsHCalBarrelTool::CellPositionsHCalBarrelTool(const std::string& type, const std::string& name,
                                                         const IInterface* parent)
    : GaudiTool(type, name, parent) {
  declareInterface<ICellPositionsTool>(this);
}

StatusCode CellPositionsHCalBarrelTool::initialize() {
  StatusCode sc = GaudiTool::initialize();
  if (sc.isFailure()) return sc;
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry service." << endmsg;
    return StatusCode::FAILURE;
  }
  // get PhiEta segmentation
  m_segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridPhiEta*>(
      m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());
  if (m_segmentation == nullptr) {
    error() << "There is no phi-eta segmentation!!!!" << endmsg;
    // return StatusCode::FAILURE;
  }
  // Take readout bitfield decoder from GeoSvc
  m_decoder = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder();
  m_volman = m_geoSvc->getDetector()->volumeManager();
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

void CellPositionsHCalBarrelTool::getPositions(const edm4hep::CalorimeterHitCollection& aCells,
                                               edm4hep::CalorimeterHitCollection& outputColl) {
  debug() << "Input collection size : " << aCells.size() << endmsg;
  // Loop through cell collection
  for (const auto& cell : aCells) {
    auto outSeg = CellPositionsHCalBarrelTool::xyzPosition(cell.getCellID());

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

dd4hep::Position CellPositionsHCalBarrelTool::xyzPosition(const uint64_t& aCellId) const {
  dd4hep::DDSegmentation::CellID volumeId = aCellId;
  m_decoder->set(volumeId,"phi",0);
  m_decoder->set(volumeId,"eta",0);

  int layer = m_decoder->get(volumeId, "layer");
    
  // radius calculated from segmenation + z postion of volumes
  auto inSeg = m_segmentation->position(aCellId);
  // get radius in cm
  double radius = m_radii[layer];
  dd4hep::Position outSeg(inSeg.x() * radius, inSeg.y() * radius, inSeg.z() * radius);
  return outSeg;
}

int CellPositionsHCalBarrelTool::layerId(const uint64_t& aCellId) {
  int layer;
  layer = m_decoder->get(aCellId,"layer");
  return layer;
}

StatusCode CellPositionsHCalBarrelTool::finalize() { return GaudiTool::finalize(); }
