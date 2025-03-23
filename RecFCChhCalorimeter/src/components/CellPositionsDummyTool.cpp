#include "CellPositionsDummyTool.h"

#include "edm4hep/CalorimeterHitCollection.h"

DECLARE_COMPONENT(CellPositionsDummyTool)

CellPositionsDummyTool::CellPositionsDummyTool(const std::string& type, const std::string& name,
                                               const IInterface* parent)
    : AlgTool(type, name, parent) {
  declareInterface<ICellPositionsTool>(this);
}

StatusCode CellPositionsDummyTool::initialize() {
  StatusCode sc = AlgTool::initialize();
  if (sc.isFailure())
    return sc;
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry service." << endmsg;
    return StatusCode::FAILURE;
  }
  return sc;
}

void CellPositionsDummyTool::getPositions(const edm4hep::CalorimeterHitCollection& aCells,
                                          edm4hep::CalorimeterHitCollection& outputColl) {
  debug() << "Input collection size : " << aCells.size() << endmsg;
  // Loop through cell collection
  for (const auto& cell : aCells) {
    auto outPos = CellPositionsDummyTool::xyzPosition(cell.getCellID());

    auto edmPos = edm4hep::Vector3f();
    edmPos.x = outPos.x() / dd4hep::mm;
    edmPos.y = outPos.y() / dd4hep::mm;
    edmPos.z = outPos.z() / dd4hep::mm;
    auto positionedHit = cell.clone();
    positionedHit.setPosition(edmPos);
    outputColl.push_back(positionedHit);
    // Debug information about cell position
    debug() << "Cell energy (GeV) : " << positionedHit.getEnergy() << "\tcellID " << positionedHit.getCellID()
            << endmsg;
    debug() << "Position of cell (mm) : \t" << outPos.x() / dd4hep::mm << "\t" << outPos.y() / dd4hep::mm << "\t"
            << outPos.z() / dd4hep::mm << endmsg;
  }
  debug() << "Output positions collection size: " << outputColl.size() << endmsg;
}

dd4hep::Position CellPositionsDummyTool::xyzPosition(const uint64_t& /*aCellId*/) const {
  dd4hep::Position outPos;
  outPos.SetCoordinates(0, 0, 0);
  return outPos;
}

int CellPositionsDummyTool::layerId(const uint64_t& /*aCellId*/) {
  int layer = 0;
  return layer;
}

StatusCode CellPositionsDummyTool::finalize() { return AlgTool::finalize(); }
