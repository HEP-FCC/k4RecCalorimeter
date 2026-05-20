#include "CellPositionsDummyTool.h"
#include "k4FWCore/GaudiChecks.h"

#include "edm4hep/CalorimeterHitCollection.h"

DECLARE_COMPONENT(CellPositionsDummyTool)

StatusCode CellPositionsDummyTool::initialize() {
  K4_GAUDI_CHECK(AlgTool::initialize());
  K4_GAUDI_CHECK(m_geoSvc.retrieve());
  return StatusCode::SUCCESS;
}

void CellPositionsDummyTool::getPositions(const edm4hep::CalorimeterHitCollection& aCells,
                                          edm4hep::CalorimeterHitCollection& outputColl) const {
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

dd4hep::Position CellPositionsDummyTool::xyzPosition(const CellID /*aCellId*/) const {
  return dd4hep::Position(0, 0, 0);
}

int CellPositionsDummyTool::layerId(const CellID /*aCellId*/) const { return 0; }
