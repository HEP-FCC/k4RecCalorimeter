#include "CreateCellPositions.h"

// FCCSW
#include "detectorCommon/DetUtils_k4geo.h"
#include "k4Interface/IGeoSvc.h"

// DD4hep
#include "DD4hep/Detector.h"

// EDM
#include "edm4hep/CalorimeterHitCollection.h"

DECLARE_COMPONENT(CreateCellPositions)

CreateCellPositions::CreateCellPositions(const std::string& name, ISvcLocator* svcLoc)
    : Gaudi::Algorithm(name, svcLoc) {
  declareProperty("hits", m_hits, "Hit collection (input)");
  declareProperty("positionsTool", m_cellPositionsTool, "Handle for tool to retrieve cell positions");
  declareProperty("positionedHits", m_positionedHits, "Output cell positions collection");
}

StatusCode CreateCellPositions::initialize() {
  StatusCode sc = Gaudi::Algorithm::initialize();
  if (sc.isFailure())
    return sc;
  if (!m_cellPositionsTool.retrieve()) {
    error() << "Unable to retrieve the cell positions tool!!!" << endmsg;
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}

StatusCode CreateCellPositions::execute(const EventContext&) const {
  // Get the input hit collection
  const auto* hits = m_hits.get();
  debug() << "Input hit collection size: " << hits->size() << endmsg;
  // Initialize output collection
  auto edmPositionedHitCollection = m_positionedHits.createAndPut();

  m_cellPositionsTool->getPositions(*hits, *edmPositionedHitCollection);

  debug() << "Output positions collection size: " << edmPositionedHitCollection->size() << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode CreateCellPositions::finalize() { return Gaudi::Algorithm::finalize(); }
