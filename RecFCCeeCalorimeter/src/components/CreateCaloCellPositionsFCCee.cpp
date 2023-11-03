#include "CreateCaloCellPositionsFCCee.h"

// k4
#include "detectorCommon/DetUtils_k4geo.h"
#include "k4Interface/IGeoSvc.h"

// DD4hep
#include "DD4hep/Detector.h"

// EDM
#include "edm4hep/CalorimeterHitCollection.h"

DECLARE_COMPONENT(CreateCaloCellPositionsFCCee)

CreateCaloCellPositionsFCCee::CreateCaloCellPositionsFCCee(const std::string& name, ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc) {
  declareProperty("hits", m_hits, "Hit collection (input)");
  declareProperty("positionsTool", m_cellPositionsTool,
                  "Handle for tool to retrieve cell positions from");
  declareProperty("positionedHits", m_positionedHits,
                  "Output cell positions collection");
}

StatusCode CreateCaloCellPositionsFCCee::initialize() {
  {
    StatusCode sc = GaudiAlgorithm::initialize();
    if (sc.isFailure()) return sc;
  }

  if (!m_cellPositionsTool) {
    error() << "CellPositionsTool is missing!" << endmsg;
    return StatusCode::FAILURE;
  }

  // Copy over the CellIDEncoding string from the input collection to the output collection
  m_positionedHitsCellIDEncoding.put(m_hitsCellIDEncoding.get());

  return StatusCode::SUCCESS;
}

StatusCode CreateCaloCellPositionsFCCee::execute() {
  // Get the input hit collection
  const auto* hits = m_hits.get();
  debug() << "Input hit collection size: " << hits->size() << endmsg;
  // Initialize output collection
  auto edmPositionedHitCollection = m_positionedHits.createAndPut();
  //edmPositionedHitCollection->reserve(hits->size()); // edm4hep uses std::deque ???

  for (const auto& hit : *hits) {
    auto positionedHit = hit.clone();
    dd4hep::DDSegmentation::CellID cellId = positionedHit.getCellID();

    auto cached_pos = m_positions_cache.find(cellId);
    if(cached_pos == m_positions_cache.end()) {
      // identify calo system
      dd4hep::Position posCell = m_cellPositionsTool->xyzPosition(cellId);

      edm4hep::Vector3f edmPos;
      edmPos.x = posCell.x() / dd4hep::mm;
      edmPos.y = posCell.y() / dd4hep::mm;
      edmPos.z = posCell.z() / dd4hep::mm;

      m_positions_cache[cellId] = edmPos;
      positionedHit.setPosition(edmPos);
    }
    else {
      positionedHit.setPosition(cached_pos->second);
    }

    // Debug information about cell position
    //debug() << "Cell energy (GeV) : " << positionedHit.getEnergy() << "\tcellID " << positionedHit.getCellID() << endmsg;
    //debug() << "Position of cell (mm) : \t" << edmPos.x << "\t" << edmPos.y << "\t" << edmPos.z << endmsg;

    edmPositionedHitCollection->push_back(positionedHit);
  }

  debug() << "Output positions collection size: " << edmPositionedHitCollection->size() << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode CreateCaloCellPositionsFCCee::finalize() { return GaudiAlgorithm::finalize(); }
