#include "CreateCaloCellPositionsFCCee.h"

// FCCSW
#include "DetCommon/DetUtils.h"
#include "k4Interface/IGeoSvc.h"

// DD4hep
#include "DD4hep/Detector.h"

// EDM
#include "edm4hep/CalorimeterHitCollection.h"

DECLARE_COMPONENT(CreateCaloCellPositionsFCCee)

CreateCaloCellPositionsFCCee::CreateCaloCellPositionsFCCee(const std::string& name, ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc), m_eventDataSvc("EventDataSvc", "CreateCaloCellPositionsFCCee") {
  declareProperty("hits", m_hits, "Hit collection (input)");
  declareProperty("positionsECalBarrelTool", m_cellPositionsECalBarrelTool,
                  "Handle for tool to retrieve cell positions in ECal Barrel");
  declareProperty("positionsHCalBarrelTool", m_cellPositionsHCalBarrelTool,
                  "Handle for tool to retrieve cell positions in HCal Barrel and ext Barrel");
  declareProperty("positionsHCalExtBarrelTool", m_cellPositionsHCalExtBarrelTool,
                  "Handle for tool to retrieve cell positions in HCal Barrel and ext Barrel");
  declareProperty("positionsEMECTool", m_cellPositionsEMECTool, "Handle for tool to retrieve cell positions in EMEC");
  declareProperty("positionsHECTool", m_cellPositionsHECTool, "Handle for tool to retrieve cell positions in HEC");
  declareProperty("positionsEMFwdTool", m_cellPositionsEMFwdTool, "Handle for tool to retrieve cell positions EM Fwd");
  declareProperty("positionsHFwdTool", m_cellPositionsHFwdTool, "Handle for tool to retrieve cell positions Had Fwd");
  declareProperty("positionedHits", m_positionedHits, "Output cell positions collection");
}

StatusCode CreateCaloCellPositionsFCCee::initialize() {
  StatusCode sc = GaudiAlgorithm::initialize();
  if (sc.isFailure()) return sc;
  StatusCode sc_dataSvc = m_eventDataSvc.retrieve();
  m_podioDataSvc = dynamic_cast<PodioDataSvc*>(m_eventDataSvc.get());
  if (sc_dataSvc == StatusCode::FAILURE) {
    error() << "Error retrieving Event Data Service" << endmsg;
    return sc_dataSvc;
  }
  return StatusCode::SUCCESS;
}

StatusCode CreateCaloCellPositionsFCCee::execute() {
  // Get the input hit collection
  const auto* hits = m_hits.get();
  debug() << "Input hit collection size: " << hits->size() << endmsg;
  // Initialize output collection
  auto edmPositionedHitCollection = m_positionedHits.createAndPut();

  for (const auto& hit : *hits) {
    auto positionedHit = hit.clone();
    dd4hep::DDSegmentation::CellID cellId = positionedHit.getCellID();
    // identify calo system
    auto systemId = m_decoder->get(cellId, "system");
    dd4hep::Position posCell;

    if (systemId == 4)  // ECAL BARREL system id
      posCell = m_cellPositionsECalBarrelTool->xyzPosition(cellId);
    else if (systemId == 10)  // HCAL BARREL system id
      posCell = m_cellPositionsHCalBarrelTool->xyzPosition(cellId);
    else if (systemId == 9)  // HCAL EXT BARREL system id
      posCell = m_cellPositionsHCalExtBarrelTool->xyzPosition(cellId);
    else if (systemId == 6)  // EMEC system id
      posCell = m_cellPositionsEMECTool->xyzPosition(cellId);
    else if (systemId == 7)  // HEC system id
      posCell = m_cellPositionsHECTool->xyzPosition(cellId);
    else if (systemId == 10)  // EMFWD system id
      posCell = m_cellPositionsEMFwdTool->xyzPosition(cellId);
    else if (systemId == 11)  // HFWD system id
      posCell = m_cellPositionsHFwdTool->xyzPosition(cellId);

    auto edmPos = edm4hep::Vector3f();
    edmPos.x = posCell.x() / dd4hep::mm;
    edmPos.y = posCell.y() / dd4hep::mm;
    edmPos.z = posCell.z() / dd4hep::mm;

    positionedHit.setPosition(edmPos);
    edmPositionedHitCollection->push_back(positionedHit);

    // Debug information about cell position
    debug() << "Cell energy (GeV) : " << positionedHit.getEnergy() << "\tcellID " << positionedHit.getCellID() << endmsg;
    debug() << "Position of cell (mm) : \t" << posCell.x() / dd4hep::mm << "\t" << posCell.y() / dd4hep::mm << "\t"
            << posCell.z() / dd4hep::mm << endmsg;
  }
  auto& coll_md = m_podioDataSvc->getProvider().getCollectionMetaData(m_positionedHits.get()->getID());
  coll_md.setValue("CellIDEncodingString", m_hits.getCollMetadataCellID(hits->getID()));

  debug() << "Output positions collection size: " << edmPositionedHitCollection->size() << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode CreateCaloCellPositionsFCCee::finalize() { return GaudiAlgorithm::finalize(); }
