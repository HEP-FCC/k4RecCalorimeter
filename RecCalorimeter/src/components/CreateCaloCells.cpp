#include "CreateCaloCells.h"

// FCCSW
#include "DetCommon/DetUtils.h"
#include "k4Interface/IGeoSvc.h"

// DD4hep
#include "DD4hep/Detector.h"

// EDM4HEP
#include "edm4hep/CalorimeterHit.h"

DECLARE_COMPONENT(CreateCaloCells)

CreateCaloCells::CreateCaloCells(const std::string& name, ISvcLocator* svcLoc) :
GaudiAlgorithm(name, svcLoc), m_geoSvc("GeoSvc", name) {
  declareProperty("hits", m_hits, "Hits from which to create cells (input)");
  declareProperty("cells", m_cells, "The created calorimeter cells (output)");

  declareProperty("calibTool", m_calibTool, "Handle for tool to calibrate Geant4 energy to EM scale tool");
  declareProperty("noiseTool", m_noiseTool, "Handle for the calorimeter cells noise tool");
  declareProperty("geometryTool", m_geoTool, "Handle for the geometry tool");
}

StatusCode CreateCaloCells::initialize() {
  StatusCode sc = GaudiAlgorithm::initialize();
  if (sc.isFailure()) return sc;

  info() << "CreateCaloCells initialized" << endmsg;
  info() << "do calibration : " << m_doCellCalibration << endmsg;
  info() << "add cell noise      : " << m_addCellNoise << endmsg;
  info() << "remove noise cells below threshold : " << m_filterCellNoise << endmsg;

  // Initialization of tools
  // Calibrate Geant4 energy to EM scale tool
  if (m_doCellCalibration) {
    if (!m_calibTool.retrieve()) {
      error() << "Unable to retrieve the calo cells calibration tool!!!" << endmsg;
      return StatusCode::FAILURE;
    }
  }
  // Cell noise tool
  if (m_addCellNoise) {
    if (!m_noiseTool.retrieve()) {
      error() << "Unable to retrieve the calo cells noise tool!!!" << endmsg;
      return StatusCode::FAILURE;
    }
    // Geometry settings
    if (!m_geoTool.retrieve()) {
      error() << "Unable to retrieve the geometry tool!!!" << endmsg;
      return StatusCode::FAILURE;
    }
    // Prepare map of all existing cells in calorimeter to add noise to all
    StatusCode sc_prepareCells = m_geoTool->prepareEmptyCells(m_cellsMap);
    if (sc_prepareCells.isFailure()) {
      error() << "Unable to create empty cells!" << endmsg;
      return StatusCode::FAILURE;
    }
  }
  return StatusCode::SUCCESS;
}

StatusCode CreateCaloCells::execute() {
  // Get the input collection with Geant4 hits
  const edm4hep::SimCalorimeterHitCollection* hits = m_hits.get();
  debug() << "Input Hit collection size: " << hits->size() << endmsg;

  // 0. Clear all cells
  if (m_addCellNoise) {
    std::for_each(m_cellsMap.begin(), m_cellsMap.end(), [](std::pair<const uint64_t, double>& p) { p.second = 0; });
  } else {
    m_cellsMap.clear();
  }

  // 1. Merge energy deposits into cells
  // If running with noise map already was prepared. Otherwise it is being
  // created below
  for (const auto& hit : *hits) {
    verbose() << "CellID : " << hit.getCellID() << endmsg;
    m_cellsMap[hit.getCellID()] += hit.getEnergy();
  }
  debug() << "Number of calorimeter cells after merging of hits: " << m_cellsMap.size() << endmsg;

  // 2. Calibrate simulation energy to EM scale
  if (m_doCellCalibration) {
    m_calibTool->calibrate(m_cellsMap);
  }

  // 3. Add noise to all cells
  if (m_addCellNoise) {
    m_noiseTool->addRandomCellNoise(m_cellsMap);
    if (m_filterCellNoise) {
      m_noiseTool->filterCellNoise(m_cellsMap);
    }
  }

  // 4. Copy information to CaloHitCollection
  edm4hep::CalorimeterHitCollection* edmCellsCollection = new edm4hep::CalorimeterHitCollection();
  for (const auto& cell : m_cellsMap) {
    if (m_addCellNoise || (!m_addCellNoise && cell.second != 0)) {
      edm4hep::CalorimeterHit newCell = edm4hep::CalorimeterHit();
      newCell.setEnergy(cell.second);
      newCell.setCellID(cell.first);
      edmCellsCollection->push_back(newCell);
    }
  }

  // push the CaloHitCollection to event store
  m_cells.put(edmCellsCollection);
  debug() << "Output Cell collection size: " << edmCellsCollection->size() << endmsg;

  return StatusCode::SUCCESS;
}

StatusCode CreateCaloCells::finalize() { return GaudiAlgorithm::finalize(); }
