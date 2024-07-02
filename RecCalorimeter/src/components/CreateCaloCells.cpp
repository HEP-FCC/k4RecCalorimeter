#include "CreateCaloCells.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"

// k4FWCore
#include "k4Interface/IGeoSvc.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DD4hep/Volumes.h"
#include "TGeoManager.h"

// edm4hep
#include "edm4hep/CalorimeterHit.h"

DECLARE_COMPONENT(CreateCaloCells)

CreateCaloCells::CreateCaloCells(const std::string& name, ISvcLocator* svcLoc) :
GaudiAlgorithm(name, svcLoc), m_geoSvc("GeoSvc", name) {
  declareProperty("hits", m_hits, "Hits from which to create cells (input)");
  declareProperty("cells", m_cells, "The created calorimeter cells (output)");

  declareProperty("crosstalksTool", m_crosstalksTool, "Handle for the cell crosstalk tool");
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
  info() << "remove cells below threshold : " << m_filterCellNoise << endmsg;
  info() << "add position information to the cell : " << m_addPosition << endmsg;

  // Initialization of tools
  // Cell crosstalk tool
  if (!m_crosstalksTool.retrieve()) {
    error() << "Unable to retrieve the cell crosstalk tool!!!" << endmsg;
    return StatusCode::FAILURE;
  }
  // Calibrate Geant4 energy to EM scale tool
  if (m_doCellCalibration) {
    if (!m_calibTool.retrieve()) {
      error() << "Unable to retrieve the calo cells calibration tool!!!" << endmsg;
      return StatusCode::FAILURE;
    }
  }
  // Cell noise tool
  if (m_addCellNoise || m_filterCellNoise) {
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
  if (m_addPosition){
    m_volman = m_geoSvc->getDetector()->volumeManager();
  }

  // Copy over the CellIDEncoding string from the input collection to the output collection
  auto hitsEncoding = m_hitsCellIDEncoding.get_optional();
  if (!hitsEncoding.has_value()) {
    error () << "Missing cellID encoding for input collection" << endmsg;
    return StatusCode::FAILURE;
  }
  m_cellsCellIDEncoding.put(hitsEncoding.value());

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

  std::unordered_map<uint64_t, double> HitCellsMap;

  // 1. Merge energy deposits into cells
  // If running with noise map already was prepared. Otherwise it is being
  // created below
  for (const auto& hit : *hits) {
    verbose() << "CellID : " << hit.getCellID() << endmsg;
    if(m_addCrosstalk) { // Crosstalk only applies to energy deposits caused by EM shower hits. Therefore, cell noise needs to be excluded.
      HitCellsMap[hit.getCellID()] += hit.getEnergy();
    }
    m_cellsMap[hit.getCellID()] += hit.getEnergy();
  }
  if(m_addCrosstalk) {
    m_CrosstalkCellsMap.clear();
    // loop over cells that get actual EM shower hits
    for (const auto& this_cell : HitCellsMap) {
      uint64_t this_cellId=this_cell.first;
      auto vec_neighbours = m_crosstalksTool->getNeighbours(this_cellId); // a vector of neighbour IDs
      auto vec_crosstalks = m_crosstalksTool->getCrosstalks(this_cellId); // a vector of crosstalk coefficients
      // loop over crosstalk neighbrous of the cell under study
      for (unsigned int i_cell=0; i_cell<vec_neighbours.size(); i_cell++) {
	// signal transfer = energy deposit brought by EM shower hits * crosstalk coefficient
        double signal_transfer = this_cell.second * vec_crosstalks[i_cell];
	// for the cell under study, record the signal transfer that will be subtracted from its final cell energy
        m_CrosstalkCellsMap[this_cellId] -= signal_transfer;
	// for the crosstalk neighbour, record the signal transfer that will be added to its final cell energy
        m_CrosstalkCellsMap[vec_neighbours[i_cell]] += signal_transfer;
      }
    }

    for (const auto& this_cell : m_CrosstalkCellsMap) {
      m_cellsMap[this_cell.first] += this_cell.second;
    }
    
  }
  debug() << "Number of calorimeter cells after merging of hits: " << m_cellsMap.size() << endmsg;

  // 2. Calibrate simulation energy to EM scale
  if (m_doCellCalibration) {
    m_calibTool->calibrate(m_cellsMap);
  }

  // 3. Add noise to all cells
  if (m_addCellNoise) {
    m_noiseTool->addRandomCellNoise(m_cellsMap);
  }

  // 4. Filter cells
  if (m_filterCellNoise) {
    m_noiseTool->filterCellNoise(m_cellsMap);
  }

  // 5. Copy information to CaloHitCollection
  edm4hep::CalorimeterHitCollection* edmCellsCollection = new edm4hep::CalorimeterHitCollection();
  for (const auto& cell : m_cellsMap) {
    if (m_addCellNoise || (!m_addCellNoise && cell.second != 0)) {
      auto newCell = edmCellsCollection->create();
      newCell.setEnergy(cell.second);
      uint64_t cellid = cell.first;
      newCell.setCellID(cellid);
      if (m_addPosition){
        auto detelement = m_volman.lookupDetElement(cellid);
        const auto& transformMatrix = detelement.nominal().worldTransformation();
        double outGlobal[3];
        double inLocal[] = {0, 0, 0};
        transformMatrix.LocalToMaster(inLocal, outGlobal);
        edm4hep::Vector3f position = edm4hep::Vector3f(outGlobal[0] / dd4hep::mm, outGlobal[1] / dd4hep::mm, outGlobal[2] / dd4hep::mm);
        newCell.setPosition(position);
      }
    }
  }

  // push the CaloHitCollection to event store
  m_cells.put(edmCellsCollection);

  debug() << "Output Cell collection size: " << edmCellsCollection->size() << endmsg;

  return StatusCode::SUCCESS;
}

StatusCode CreateCaloCells::finalize() { return GaudiAlgorithm::finalize(); }
