#include "CreatePositionedCaloCells.h"

// dd4hep
#include "DD4hep/Detector.h"
#include "DD4hep/DetType.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"

// edm4hep
#include "edm4hep/CalorimeterHit.h"

DECLARE_COMPONENT(CreatePositionedCaloCells)

CreatePositionedCaloCells::CreatePositionedCaloCells(const std::string& name, ISvcLocator* svcLoc) :
Gaudi::Algorithm(name, svcLoc) {
  declareProperty("hits", m_hits, "Hits from which to create cells (input)");
  declareProperty("cells", m_cells, "The created calorimeter cells (output)");
  declareProperty("links", m_links, "The links between hits and cells (output)");

  declareProperty("positionsTool", m_cellPositionsTool, "Handle for cell positions tool");
  declareProperty("crosstalkTool", m_crosstalkTool, "Handle for the cell crosstalk tool");
  declareProperty("calibTool", m_calibTool, "Handle for tool to calibrate Geant4 energy to EM scale tool");
  declareProperty("noiseTool", m_noiseTool, "Handle for the calorimeter cells noise tool");
  declareProperty("geometryTool", m_geoTool, "Handle for the geometry tool");

  m_decoder = nullptr;
}

CreatePositionedCaloCells::~CreatePositionedCaloCells() {
  delete m_decoder;
}

StatusCode CreatePositionedCaloCells::initialize() {
  StatusCode sc = Gaudi::Algorithm::initialize();
  if (sc.isFailure()) return sc;

  info() << "CreatePositionedCaloCells initialized" << endmsg;
  info() << "do calibration : " << m_doCellCalibration << endmsg;
  info() << "add cell noise : " << m_addCellNoise << endmsg;
  info() << "remove cells below threshold : " << m_filterCellNoise << endmsg;
  info() << "emulate crosstalk : " << m_addCrosstalk << endmsg;

  // Initialization of tools

  // Cell position tool
  if (!m_cellPositionsTool.retrieve()) {
    error() << "Unable to retrieve the cell positions tool!!!" << endmsg;
    return StatusCode::FAILURE;
  }
  // Cell crosstalk tool
  if (m_addCrosstalk) {
    if (!m_crosstalkTool.retrieve()) {
      error() << "Unable to retrieve the cell crosstalk tool!!!" << endmsg;
      return StatusCode::FAILURE;
    }
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
  }
  if (m_addCellNoise) {
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
    verbose() << "Initialised empty cell map with size " << m_cellsMap.size() << endmsg;
    // noise filtering erases cells from the cell map after each event, so we need
    // to backup the empty cell map for later reuse
    if (m_filterCellNoise) {
      m_emptyCellsMap = m_cellsMap;
    }
  }

  // Copy over the CellIDEncoding string from the input collection to the output collection
  auto hitsEncoding = m_hitsCellIDEncoding.get_optional();
  if (!hitsEncoding.has_value()) {
    error () << "Missing cellID encoding for input collection" << endmsg;
    return StatusCode::FAILURE;
  }
  m_cellsCellIDEncoding.put(hitsEncoding.value());
  m_decoder = new dd4hep::DDSegmentation::BitFieldCoder(hitsEncoding.value());

  // these variables will be initilized in the execute() method
  // once the system ID is extracted from one cell
  m_calotype = -99; // -99 not initialized, 1 unknown, 0 em, 1 had, 2 muon
  m_caloid = 0 ;  // 0 unknown, 1 ecal, 2 hcal, 3 yoke
  m_layout = 0 ; // 0 any, 1 barrel, 2 endcap

  return StatusCode::SUCCESS;
}

StatusCode CreatePositionedCaloCells::execute(const EventContext&) const {
  // Get the input collection with Geant4 hits
  const edm4hep::SimCalorimeterHitCollection* hits = m_hits.get();
  debug() << "Input Hit collection size: " << hits->size() << endmsg;

  // 0. Clear all cells
  if (m_addCellNoise) {
    // if cells are not filtered, the map has same size in each event, equal to the total number
    // of cells in the calorimeter, so we can just reset the values to 0
    // if cells are filtered, during each event they are removed from the cellsMap, so one has to
    // restore the initial map of all empty cells
    if (!m_filterCellNoise)
      std::for_each(m_cellsMap.begin(), m_cellsMap.end(), [](std::pair<const uint64_t, double>& p) { p.second = 0; });
    else
      m_cellsMap = m_emptyCellsMap;
  } else {
    m_cellsMap.clear();
  }

  // 1. Merge energy deposits into cells
  // If running with noise, map was already prepared in initialize().
  // Otherwise it is being created below
  for (const auto& hit : *hits) {
    auto id = hit.getCellID();
    verbose() << "CellID : " << id << endmsg;
    m_cellsMap[id] += hit.getEnergy();
  }
  debug() << "Number of calorimeter cells after merging of hits: " << m_cellsMap.size() << endmsg;

  // 2. Emulate cross-talk (if asked)
  if (m_addCrosstalk) {
    // Derive the cross-talk contributions without affecting yet the nominal energy
    // (one has to emulate crosstalk based on cells free from any cross-talk contributions)
    m_crosstalkCellsMap.clear(); // this is a temporary map to hold energy exchange due to cross-talk, without affecting yet the nominal energy
    // loop over cells with nominal energies
    for (const auto& this_cell : m_cellsMap) {
      uint64_t this_cellId = this_cell.first;
      auto vec_neighbours = m_crosstalkTool->getNeighbours(this_cellId); // a vector of neighbour IDs
      auto vec_crosstalks = m_crosstalkTool->getCrosstalks(this_cellId); // a vector of crosstalk coefficients
      // loop over crosstalk neighbours of the cell under study
      for (unsigned int i_cell=0; i_cell<vec_neighbours.size(); i_cell++) {
        // signal transfer = energy deposit brought by EM shower hits * crosstalk coefficient
        double signal_transfer = this_cell.second * vec_crosstalks[i_cell];
        // for the cell under study, record the signal transfer that will be subtracted from its final cell energy
        m_crosstalkCellsMap[this_cellId] -= signal_transfer;
        // for the crosstalk neighbour, record the signal transfer that will be added to its final cell energy
        m_crosstalkCellsMap[vec_neighbours[i_cell]] += signal_transfer;
      }
    }

    // apply the cross-talk contributions on the nominal cell-energy map
    for (const auto& this_cell : m_crosstalkCellsMap) {
      m_cellsMap[this_cell.first] += this_cell.second;
    }
    
  }

  // 3. Calibrate simulation energy to EM scale
  if (m_doCellCalibration) {
    m_calibTool->calibrate(m_cellsMap);
  }

  // 4. Add noise to all cells
  if (m_addCellNoise) {
    m_noiseTool->addRandomCellNoise(m_cellsMap);
  }

  // 5. Filter cells
  if (m_filterCellNoise) {
    m_noiseTool->filterCellNoise(m_cellsMap);
  }

  // determine detector type (only once)
  if (m_calotype==-99 && m_cellsMap.size()>0) {
    info() << "Determining calorimeter type for input collection " << m_hits.objKey() << endmsg;
    uint cellid = m_cellsMap.begin()->first;
    int system = m_decoder->get(cellid, "system");
    debug() << "System: " << system << endmsg;
    dd4hep::Detector* dd4hepgeo = &(dd4hep::Detector::getInstance());
    const dd4hep::Detector::HandleMap& children = dd4hepgeo->detectors();
    for (dd4hep::Detector::HandleMap::const_iterator i=children.begin(); i!=children.end(); ++i) {
      dd4hep::DetElement det((*i).second);
      int id = det.id();
      if (id==system) {
        dd4hep::DetType detType(det.typeFlag());
        if ( detType.is( dd4hep::DetType::CALORIMETER ) ) {
          if ( detType.is( dd4hep::DetType::ELECTROMAGNETIC ) ) {
            m_calotype = 0;
            m_caloid = 1;
          }
          else if ( detType.is( dd4hep::DetType::HADRONIC ) ) {
            m_calotype = 1;
            m_caloid = 2;
          }
          else if ( detType.is( dd4hep::DetType::MUON ) ) {
            m_calotype = 2;
            m_caloid = 3;
          }
          else {
            warning() << "Detector type is neither ELECTROMAGNETIC, HADRONIC nor MUON" << endmsg;
          }
          if ( detType.is( dd4hep::DetType::BARREL ) ) {
            m_layout = 1;
          }
          else if ( detType.is( dd4hep::DetType::ENDCAP ) ) {
            m_layout = 2;
          }
          else {
            warning() << "Detector type is neither BARREL nor ENDCAP" << endmsg;
          }
        }
        else {
          warning() << "Detector type is not CALORIMETER" << endmsg;
          m_calotype = -1;
          m_caloid = 0;
        }
        break;
      }
    }
    info() << "CaloType: " << m_calotype << endmsg;
    info() << "CaloId: " << m_caloid << endmsg;
    info() << "Layout: " << m_layout << endmsg;
  }

  // 6. Copy information to CaloHitCollection
  edm4hep::CalorimeterHitCollection* edmCellsCollection = new edm4hep::CalorimeterHitCollection();
  for (const auto& cell : m_cellsMap) {
    if (m_addCellNoise || (!m_addCellNoise && cell.second != 0.)) {
      auto newCell = edmCellsCollection->create();
      newCell.setEnergy(cell.second);
      uint64_t cellid = cell.first;
      newCell.setCellID(cellid);

      // add cell position
      auto cached_pos = m_positions_cache.find(cellid);
      if(cached_pos == m_positions_cache.end()) {
        // retrieve position from tool
        dd4hep::Position posCell = m_cellPositionsTool->xyzPosition(cellid);
        edm4hep::Vector3f edmPos;
        edmPos.x = posCell.x() / dd4hep::mm;
        edmPos.y = posCell.y() / dd4hep::mm;
        edmPos.z = posCell.z() / dd4hep::mm;
        m_positions_cache[cellid] = edmPos;
        newCell.setPosition(edmPos);
      }
      else {
        newCell.setPosition(cached_pos->second);
      }

      // add cell type (for Pandora) - see iLCSoft/MarlinUtil/source/include/CalorimeterHitType.h
      int layer = m_decoder->get(cellid, "layer");
      newCell.setType(m_calotype + 10*m_caloid + 1000*m_layout + 10000*layer);

      debug() << "Cell energy (GeV) : " << newCell.getEnergy() << "\tcellID " << newCell.getCellID()  << "\tcellType " << newCell.getType() << endmsg;
      debug() << "Position of cell (mm) : \t" << newCell.getPosition().x
                                      << "\t" << newCell.getPosition().y
                                      << "\t" << newCell.getPosition().z << endmsg;
    }
  }

  // create hits<->cell links
  edm4hep::CaloHitSimCaloHitLinkCollection* edmCellHitLinksCollection = new edm4hep::CaloHitSimCaloHitLinkCollection();
  for (const auto& cell : *edmCellsCollection) {
    auto cellID = cell.getCellID();
    for (const auto& hit : *hits) {
      auto hitID = hit.getCellID();
      if (hitID == cellID) {
        // create Sim<->Reco hit associations
        auto link = edmCellHitLinksCollection->create();
        link.setFrom(cell);
        link.setTo(hit);
      }
    }
  }

  // push the CaloHitCollection to event store
  m_cells.put(edmCellsCollection);
  m_links.put(edmCellHitLinksCollection);

  debug() << "Output Cell collection size: " << edmCellsCollection->size() << endmsg;

  return StatusCode::SUCCESS;
}

StatusCode CreatePositionedCaloCells::finalize() { return Gaudi::Algorithm::finalize(); }
