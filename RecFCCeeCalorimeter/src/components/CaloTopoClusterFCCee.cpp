#include "CaloTopoClusterFCCee.h"

// std
#include <algorithm>
#include <memory>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"

// EDM4hep
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ClusterCollection.h"

// DD4hep
#include "DD4hep/Readout.h"

DECLARE_COMPONENT(CaloTopoClusterFCCee)

CaloTopoClusterFCCee::CaloTopoClusterFCCee(const std::string& name, ISvcLocator* svcLoc)
    : Gaudi::Algorithm(name, svcLoc) {
  declareProperty("noiseTool", m_noiseTool, "Handle for the cells noise tool");
  declareProperty("neigboursTool", m_neighboursTool, "Handle for tool to retrieve cell neighbours");
  declareProperty("clusters", m_clusterCollection, "Handle for calo clusters (output collection)");
  declareProperty("clusterCells", m_clusterCellsCollection, "Handle for clusters (output collection)");
}

StatusCode CaloTopoClusterFCCee::initialize() {

  if (Gaudi::Algorithm::initialize().isFailure()) {
    return StatusCode::FAILURE;
  }

  // create handles for input cell collections
  for (const auto& col : m_cellCollections) {
    debug() << "Creating handle for input cell (CalorimeterHit) collection : " << col << endmsg;
    try {
      m_cellCollectionHandles.push_back(
          new k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection>(col, Gaudi::DataHandle::Reader, this));
    } catch (...) {
      error() << "Error creating handle for input collection: " << col << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // use pre-calculated neighbor map i.e. TTree to retrieve neighbors
  if (m_useNeighborMap) {
    // retrieve cells neighbours tool
    if (!m_neighboursTool.retrieve()) {
      error() << "Unable to retrieve the cells neighbours tool!!!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // use DDSegmentation to retrieve neighbors
  if (!m_useNeighborMap) {
    m_geoSvc = service("GeoSvc");

    if (!m_geoSvc) {
      error() << "Unable to locate Geometry Service. "
              << "Make sure you have GeoSvc in the configuration." << endmsg;

      return StatusCode::FAILURE;
    }

    if (m_geoSvc->getDetector()->readouts().find(m_readoutName) == m_geoSvc->getDetector()->readouts().end()) {
      error() << "Readout <<" << m_readoutName << ">> does not exist." << endmsg;

      return StatusCode::FAILURE;
    }

    // get segmentation
    m_segmentation = m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation();
  }

  // retrieve cells noise tool
  if (!m_noiseTool.retrieve()) {
    error() << "Unable to retrieve the cells noise tool!!!" << endmsg;
    return StatusCode::FAILURE;
  }

  // setup system decoder
  m_decoder = new dd4hep::DDSegmentation::BitFieldCoder(m_systemEncoding);
  m_indexSystem = m_decoder->index("system");

  // initialise the list of metadata for the clusters
  std::vector<std::string> shapeParameterNames = {"dR_over_E"};
  m_shapeParametersHandle.put(shapeParameterNames);

  if (m_createClusterCellCollection) {
    std::vector<int> IDs;
    for (auto ID: m_caloIDs) {
      IDs.push_back(ID);
    }

    std::vector<std::string> colls;
    for (auto coll: m_cellCollections) {
      colls.push_back(coll);
    }

    if (IDs.size()==colls.size()) {
      m_caloIDsMetaData.put(IDs);
      m_cellsMetaData.put(colls);
    } else {
      warning() << "Sizes of input cell and systemID collections of tower tool are different, no metadata written" << endmsg;
    }
  }

  return StatusCode::SUCCESS;
}

StatusCode CaloTopoClusterFCCee::execute(const EventContext&) const {

  // Create output collections
  edm4hep::ClusterCollection* outClusters = m_clusterCollection.createAndPut();
  edm4hep::CalorimeterHitCollection* outClusterCells = nullptr;
  if (m_createClusterCellCollection) {
    outClusterCells = m_clusterCellsCollection.createAndPut();
  }

  // Get input collection with calorimeter cells
  edm4hep::CalorimeterHitCollection* inCells = new edm4hep::CalorimeterHitCollection();
  for (size_t ih = 0; ih < m_cellCollectionHandles.size(); ih++) {
    verbose() << "Processing collection " << ih << endmsg;
    const edm4hep::CalorimeterHitCollection* coll = m_cellCollectionHandles[ih]->get();
    for (const auto& hit : *coll) {
      auto newCell = hit.clone();
      newCell.setType(0); // topoclustering will set the type of clustered cells to 1-2-3 depending whether they are
                          // seed/neighbours/last neighbours
      // note that this overwrites the type information from the digitiser, which encodes calorimeter type / layout /
      // layer
      inCells->push_back(newCell);
    }
  }
  if (inCells->empty()) {
    debug() << "No active cells, skipping event..." << endmsg;
    return StatusCode::SUCCESS;
  }

  debug() << "Number of active cells                               : " << inCells->size() << endmsg;

  // Find seeds
  edm4hep::CalorimeterHitCollection seedCells = findSeeds(inCells);
  debug() << "Number of seeds found                                : " << seedCells.size() << endmsg;

  // Build protoclusters (find neighbouring cells)
  std::map<uint32_t, edm4hep::CalorimeterHitCollection> protoClusters;
  {
    StatusCode sc = buildProtoClusters(seedCells, inCells, protoClusters);
    if (sc.isFailure()) {
      error() << "Unable to build the protoclusters!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // Build clusters
  debug() << "Building " << protoClusters.size() << " clusters" << endmsg;
  double checkTotEnergy = 0.;
  double checkTotEnergyAboveThreshold = 0.;
  int clusterWithMixedCells = 0;
  for (const auto& protoCluster : protoClusters) {

    // calculate cluster energy and decide whether to keep it
    double clusterEnergy = 0.;
    for (const auto& protoCell : protoCluster.second) {
      clusterEnergy += protoCell.getEnergy();
    }
    verbose() << "Cluster energy:     " << clusterEnergy << endmsg;
    checkTotEnergy += clusterEnergy;
    if (clusterEnergy < m_minClusterEnergy) {
      continue;
    }

    // build cluster
    debug() << "Building cluster with ID: " << protoCluster.first << endmsg;
    edm4hep::MutableCluster cluster;

    // set cluster energy
    cluster.setEnergy(clusterEnergy);
    checkTotEnergyAboveThreshold += cluster.getEnergy();

    // loop over the cells attached to the cluster to calculate cluster barycenter and attach cells to cluster
    double clusterPosX = 0.;
    double clusterPosY = 0.;
    double clusterPosZ = 0.;
    double deltaR = 0.;
    std::vector<double> cellPosPhi(protoCluster.second.size(), 0);
    std::vector<double> cellPosTheta(protoCluster.second.size(), 0);
    std::vector<double> cellEnergy(protoCluster.second.size(), 0);
    double sumCellPhi = 0.;
    double sumCellTheta = 0.;
    std::map<int, int> system;
    for (const auto& protoCell : protoCluster.second) {
      // identify calo system
      auto systemId = m_decoder->get(protoCell.getCellID(), m_indexSystem);
      system[int(systemId)]++;
      auto cellPos = dd4hep::Position(protoCell.getPosition().x, protoCell.getPosition().y, protoCell.getPosition().z);

      clusterPosX += protoCell.getPosition().x * protoCell.getEnergy();
      clusterPosY += protoCell.getPosition().y * protoCell.getEnergy();
      clusterPosZ += protoCell.getPosition().z * protoCell.getEnergy();
      cellPosPhi.push_back(cellPos.Phi());
      cellPosTheta.push_back(cellPos.Theta());
      cellEnergy.push_back(protoCell.getEnergy());
      sumCellPhi += cellPos.Phi() * protoCell.getEnergy();
      sumCellTheta += cellPos.Theta() * protoCell.getEnergy();

      if(m_createClusterCellCollection){
        auto cell = protoCell.clone();
        outClusterCells->push_back(cell);
        cluster.addToHits(cell);
      } else {
        for (size_t ih = 0; ih < m_cellCollectionHandles.size(); ih++) {
          const edm4hep::CalorimeterHitCollection* coll = m_cellCollectionHandles[ih]->get();
          for (const auto& hit : *coll) {
            if(hit.id() == protoCell.id()){
              cluster.addToHits(hit);
            }
          }
        }
      }
    }

    // set cluster position (weighted barycentre of cell positions)
    cluster.setPosition(
        edm4hep::Vector3f(clusterPosX / clusterEnergy, clusterPosY / clusterEnergy, clusterPosZ / clusterEnergy));

    // store deltaR of cluster in time for the moment..
    sumCellPhi = sumCellPhi / clusterEnergy;
    sumCellTheta = sumCellTheta / clusterEnergy;
    for (size_t i = 0; i < cellEnergy.size(); ++i) {
      deltaR += std::sqrt(std::pow(cellPosTheta[i] - sumCellTheta, 2) + std::pow(cellPosPhi[i] - sumCellPhi, 2)) *
                cellEnergy[i];
    }
    cluster.addToShapeParameters(deltaR / clusterEnergy);

    outClusters->push_back(cluster);
    if (system.size() > 1)
      clusterWithMixedCells++;

    cellPosPhi.clear();
    cellPosTheta.clear();
    cellEnergy.clear();
  }

  debug() << "Number of clusters with cells in E and HCal:        " << clusterWithMixedCells << endmsg;
  debug() << "Total energy of clusters:                           " << checkTotEnergy << endmsg;
  debug() << "Total energy of clusters above threshold:                           " << checkTotEnergyAboveThreshold
          << endmsg;
  if(m_createClusterCellCollection){
    debug() << "Leftover cells :                                    " << inCells->size() - outClusterCells->size()
            << endmsg;
  }

  delete inCells;
  return StatusCode::SUCCESS;
}

edm4hep::CalorimeterHitCollection
CaloTopoClusterFCCee::findSeeds(const edm4hep::CalorimeterHitCollection* allCells) const {

  std::vector<edm4hep::CalorimeterHit> seedCellsVec;

  for (const auto& cell : *allCells) {

    verbose() << "cellID   = " << cell.getCellID() << endmsg;

    // retrieve the noise const and offset assigned to cell
    double offset = m_noiseTool->getNoiseOffsetPerCell(cell.getCellID());
    double rms = m_noiseTool->getNoiseRMSPerCell(cell.getCellID());
    double threshold = offset + rms * m_seedSigma;

    debug() << "======================================" << endmsg;
    debug() << "noise offset    = " << offset << " GeV " << endmsg;
    debug() << "noise rms       = " << rms << " GeV " << endmsg;
    debug() << "seed threshold  = " << threshold << " GeV " << endmsg;
    debug() << "======================================" << endmsg;
    if (std::fabs(cell.getEnergy()) > threshold) {
      debug() << "Found seed" << endmsg;
      seedCellsVec.emplace_back(cell);
    }
  }

  // Sort the seeds in decending order of their energy
  std::sort(seedCellsVec.begin(), seedCellsVec.end(),
            [](const auto& lhs, const auto& rhs) { return lhs.getEnergy() > rhs.getEnergy(); });

  edm4hep::CalorimeterHitCollection seedCells;
  seedCells.setSubsetCollection();
  for (const auto& cell : seedCellsVec) {
    seedCells.push_back(cell);
  }

  return seedCells;
}

StatusCode
CaloTopoClusterFCCee::buildProtoClusters(const edm4hep::CalorimeterHitCollection& seedCells,
                                         const edm4hep::CalorimeterHitCollection* allCells,
                                         std::map<uint32_t, edm4hep::CalorimeterHitCollection>& protoClusters) const {

  verbose() << "Initial number of seeds to loop over: " << seedCells.size() << endmsg;

  std::map<uint64_t, const edm4hep::CalorimeterHit> allCellsMap;
  for (const auto& cell : *allCells) {
    allCellsMap.emplace(cell.getCellID(), cell);
  }
  std::map<uint64_t, uint32_t> alreadyUsedCells;

  // Loop over every seed in Calo to create first cluster
  uint32_t seedCounter = 0;
  for (const auto& seedCell : seedCells) {
    seedCounter++;
    verbose() << "Looking at seed: " << seedCounter << endmsg;
    auto seedId = seedCell.getCellID();
    auto cellInCluster = alreadyUsedCells.find(seedId);
    if (cellInCluster != alreadyUsedCells.end()) {
      verbose() << "Seed is already assigned to another cluster!" << endmsg;
      continue;
    }

    uint32_t clusterId = seedCounter;
    // new cluster starts with seed
    // set cell type to 1 for seed cell
    edm4hep::MutableCalorimeterHit clusteredCell = seedCell.clone();
    clusteredCell.setType(1);
    protoClusters[clusterId].push_back(clusteredCell);
    alreadyUsedCells[seedId] = clusterId;

    std::vector<std::vector<std::pair<uint64_t, uint32_t>>> nextNeighbours(100);
    nextNeighbours[0] =
        searchForNeighbours(seedId, clusterId, m_neighbourSigma, allCellsMap, alreadyUsedCells, protoClusters, true);

    // first loop over seeds neighbours
    verbose() << "Found " << nextNeighbours[0].size() << " neighbours.." << endmsg;
    int it = 0;
    while (nextNeighbours[it].size() > 0) {
      it++;
      nextNeighbours.emplace_back(std::vector<std::pair<uint64_t, uint>>{});
      for (auto& id : nextNeighbours[it - 1]) {
        if (id.first == 0) {
          error() << "Building of cluster is stopped due to missing cell ID "
                     "in neighbours map!"
                  << endmsg;
          return StatusCode::FAILURE;
        }
        verbose() << "Next neighbours assigned to cluster ID: " << clusterId << endmsg;
        auto additionalNeighbours = searchForNeighbours(id.first, clusterId, m_neighbourSigma, allCellsMap,
                                                        alreadyUsedCells, protoClusters, true);
        nextNeighbours[it].insert(nextNeighbours[it].end(), additionalNeighbours.begin(), additionalNeighbours.end());
      }
      verbose() << "Found " << nextNeighbours[it].size() << " more neighbours.." << endmsg;
    }

    // last try with different condition on neighbours
    if (nextNeighbours[it].size() == 0) {
      // loop over all clustered cells
      for (const auto& cell : protoClusters[clusterId]) {
        if (cell.getType() <= 2) {
          verbose() << "Add neighbours of " << cell.getCellID()
                    << " in last round with thr = " << m_lastNeighbourSigma.value() << " x sigma." << endmsg;
          auto lastNeighours = searchForNeighbours(cell.getCellID(), clusterId, m_lastNeighbourSigma, allCellsMap,
                                                   alreadyUsedCells, protoClusters, false);
        }
      }
    }
  }

  return StatusCode::SUCCESS;
}

std::vector<std::pair<uint64_t, uint32_t>> CaloTopoClusterFCCee::searchForNeighbours(
    const uint64_t aCellId, uint& aClusterID, int aNumSigma,
    std::map<uint64_t, const edm4hep::CalorimeterHit>& allCellsMap, std::map<uint64_t, uint32_t>& alreadyUsedCells,
    std::map<uint32_t, edm4hep::CalorimeterHitCollection>& protoClusters, bool allowClusterMerge) const {

  // Fill vector to be returned, next cell ids and cluster id for which
  // neighbours are found
  std::vector<std::pair<uint64_t, uint32_t>> additionalNeighbours;
  std::vector<uint64_t> neighboursVec;

  // Retrieve cellIDs of neighbours
  if (m_useNeighborMap) {
    neighboursVec = m_neighboursTool->neighbours(aCellId);

    if (neighboursVec.size() == 0) {
      error() << "No neighbours for cellID found! " << endmsg;
      error() << "to cellID :  " << aCellId << endmsg;
      error() << "in system:   " << m_decoder->get(aCellId, m_indexSystem) << endmsg;
      additionalNeighbours.resize(0);
      additionalNeighbours.push_back(std::make_pair(0, 0));

      return additionalNeighbours;
    }
  }

  if (!m_useNeighborMap) {
    // DDSegmentation returns std::set
    std::set<dd4hep::DDSegmentation::CellID> outputNeighbors;
    m_segmentation->neighbours(aCellId, outputNeighbors);
    neighboursVec = std::vector<uint64_t>(outputNeighbors.begin(), outputNeighbors.end());
  }

  verbose() << "For cluster: " << aClusterID << endmsg;
  // loop over neighbours
  for (const auto& neighbourID : neighboursVec) {
    // Find the neighbour in the Calo cells list
    auto itAllCells = allCellsMap.find(neighbourID);
    auto itAllUsedCells = alreadyUsedCells.find(neighbourID);

    // If cell is hit.. and is not assigned to a cluster
    if (itAllCells != allCellsMap.end() && itAllUsedCells == alreadyUsedCells.end()) {
      verbose() << "Found neighbour with CellID: " << neighbourID << endmsg;
      auto neighbouringCellEnergy = allCellsMap[neighbourID].getEnergy();
      bool addNeighbour = false;
      int cellType = 2;
      // retrieve the cell noise level [GeV]
      double offset = m_noiseTool->getNoiseOffsetPerCell(neighbourID);
      double rms = m_noiseTool->getNoiseRMSPerCell(neighbourID);
      double thr = offset + rms * aNumSigma;
      if (std::fabs(neighbouringCellEnergy) > thr)
        addNeighbour = true;
      else
        addNeighbour = false;

      // give cell type according to threshold
      if (aNumSigma == m_lastNeighbourSigma) {
        cellType = 3;
      }
      // if threshold is 0, collect the cell independent on its energy
      if (aNumSigma == 0) {
        addNeighbour = true;
      }
      // if neighbour is validated
      if (addNeighbour) {
        // retrieve the cell
        // add neighbour to cells for cluster
        edm4hep::MutableCalorimeterHit clusteredCell = allCellsMap[neighbourID].clone();
        clusteredCell.setType(cellType);
        protoClusters[aClusterID].push_back(clusteredCell);
        alreadyUsedCells[neighbourID] = aClusterID;
        additionalNeighbours.push_back(std::make_pair(neighbourID, aClusterID));
      }
    }
    // If cell is hit.. but is assigned to another cluster
    else if (itAllUsedCells != alreadyUsedCells.end() && itAllUsedCells->second != aClusterID && allowClusterMerge) {
      uint32_t clusterIDToMergeTo = itAllUsedCells->second;
      if (msgLevel() <= MSG::VERBOSE) {
        verbose() << "This neighbour was found in cluster " << clusterIDToMergeTo << ", cluster " << aClusterID
                  << " will be merged!" << endmsg;
        verbose() << "Assigning all cells ( " << protoClusters[aClusterID].size() << " ) to Cluster "
                  << clusterIDToMergeTo << " with ( " << protoClusters[clusterIDToMergeTo].size() << " ). " << endmsg;
      }
      // Fill all cells into cluster, and assigned cells to new cluster
      alreadyUsedCells[neighbourID] = clusterIDToMergeTo;
      for (const auto& cell : protoClusters[aClusterID]) {
        alreadyUsedCells[cell.getCellID()] = clusterIDToMergeTo;
        // make sure that already assigned cells are not added
        if (cellIdInColl(cell.getCellID(), protoClusters[clusterIDToMergeTo])) {
          continue;
        }
        protoClusters[clusterIDToMergeTo].push_back(cell.clone());
      }
      protoClusters.erase(aClusterID);
      // changed clusterId -> if more neighbours are found, correct assignment
      verbose() << "Cluster Id changed to " << clusterIDToMergeTo << endmsg;
      aClusterID = clusterIDToMergeTo;
      // found neighbour for next search
      additionalNeighbours.push_back(std::make_pair(neighbourID, aClusterID));
      // end loop to ensure correct cluster assignment
      break;
    }
  }

  return additionalNeighbours;
}

StatusCode CaloTopoClusterFCCee::finalize() {
  delete m_decoder;
  for (size_t ih = 0; ih < m_cellCollectionHandles.size(); ih++)
    delete m_cellCollectionHandles[ih];

  return Gaudi::Algorithm::finalize();
}

inline bool CaloTopoClusterFCCee::cellIdInColl(const uint64_t cellId,
                                               const edm4hep::CalorimeterHitCollection& coll) const {
  for (const auto& cell : coll) {
    if (cell.getCellID() == cellId) {
      return true;
    }
  }
  return false;
}
