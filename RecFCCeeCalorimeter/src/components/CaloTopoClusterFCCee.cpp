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

#include "k4FWCore/MetadataUtils.h"

#include "RecCaloCommon/phihelper.h"

// EDM4hep
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/Constants.h"

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
  k4FWCore::putCollectionParameter(m_clusterCollection.objKey(), edm4hep::labels::ShapeParameterNames,
                                   shapeParameterNames, this);

  if (m_createClusterCellCollection) {
    std::vector<int> IDs;
    for (auto ID : m_caloIDs) {
      IDs.push_back(ID);
    }

    std::vector<std::string> colls;
    for (auto coll : m_cellCollections) {
      colls.push_back(coll);
    }

    if (IDs.size() == colls.size()) {
      k4FWCore::putCollectionParameter(m_clusterCollection.objKey(), "inputSystemIDs", IDs, this);
      k4FWCore::putCollectionParameter(m_clusterCollection.objKey(), "inputCellCollections", colls, this);
    } else {
      warning() << "Sizes of input cell and systemID collections of tower tool are different, no metadata written"
                << endmsg;
    }
  }

  return StatusCode::SUCCESS;
}

StatusCode CaloTopoClusterFCCee::execute(const EventContext&) const {

  // create output collections
  edm4hep::ClusterCollection* outClusters = m_clusterCollection.createAndPut();
  edm4hep::CalorimeterHitCollection* outClusterCells = nullptr;
  if (m_createClusterCellCollection) {
    outClusterCells = m_clusterCellsCollection.createAndPut();
  }

  // get input collection with calorimeter cells and build cell cache and flat cell map
  std::unordered_map<uint64_t, const edm4hep::CalorimeterHit> allCellsMap;
  allCellsMap.reserve(2000000);
  std::unordered_map<uint64_t, FastCell> allCells;
  allCells.reserve(2000000);

  for (auto& hdl : m_cellCollectionHandles) {
    const auto* coll = hdl->get();
    for (const auto& hit : *coll) {
      // cache EDM hit
      const uint64_t cID = hit.getCellID();
      allCellsMap.emplace(cID, hit);

      // create fast flat cell
      float energy = hit.getEnergy();
      auto pos = hit.getPosition();
      auto [rms, offset] = m_noiseTool->getNoisePerCell(cID);
      float sovern = (rms > 0.) ? (std::fabs(energy - offset) / rms) : 999999.;
      allCells.emplace(cID, FastCell{cID, energy, (float)pos.x, (float)pos.y, (float)pos.z, 0, sovern});
    }
  }

  // skip event if no cells to cluster
  if (allCells.empty()) {
    debug() << "No active cells, skipping event..." << endmsg;
    return StatusCode::SUCCESS;
  }
  debug() << "Number of active cells                               : " << allCells.size() << endmsg;

  // find seeds (cells with S/N > seedSigma)
  // and sort by energy in reversed order
  std::vector<FastCell> seedCellsVec;
  seedCellsVec.reserve(allCells.size() / 10);
  for (const auto& [cID, cell] : allCells) {
    if (msgLevel() <= MSG::VERBOSE)
      verbose() << "cellID   = " << cID << endmsg;
    if (cell.SoverN > m_seedSigma) {
      if (msgLevel() <= MSG::VERBOSE)
        verbose() << "Found seed" << endmsg;
      seedCellsVec.push_back(cell);
    }
  }
  std::sort(seedCellsVec.begin(), seedCellsVec.end(),
            [](const FastCell& a, const FastCell& b) { return a.energy > b.energy; });

  debug() << "Number of seeds found                                : " << seedCellsVec.size() << endmsg;

  // build clusters (find neighbouring cells)
  // cluster maps clusterID to a FastCluster (vector of FastCells)
  debug() << "Building clusters" << endmsg;

  FastClusterMap clusters;
  StatusCode sc = buildClusters(seedCellsVec, allCells, clusters);

  if (sc.isFailure()) {
    error() << "Unable to build the clusters!" << endmsg;
    return StatusCode::FAILURE;
  }

  // keep only clusters with sufficient energy and build EDM output clusters
  debug() << "Building EDM clusters from " << clusters.size() << " clusters" << endmsg;

  double checkTotEnergy = 0.;
  double checkTotEnergyAboveThreshold = 0.;
  int clusterWithMixedCells = 0;

  for (const auto& [clusterId, cluster] : clusters) {

    double clusterEnergy = 0.;
    std::unordered_map<int, int> system;
    system.reserve(4);

    for (const auto& fastcell : cluster) {

      clusterEnergy += fastcell.energy;
      auto systemId = m_decoder->get(fastcell.cellID, m_indexSystem);
      system[int(systemId)]++;
    }

    checkTotEnergy += clusterEnergy;
    if (msgLevel() <= MSG::VERBOSE)
      verbose() << "Cluster energy: " << clusterEnergy << endmsg;
    if (clusterEnergy < m_minClusterEnergy) {
      continue;
    }

    // build cluster
    debug() << "Building cluster with ID: " << clusterId << endmsg;
    edm4hep::MutableCluster outCluster;

    // set cluster energy
    outCluster.setEnergy(clusterEnergy);
    checkTotEnergyAboveThreshold += clusterEnergy;

    // loop over the cells attached to the cluster to calculate cluster barycenter and attach cells to cluster
    double clusterPosX = 0.;
    double clusterPosY = 0.;
    double clusterPosZ = 0.;

    double sumCellPhi = 0.;
    double sumCellTheta = 0.;
    double phi0 = 0;

    double deltaR = 0.;

    std::vector<double> cellPhi;
    std::vector<double> cellTheta;
    std::vector<double> cellEnergy;

    cellPhi.reserve(cluster.size());
    cellTheta.reserve(cluster.size());
    cellEnergy.reserve(cluster.size());

    using k4::recCalo::deltaPhi, k4::recCalo::wrapToPi;

    for (const auto& fastcell : cluster) {

      const auto& cell = allCellsMap.at(fastcell.cellID);

      double energy = fastcell.energy;

      clusterPosX += fastcell.x * energy;
      clusterPosY += fastcell.y * energy;
      clusterPosZ += fastcell.z * energy;

      double phi = std::atan2(fastcell.y, fastcell.x);
      double theta = std::atan2(std::hypot(fastcell.x, fastcell.y), fastcell.z);

      cellPhi.push_back(phi);
      cellTheta.push_back(theta);
      cellEnergy.push_back(energy);

      if (sumCellPhi == 0)
        phi0 = phi;
      sumCellPhi += deltaPhi(phi, phi0) * energy;
      sumCellTheta += theta * energy;

      // attach cell
      if (m_createClusterCellCollection) {

        auto newcell = cell.clone();
        newcell.setType(fastcell.type);
        outClusterCells->push_back(newcell);
        outCluster.addToHits(newcell);

      } else {
        outCluster.addToHits(cell);
      }
    }

    // calculate cluster barycenter
    if (clusterEnergy > 0.0) {
      outCluster.setPosition(
          edm4hep::Vector3f(clusterPosX / clusterEnergy, clusterPosY / clusterEnergy, clusterPosZ / clusterEnergy));

      sumCellPhi = wrapToPi(sumCellPhi / clusterEnergy + phi0);
      sumCellTheta /= clusterEnergy;

      for (size_t i = 0; i < cellEnergy.size(); ++i) {
        deltaR += std::sqrt(std::pow(cellTheta[i] - sumCellTheta, 2) + std::pow(deltaPhi(cellPhi[i], sumCellPhi), 2)) *
                  cellEnergy[i];
      }
      outCluster.addToShapeParameters(deltaR / clusterEnergy);
    } else {
      outCluster.setPosition(edm4hep::Vector3f(0., 0., 0.));
      outCluster.addToShapeParameters(0.0);
    }

    outClusters->push_back(outCluster);

    if (system.size() > 1)
      clusterWithMixedCells++;
  }

  if (msgLevel() <= MSG::DEBUG) {
    debug() << "Number of clusters:                                 " << outClusters->size() << endmsg;
    debug() << "Number of clusters with cells in multiple systems:  " << clusterWithMixedCells << endmsg;
    debug() << "Total energy of clusters:                           " << checkTotEnergy << endmsg;
    debug() << "Total energy of clusters above threshold:           " << checkTotEnergyAboveThreshold << endmsg;
    if (m_createClusterCellCollection) {
      debug() << "Leftover cells :                                    " << allCells.size() - outClusterCells->size()
              << endmsg;
    }
  }

  return StatusCode::SUCCESS;
}

StatusCode CaloTopoClusterFCCee::buildClusters(const std::vector<FastCell>& seedCells,
                                               const std::unordered_map<uint64_t, FastCell>& allCells,
                                               FastClusterMap& clusters) const {

  if (msgLevel() <= MSG::VERBOSE)
    verbose() << "Initial number of seeds to loop over: " << seedCells.size() << endmsg;

  std::unordered_map<uint64_t, uint32_t> usedCells;
  usedCells.reserve(allCells.size());

  std::unordered_map<uint32_t, std::unordered_set<uint64_t>> clusterMembers;
  clusterMembers.reserve(seedCells.size());

  // loop over every seeds in calo to build a cluster (or merge with another cluster if appropriate)
  uint32_t seedCounter = 0;
  for (const auto& seedCell : seedCells) {
    seedCounter++;
    if (msgLevel() <= MSG::VERBOSE)
      verbose() << "Looking at seed: " << seedCounter << endmsg;

    auto seedId = seedCell.cellID;
    if (usedCells.find(seedId) != usedCells.end()) {
      if (msgLevel() <= MSG::VERBOSE)
        verbose() << "Seed already assigned to another cluster" << endmsg;
      continue;
    }

    uint32_t clusterId = seedCounter;

    // seed insertion (type = 1)
    auto& cluster = clusters[clusterId];
    cluster.reserve(128);
    cluster.push_back(seedCell);
    cluster.back().type = 1;
    usedCells[seedId] = clusterId;
    clusterMembers[clusterId].insert(seedId);

    std::vector<std::vector<uint64_t>> nextNeighbours(100);
    nextNeighbours[0] =
        searchForNeighbours(seedId, clusterId, m_neighbourSigma, allCells, usedCells, clusters, clusterMembers, true);
    if (msgLevel() <= MSG::VERBOSE)
      verbose() << "Found " << nextNeighbours[0].size() << " neighbours.." << endmsg;

    // first loop over seeds neighbours
    int it = 0;
    while (nextNeighbours[it].size() > 0) {
      it++;
      if (msgLevel() <= MSG::VERBOSE)
        verbose() << "it: " << it << endmsg;
      nextNeighbours.emplace_back(std::vector<uint64_t>{});
      for (auto& id : nextNeighbours[it - 1]) {
        if (id == 0) {
          error() << "Building of cluster is stopped due to missing cell ID "
                     "in neighbours map!"
                  << endmsg;
          return StatusCode::FAILURE;
        }
        if (msgLevel() <= MSG::VERBOSE)
          verbose() << "Next neighbours assigned to cluster ID: " << clusterId << endmsg;
        auto additionalNeighbours =
            searchForNeighbours(id, clusterId, m_neighbourSigma, allCells, usedCells, clusters, clusterMembers, true);
        nextNeighbours[it].insert(nextNeighbours[it].end(), additionalNeighbours.begin(), additionalNeighbours.end());
      }
      if (msgLevel() <= MSG::VERBOSE)
        verbose() << "Found " << nextNeighbours[it].size() << " more neighbours.." << endmsg;
    }

    // last try with different condition on neighbours
    if (nextNeighbours[it].size() == 0) {
      // loop over all clustered cells
      auto& aCluster = clusters[clusterId];
      for (size_t i = 0; i < aCluster.size(); ++i) {
        const auto& cell = aCluster[i];
        if (cell.type <= 2) {
          uint64_t cID = cell.cellID;
          if (msgLevel() <= MSG::VERBOSE)
            verbose() << "Add neighbours of " << cID << " in last round with thr = " << m_lastNeighbourSigma.value()
                      << " x sigma." << endmsg;
          auto lastNeighbours = searchForNeighbours(cID, clusterId, m_lastNeighbourSigma, allCells, usedCells, clusters,
                                                    clusterMembers, false);
        }
      }
    }
  }

  return StatusCode::SUCCESS;
}

std::vector<uint64_t> CaloTopoClusterFCCee::searchForNeighbours(
    const uint64_t cellID, uint& clusterID, int nSigma, const std::unordered_map<uint64_t, FastCell>& allCells,
    std::unordered_map<uint64_t, uint32_t>& usedCells, FastClusterMap& clusters,
    std::unordered_map<uint32_t, std::unordered_set<uint64_t>>& clusterMembers, bool allowClusterMerge) const {

  std::vector<uint64_t> additionalNeighbours;

  assert(allCells.find(cellID) != allCells.end());

  // retrieve neighbours
  std::vector<uint64_t> neighboursVec;
  if (m_useNeighborMap) {

    neighboursVec = m_neighboursTool->neighbours(cellID);

  } else {

    std::set<dd4hep::DDSegmentation::CellID> outputNeighbors;
    m_segmentation->neighbours(cellID, outputNeighbors);
    neighboursVec.assign(outputNeighbors.begin(), outputNeighbors.end());
  }

  if (neighboursVec.empty()) {
    error() << "No neighbours for cellID " << cellID << endmsg;
    return {0};
  }

  // loop over neighbours
  if (msgLevel() <= MSG::VERBOSE)
    verbose() << "For cluster: " << clusterID << " , cell " << cellID << endmsg;
  for (const auto& neighbourID : neighboursVec) {

    auto itCell = allCells.find(neighbourID);
    auto itUsed = usedCells.find(neighbourID);

    // CASE 1: unused cell -> candidate addition
    if (itCell != allCells.end() && itUsed == usedCells.end()) {
      if (msgLevel() <= MSG::VERBOSE)
        verbose() << "Found neighbour with CellID: " << neighbourID << endmsg;

      // const auto& hit = *(itCell->second);
      const auto& hit = (itCell->second);
      bool addNeighbour = (hit.SoverN > nSigma) || (nSigma == 0);

      if (addNeighbour) {
        if (msgLevel() <= MSG::VERBOSE)
          verbose() << "Neighbour kept, hit = " << hit.cellID << endmsg;
        int cellType = (nSigma == m_lastNeighbourSigma) ? 3 : 2;
        auto& cluster = clusters[clusterID];
        cluster.push_back(hit);
        cluster.back().type = cellType;
        usedCells[neighbourID] = clusterID;
        clusterMembers[clusterID].insert(neighbourID);
        additionalNeighbours.emplace_back(neighbourID);
      } else {
        if (msgLevel() <= MSG::VERBOSE)
          verbose() << "Neighbour NOT kept, hit = " << hit.cellID << endmsg;
      }
    }

    // CASE 2: already used -> possible merge
    else if (itUsed != usedCells.end() && itUsed->second != clusterID && allowClusterMerge) {

      uint32_t targetCluster = itUsed->second;

      auto& src = clusters[clusterID];
      auto& dst = clusters[targetCluster];
      auto& dstMembers = clusterMembers[targetCluster];

      if (msgLevel() <= MSG::VERBOSE) {
        verbose() << "Neighbour " << neighbourID << " was found in cluster " << targetCluster << ", cluster "
                  << clusterID << " will be merged!" << endmsg;
        verbose() << "Assigning all cells ( " << clusters[clusterID].size() << " ) to Cluster " << targetCluster
                  << " with ( " << clusters[targetCluster].size() << " ). " << endmsg;
      }

      // merge all cells
      for (const auto& c : src) {

        usedCells[c.cellID] = targetCluster;

        if (dstMembers.insert(c.cellID).second) {
          dst.push_back(c);
        }
      }

      clusters.erase(clusterID);
      clusterMembers.erase(clusterID);
      // changed clusterId -> if more neighbours are found, correct assignment
      clusterID = targetCluster;
      // found neighbour for next search
      additionalNeighbours.emplace_back(neighbourID);
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
