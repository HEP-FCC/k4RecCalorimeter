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

  // get input collection with calorimeter cells and build flat cell cache
  m_cellCache.clear();
  m_cellCache.reserve(2000000);
  std::vector<FastCell> allCells;
  allCells.reserve(2000000);

  for (size_t ih = 0; ih < m_cellCollectionHandles.size(); ih++) {

    const auto* coll = m_cellCollectionHandles[ih]->get();
    for (const auto& hit : *coll) {

      // cache EDM hit
      const uint64_t cID = hit.getCellID();
      m_cellCache.emplace(cID, hit);

      // create fast flat cell
      float energy = hit.getEnergy();
      auto pos = hit.getPosition();
      auto [rms, offset] = m_noiseTool->getNoisePerCell(cID);
      float sovern = (rms > 0.) ? (std::fabs(energy - offset) / rms) : 999999.;
      allCells.emplace_back(FastCell{cID, energy, (float)pos.x, (float)pos.y, (float)pos.z, 0, sovern});
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

  for (const auto& cell : allCells) {
    if (cell.SoverN > m_seedSigma) {
      seedCellsVec.push_back(cell);
    }
  }

  std::sort(seedCellsVec.begin(), seedCellsVec.end(),
            [](const FastCell& a, const FastCell& b) { return a.energy > b.energy; });

  debug() << "Number of seeds found                                : " << seedCellsVec.size() << endmsg;

  // build clusters (find neighbouring cells)
  FastClusterMap clusters;
  clusters.reserve(seedCellsVec.size());

  StatusCode sc = buildClusters(seedCellsVec, allCells, clusters);
  if (sc.isFailure()) {
    error() << "Unable to build clusters!" << endmsg;
    return StatusCode::FAILURE;
  }

  // keep only clusters with sufficient energy and build EDM output clusters
  debug() << "Building " << clusters.size() << " EDM clusters" << endmsg;
  double checkTotEnergy = 0.;
  double checkTotEnergyAboveThreshold = 0.;
  int clusterWithMixedCells = 0;

  for (auto& [clusterId, cluster] : clusters) {

    double clusterEnergy = 0.;
    for (const auto& c : cluster) {
      clusterEnergy += c.energy;
    }
    checkTotEnergy += clusterEnergy;
    verbose() << "Cluster energy:     " << clusterEnergy << endmsg;
    if (clusterEnergy < m_minClusterEnergy) {
      continue;
    }

    // build cluster
    debug() << "Building cluster with ID: " << clusterId << endmsg;
    edm4hep::MutableCluster outCluster;

    // set cluster energy
    outCluster.setEnergy(clusterEnergy);

    // loop over the cells attached to the cluster to calculate cluster barycenter and attach cells to cluster
    double cx = 0, cy = 0, cz = 0;

    std::vector<double> phiVec;
    std::vector<double> thetaVec;
    std::vector<double> eVec;

    phiVec.reserve(cluster.size());
    thetaVec.reserve(cluster.size());
    eVec.reserve(cluster.size());

    double sumPhi = 0;
    double sumTheta = 0;

    std::map<int, int> system;

    for (const auto& c : cluster) {

      auto it = m_cellCache.find(c.cellID);
      if (it == m_cellCache.end())
        continue;

      const auto& hit = it->second;

      // identify calo system
      auto sys = m_decoder->get(hit.getCellID(), m_indexSystem);
      system[(int)sys]++;

      // get position, calculate phi/theta
      cx += c.x * c.energy;
      cy += c.y * c.energy;
      cz += c.z * c.energy;
      dd4hep::Position pos(c.x, c.y, c.z);

      phiVec.push_back(pos.Phi());
      thetaVec.push_back(pos.Theta());
      eVec.push_back(c.energy);

      sumPhi += pos.Phi() * c.energy;
      sumTheta += pos.Theta() * c.energy;

      if (m_createClusterCellCollection) {
        auto edmCell = hit.clone();
        outClusterCells->push_back(edmCell);
        outCluster.addToHits(edmCell);
      } else {
        outCluster.addToHits(hit);
      }
    }

    outCluster.setPosition(edm4hep::Vector3f(cx / clusterEnergy, cy / clusterEnergy, cz / clusterEnergy));

    sumPhi /= clusterEnergy;
    sumTheta /= clusterEnergy;

    double deltaR = 0.;
    for (size_t i = 0; i < eVec.size(); i++) {
      deltaR += std::sqrt(std::pow(thetaVec[i] - sumTheta, 2) + std::pow(phiVec[i] - sumPhi, 2)) * eVec[i];
    }
    outCluster.addToShapeParameters(deltaR / clusterEnergy);

    outClusters->push_back(outCluster);

    if (system.size() > 1)
      clusterWithMixedCells++;

    checkTotEnergyAboveThreshold += clusterEnergy;
  }

  debug() << "Number of clusters with cells in E and HCal:        " << clusterWithMixedCells << endmsg;
  debug() << "Total energy of clusters:                           " << checkTotEnergy << endmsg;
  debug() << "Total energy of clusters above threshold:           " << checkTotEnergyAboveThreshold << endmsg;
  if (m_createClusterCellCollection) {
    debug() << "Leftover cells :                                    " << allCells.size() - outClusterCells->size()
            << endmsg;
  }

  return StatusCode::SUCCESS;
}

StatusCode CaloTopoClusterFCCee::buildClusters(const std::vector<FastCell>& seedCells,
                                               const std::vector<FastCell>& allCells, FastClusterMap& clusters) const {

  // will create clusters finding in allCells the neighbours of seedCells and
  // of their own neighbours
  // clusters will be saved in a map (clusters) of clusterId -> FastClusters

  verbose() << "Initial number of seeds to loop over: " << seedCells.size() << endmsg;

  // build fast lookup (CellID -> FastCell)
  std::unordered_map<uint64_t, const FastCell*> cellMap;
  cellMap.reserve(allCells.size());
  for (const auto& c : allCells) {
    cellMap.emplace(c.cellID, &c);
  }

  // initialise map of already used cells
  std::unordered_map<uint64_t, uint32_t> used;
  used.reserve(allCells.size());

  // initialise map of clusterId -> set(cellIDs)
  ClusterMaskMap clusterMembers;
  clusterMembers.reserve(seedCells.size());

  // loop over every seeds in calo to build a cluster (or merge with another cluster if appropriate)
  uint32_t seedCounter = 0;
  for (const auto& seed : seedCells) {
    ++seedCounter;
    if (msgLevel() <= MSG::VERBOSE) {
      verbose() << "Looking at seed: " << seedCounter << endmsg;
    }

    // skip already used seeds
    if (used.count(seed.cellID)) {
      if (msgLevel() <= MSG::VERBOSE) {
        verbose() << "Seed is already assigned to another cluster!" << endmsg;
      }
      continue;
    }

    uint32_t clusterId = seedCounter;

    // create cluster
    clusterMembers[clusterId].clear();
    clusters[clusterId].clear();
    clusterMembers[clusterId].reserve(256);
    clusters[clusterId].reserve(256);

    // insert seed and set its type to 1
    clusters[clusterId].push_back(seed);
    clusters[clusterId].back().type = 1;
    used[seed.cellID] = clusterId;
    clusterMembers[clusterId].insert(seed.cellID);

    // attach neighbours
    // recursively find neighbours (nextLayer) of cells in current layer
    std::vector<uint64_t> currentLayer;
    std::vector<uint64_t> nextLayer;
    currentLayer.reserve(1024);
    nextLayer.reserve(1024);

    // start with neighbours of seed
    currentLayer.push_back(seed.cellID);

    // the flag restartBFS is needed to restart the neighbour search
    // after two clusters are merged and currentLayer is updated
    bool restartBFS = false;
    while (!currentLayer.empty()) {

      nextLayer.clear();
      restartBFS = false;

      // loop over cells in currentLayer
      for (uint64_t cid : currentLayer) {

        if (msgLevel() <= MSG::VERBOSE) {
          verbose() << "For cluster: " << clusterId << endmsg;
        }

        auto itCell = cellMap.find(cid);
        if (itCell == cellMap.end())
          continue;

        // retrieve neighbours of cell
        std::vector<uint64_t> neighboursVec;
        if (m_useNeighborMap) {
          neighboursVec = m_neighboursTool->neighbours(cid);
        } else {
          std::set<dd4hep::DDSegmentation::CellID> tmp;
          m_segmentation->neighbours(cid, tmp);
          neighboursVec.assign(tmp.begin(), tmp.end());
        }
        if (neighboursVec.size() == 0) {
          error() << "No neighbours for cellID found! " << endmsg;
          error() << "to cellID :  " << cid << endmsg;
          error() << "in system:   " << m_decoder->get(cid, m_indexSystem) << endmsg;
        }

        // loop over neighbours
        for (uint64_t nid : neighboursVec) {

          auto itN = cellMap.find(nid);
          if (itN == cellMap.end())
            continue;

          // check if neighbour is used
          const FastCell& ncell = *(itN->second);
          auto itUsed = used.find(nid);

          // new cell (not used yet)
          if (itUsed == used.end()) {
            if (msgLevel() <= MSG::VERBOSE) {
              verbose() << "Found neighbour with CellID: " << nid << endmsg;
            }

            bool accept = (ncell.SoverN > m_neighbourSigma) || (m_neighbourSigma == 0);
            if (!accept)
              continue;

            used[nid] = clusterId;
            clusters[clusterId].push_back(FastCell{nid, ncell.energy, ncell.x, ncell.y, ncell.z, 2, ncell.SoverN});
            clusterMembers[clusterId].insert(nid);
            nextLayer.push_back(nid);
          } else if (itUsed->second != clusterId) {
            // cell already used for a different cluster -> merge
            const uint32_t source = clusterId;
            const uint32_t target = itUsed->second;
            auto& srcCluster = clusters[source];
            auto& dstCluster = clusters[target];
            auto& dstMask = clusterMembers[target];
            if (msgLevel() <= MSG::VERBOSE) {
              verbose() << "This neighbour was found in cluster " << target << ", cluster " << source
                        << " will be merged!" << endmsg;
              verbose() << "Assigning all cells ( " << srcCluster.size() << " ) to Cluster " << target << " with ( "
                        << dstCluster.size() << " ). " << endmsg;
            }

            // move ALL cells from source -> target
            for (const auto& c : srcCluster) {

              // update global ownership FIRST
              used[c.cellID] = target;

              // avoid duplicates
              if (!dstMask.insert(c.cellID).second)
                continue;

              dstCluster.push_back(c);
            }

            // erase source cluster completely
            srcCluster.clear();
            clusters.erase(source);
            clusterMembers.erase(source);

            // update current clusterId
            clusterId = target;

            // continue search from merged context
            restartBFS = true;
            currentLayer.clear();
            for (const auto& c : clusters[clusterId]) {
              currentLayer.push_back(c.cellID);
            }
            break; // leave neighbour loop
          }
        } // end of neighbour loop
        if (msgLevel() <= MSG::VERBOSE) {
          verbose() << "Found " << nextLayer.size() << " neighbours.." << endmsg;
        }
        if (restartBFS)
          break; // leave currentLayer loop
      } // end of currentLayer loop

      if (restartBFS)
        continue; // restart search for currentLayer since we rebuilt it after merge

      // go to next layer
      currentLayer.swap(nextLayer);
    }

    // add last set of neighbours
    // loop over the cells of the cluster
    int lastNeighbours(0);
    for (const auto& c : clusters[clusterId]) {
      int lastNeighboursOfCell(0);

      // skip type-3 cells (last neighbours)
      if (c.type > 2)
        continue;

      // retrieve neighbours of cell
      std::vector<uint64_t> lastVec;
      if (m_useNeighborMap) {
        lastVec = m_neighboursTool->neighbours(c.cellID);
      } else {
        std::set<dd4hep::DDSegmentation::CellID> tmp;
        m_segmentation->neighbours(c.cellID, tmp);
        lastVec.assign(tmp.begin(), tmp.end());
      }

      // loop over neighbours
      // verbose() << "Number of neighbours for cell " << c.cellID << " : " << lastVec.size() << endmsg;
      for (uint64_t nid : lastVec) {
        // skip already used neighbours
        if (used.count(nid))
          continue;

        // find corresponding cell
        auto itN = cellMap.find(nid);
        if (itN == cellMap.end())
          continue;
        const FastCell& ncell = *(itN->second);

        // skip cells with S/N below threshold
        if (m_lastNeighbourSigma != 0 && ncell.SoverN < m_lastNeighbourSigma)
          continue;

        // add cell to cluster
        used[nid] = clusterId;
        clusters[clusterId].push_back(FastCell{nid, ncell.energy, ncell.x, ncell.y, ncell.z, 3, ncell.SoverN});
        clusterMembers[clusterId].insert(nid);
        lastNeighboursOfCell++;
        lastNeighbours++;
      }
      // verbose() << "Added " << lastNeighboursOfCell << " last neighbours for cell " << c.cellID << endmsg;
    }
    if (msgLevel() <= MSG::VERBOSE) {
      verbose() << "Found " << lastNeighbours << " last neighbours.." << endmsg;
    }
  }

  return StatusCode::SUCCESS;
}

StatusCode CaloTopoClusterFCCee::finalize() {
  delete m_decoder;
  for (size_t ih = 0; ih < m_cellCollectionHandles.size(); ih++)
    delete m_cellCollectionHandles[ih];

  return Gaudi::Algorithm::finalize();
}
