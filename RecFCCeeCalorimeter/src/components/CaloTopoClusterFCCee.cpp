#include "CaloTopoClusterFCCee.h"
#include "../../../RecCalorimeter/src/components/NoiseCaloCellsFromFileTool.h"

// std
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <edm4hep/CalorimeterHitCollectionData.h>
#include <edm4hep/MutableCalorimeterHit.h>
#include <edm4hep/MutableCluster.h>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Gaudi
#include <GaudiKernel/MsgStream.h>

// k4FWCore
#include "k4Interface/IGeoSvc.h"

// Datamodel
#include "edm4hep/Cluster.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DD4hep/Readout.h"


DECLARE_COMPONENT(CaloTopoClusterFCCee)


CaloTopoClusterFCCee::CaloTopoClusterFCCee(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
  declareProperty("cells", m_inCells,
                  "Handle for cells from ECal Barrel (input collection)");
  declareProperty("noiseTool", m_noiseTool, "Handle for the cells noise tool");
  declareProperty("neigboursTool", m_neighboursTool, "Handle for tool to retrieve cell neighbors");
  declareProperty("clusters", m_clusterCollection,
                  "Handle for calo clusters (output collection)");
  declareProperty("clusterCells", m_clusterCellsCollection,
                  "Handle for clusters (output collection)");
}
StatusCode CaloTopoClusterFCCee::initialize() {
  if (GaudiAlgorithm::initialize().isFailure()) return StatusCode::FAILURE;
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service. "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }
  if (!m_neighboursTool.retrieve()) {
    error() << "Unable to retrieve the cells neighbours tool!!!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (!m_noiseTool.retrieve()) {
    error() << "Unable to retrieve the cells noise tool!!!" << endmsg;
    return StatusCode::FAILURE;
  }
  m_decoder_ecal = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder();
  m_index_layer_ecal = m_decoder_ecal->index("layer");

  return StatusCode::SUCCESS;
}

StatusCode CaloTopoClusterFCCee::execute() {
  // Get the input collection with calorimeter cells
  const edm4hep::CalorimeterHitCollection* inCells = m_inCells.get();
  debug() << "Input calo cell collection size: " << inCells->size() << endmsg;

  // On first event, create cache
  if (m_min_noise.empty()) {
    createCache(inCells);
  }

  // Create output collections
  edm4hep::ClusterCollection* outClusters = m_clusterCollection.createAndPut();
  edm4hep::CalorimeterHitCollection* outClusterCells = m_clusterCellsCollection.createAndPut();

  // Find seeds
  edm4hep::CalorimeterHitCollection seedCells = findSeeds(inCells);
  debug() << "Number of seeds found :    " << seedCells.size() << endmsg;

  // Build protoclusters (find neighbouring cells)
  std::map<uint32_t, edm4hep::CalorimeterHitCollection> protoClusters;
  {
    StatusCode sc = buildProtoClusters(seedCells,
                                       inCells,
                                       protoClusters);
    if (sc.isFailure()) {
      error() << "Unable to build the protoclusters!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // Build clusters
  debug() << "Building " << protoClusters.size() << " clusters." << endmsg;
  double checkTotEnergy = 0.;
  int clusterWithMixedCells = 0;
  for (const auto& protoCluster: protoClusters) {
    edm4hep::MutableCluster cluster;
    //auto& clusterCore = cluster.core();
    double clusterPosX = 0.;
    double clusterPosY = 0.;
    double clusterPosZ = 0.;
    double clusterEnergy = 0.;
    double deltaR = 0.;
    std::vector<double> cellPosPhi(protoCluster.second.size(), 0);
    std::vector<double> cellPosTheta(protoCluster.second.size(), 0);
    std::vector<double> cellEnergy(protoCluster.second.size(), 0);
    double sumCellPhi = 0.;
    double sumCellTheta = 0.;
    std::map<int, int> system;

    debug() << "Building cluster with ID: " << protoCluster.first << endmsg;

    for (const auto& cell: protoCluster.second) {
      clusterEnergy += cell.getEnergy();
      // identify calo system
      auto systemId = m_decoder->get(cell.getCellID(), m_index_system);
      system[int(systemId)]++;
      auto cellPos = dd4hep::Position(cell.getPosition().x,
                                      cell.getPosition().y,
                                      cell.getPosition().z);

      clusterPosX += cell.getPosition().x * cell.getEnergy();
      clusterPosY += cell.getPosition().y * cell.getEnergy();
      clusterPosZ += cell.getPosition().z * cell.getEnergy();
      cellPosPhi.push_back(cellPos.Phi());
      cellPosTheta.push_back(cellPos.Theta());
      cellEnergy.push_back(cell.getEnergy());
      sumCellPhi += cellPos.Phi() * cell.getEnergy();
      sumCellTheta += cellPos.Theta() * cell.getEnergy();

      cluster.addToHits(cell);
      outClusterCells->push_back(cell);
    }
    cluster.setEnergy(clusterEnergy);
    cluster.setPosition(edm4hep::Vector3f(clusterPosX / clusterEnergy,
                                          clusterPosY / clusterEnergy,
                                          clusterPosZ / clusterEnergy));
    // store deltaR of cluster in time for the moment..
    sumCellPhi = sumCellPhi / clusterEnergy;
    sumCellTheta = sumCellTheta / clusterEnergy;
    for (size_t i = 0; i < cellEnergy.size(); ++i) {
      deltaR += std::sqrt(std::pow(cellPosTheta[i] - sumCellTheta, 2) +
                          std::pow(cellPosPhi[i] - sumCellPhi, 2))
                * cellEnergy[i];
    }
    cluster.addToShapeParameters(deltaR / clusterEnergy);
    verbose() << "Cluster energy:     " << cluster.getEnergy() << endmsg;
    checkTotEnergy += cluster.getEnergy();

    outClusters->push_back(cluster);
    if (system.size() > 1)
      clusterWithMixedCells++;

    cellPosPhi.clear();
    cellPosTheta.clear();
    cellEnergy.clear();
  }

  debug() << "Number of clusters with cells in E and HCal:        "
          << clusterWithMixedCells << endmsg;
  debug() << "Total energy of clusters:                           "
          << checkTotEnergy << endmsg;
  debug() << "Leftover cells :                                    "
          << outClusterCells->size() - inCells->size() << endmsg;

  return StatusCode::SUCCESS;
}


StatusCode CaloTopoClusterFCCee::finalize() { return GaudiAlgorithm::finalize(); }


edm4hep::CalorimeterHitCollection CaloTopoClusterFCCee::findSeeds(
    const edm4hep::CalorimeterHitCollection* allCells) {

  std::vector<edm4hep::CalorimeterHit> seedCellsVec;

  for (const auto& cell : *allCells) {
    // retrieve the noise const and offset assigned to cell
    // first try to use the cache
    int system = m_decoder->get(cell.getCellID(), m_index_system);
    if (system == 4) { // ECal barrel
      int layer = m_decoder_ecal->get(cell.getCellID(), m_index_layer_ecal);

      double min_threshold = m_min_offset[layer] +
                             m_min_noise[layer] * m_seedSigma;

      debug() << "m_min_offset[layer]   = " << m_min_offset[layer] << endmsg;
      debug() << "m_min_noise[layer]   = " << m_min_noise[layer] << endmsg;
      debug() << "aNumSigma   = " << m_seedSigma << endmsg;
      debug() << "min_threshold   = " << min_threshold << endmsg;
      debug() << "abs(cell.second)   = " << abs(cell.getEnergy()) << endmsg;

      if (std::fabs(cell.getEnergy()) < min_threshold) {
        // if we are below the minimum threshold for the full layer, no need to
        // retrieve the exact value
        continue;
      }
    }

    // we are above the minimum threshold of the layer. Let's see if we are
    // above the threshold for this cell.
    double threshold = m_noiseTool->noiseOffset(cell.getCellID()) +
                       (m_noiseTool->noiseRMS(cell.getCellID()) * m_seedSigma);
    if (msgLevel() <= MSG::VERBOSE) {
      debug() << "noise offset    = " << m_noiseTool->noiseOffset(cell.getCellID()) << "GeV " << endmsg;
      debug() << "noise rms       = " << m_noiseTool->noiseRMS(cell.getCellID()) << "GeV " << endmsg;
      debug() << "seed threshold  = " << threshold << "GeV " << endmsg;
      debug() << "======================================" << endmsg;
    }
    if (std::fabs(cell.getEnergy()) > threshold) {
      seedCellsVec.emplace_back(cell);
    }
  }

  // descending order of seeds
  std::sort(seedCellsVec.begin(), seedCellsVec.end(),
            [](const auto& lhs, const auto& rhs) {
              return lhs.getEnergy() < rhs.getEnergy();
            });

  edm4hep::CalorimeterHitCollection seedCells;
  seedCells.setSubsetCollection();
  for (const auto& cell: seedCellsVec) {
    seedCells.push_back(cell);
  }

  return seedCells;
}


StatusCode CaloTopoClusterFCCee::buildProtoClusters(
    const edm4hep::CalorimeterHitCollection& seedCells,
    const edm4hep::CalorimeterHitCollection* allCells,
    std::map<uint32_t, edm4hep::CalorimeterHitCollection>& protoClusters) {

  verbose() << "Initial number of seeds to loop over: " << seedCells.size()
            << endmsg;

  std::map<uint64_t, const edm4hep::CalorimeterHit> allCellsMap;
  for (const auto& cell: *allCells) {
    allCellsMap.emplace(cell.getCellID(), cell);
  }
  std::map<uint64_t, uint32_t> alreadyUsedCells;

  // Loop over every seed in Calo to create first cluster
  uint32_t clusterId = 0;
  for (const auto& seedCell: seedCells) {
    verbose() << "Clusters so far: " << clusterId << endmsg;
    auto seedId = seedCell.getObjectID().index;
    auto cellInCluster = alreadyUsedCells.find(seedId);
    if (cellInCluster != alreadyUsedCells.end()) {
      verbose() << "Seed is already assigned to another cluster!" << endmsg;
      continue;
    }

    clusterId++;
    // new cluster starts with seed
    // set cell type to 1 for seed cell
    edm4hep::MutableCalorimeterHit clusteredCell = seedCell.clone();
    clusteredCell.setType(1);
    protoClusters[clusterId].push_back(clusteredCell);
    alreadyUsedCells[seedId] = clusterId;

    std::vector<std::vector<std::pair<uint64_t, uint32_t>>> nextNeighbours(100);
    nextNeighbours[0] = searchForNeighbours(seedId,
                                            clusterId,
                                            m_neighbourSigma,
                                            allCellsMap,
                                            alreadyUsedCells,
                                            protoClusters,
                                            true);

    // first loop over seeds neighbours
    verbose() << "Found " << nextNeighbours[0].size() << " neighbours.."
              << endmsg;
    int it = 0;
    while (nextNeighbours[it].size() > 0) {
      it++;
      nextNeighbours.emplace_back(std::vector<std::pair<uint64_t, uint>>{});
      for (auto& id : nextNeighbours[it - 1]) {
        if (id.first == 0) {
          error() << "Building of cluster is stopped due to missing cell ID "
                     "in neighbours map!" << endmsg;
          return StatusCode::FAILURE;
        }
        verbose() << "Next neighbours assigned to cluster ID: " << clusterId
                  << endmsg;
        auto additionalNeighbours = searchForNeighbours(id.first,
                                                        clusterId,
                                                        m_neighbourSigma,
                                                        allCellsMap,
                                                        alreadyUsedCells,
                                                        protoClusters,
                                                        true);
        nextNeighbours[it].insert(nextNeighbours[it].end(),
                                  additionalNeighbours.begin(),
                                  additionalNeighbours.end());
      }
      verbose() << "Found " << nextNeighbours[it].size()
                << " more neighbours.." << endmsg;
    }

    // last try with different condition on neighbours
    if (nextNeighbours[it].size() == 0) {
      // loop over all clustered cells
      for (const auto& cell: protoClusters[clusterId]) {
        if (cell.getType() <= 2) {
          verbose() << "Add neighbours of " << cell.getCellID()
                    << " in last round with thr = " << m_lastNeighbourSigma
                    << " x sigma." << endmsg;
          auto lastNeighours = searchForNeighbours(cell.getCellID(),
                                                   clusterId,
                                                   m_lastNeighbourSigma,
                                                   allCellsMap,
                                                   alreadyUsedCells,
                                                   protoClusters,
                                                   false);
        }
      }
    }
  }

  return StatusCode::SUCCESS;
}

std::vector<std::pair<uint64_t, uint32_t>>
CaloTopoClusterFCCee::searchForNeighbours (
    const uint64_t aCellId,
    uint32_t& aClusterID,
    const int aNumSigma,
    std::map<uint64_t, const edm4hep::CalorimeterHit>& allCellsMap,
    std::map<uint64_t, uint32_t>& alreadyUsedCells,
    std::map<uint32_t, edm4hep::CalorimeterHitCollection>& protoClusters,
    const bool aAllowClusterMerge) {
  // Fill vector to be returned, next cell ids and cluster id for which
  // neighbours are found
  std::vector<std::pair<uint64_t, uint32_t>> additionalNeighbours;

  // Retrieve cellIDs of neighbours
  auto neighboursVec = m_neighboursTool->neighbours(aCellId);
  if (neighboursVec.size() == 0) {
    error() << "No neighbours for cellID found! " << endmsg;
    error() << "to cellID :  " << aCellId << endmsg;
    error() << "in system:   " << m_decoder->get(aCellId, m_index_system) << endmsg;
    additionalNeighbours.resize(0);
    additionalNeighbours.push_back(std::make_pair(0, 0));
    return additionalNeighbours;
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
      //
      // first try to use the cache
      bool isBelow = false;
      int system = m_decoder->get(neighbourID, m_index_system);
      if (system == 4) { // ECal barrel
        int layer = m_decoder_ecal->get(neighbourID, m_index_layer_ecal);
        double min_threshold = m_min_offset[layer] +
                               m_min_noise[layer] * aNumSigma;
        if (std::fabs(neighbouringCellEnergy) < min_threshold) {
          // if we are below the minimum threshold for the full layer, no
          // need to retrieve the exact value
          isBelow = true;
        }
      }

      if (isBelow) {
        addNeighbour = false;
      }
      else {
        double thr = m_noiseTool->noiseOffset(neighbourID) +
                     (aNumSigma * m_noiseTool->noiseRMS(neighbourID));
        if (abs(neighbouringCellEnergy) > thr)
          addNeighbour = true;
        else
          addNeighbour = false;
      }
      // give cell type according to threshold
      if (aNumSigma == m_lastNeighbourSigma){
        cellType = 3;
      }
      // if threshold is 0, collect the cell independent on its energy
      if (aNumSigma == 0){
        addNeighbour = true;
      }
      // if neighbour is validated
      if (addNeighbour) {
        // retrieve the cell
        // add neighbour to cells for cluster
        edm4hep::MutableCalorimeterHit clusteredCell =
            allCellsMap[neighbourID].clone();
        clusteredCell.setType(cellType);
        protoClusters[aClusterID].push_back(clusteredCell);
        alreadyUsedCells[neighbourID] = aClusterID;
        additionalNeighbours.push_back(std::make_pair(neighbourID, aClusterID));
      }
    }
    // If cell is hit.. but is assigned to another cluster
    else if (itAllUsedCells != alreadyUsedCells.end() && itAllUsedCells->second != aClusterID && aAllowClusterMerge) {
      uint32_t clusterIDToMerge = itAllUsedCells->second;
      if (msgLevel() <= MSG::VERBOSE){
        verbose() << "This neighbour was found in cluster " << clusterIDToMerge
                  << ", cluster " << aClusterID
                  << " will be merged!" << endmsg;
        verbose() << "Assigning all cells ( "
                  << protoClusters[aClusterID].size() << " ) to Cluster "
                  << clusterIDToMerge << " with ( "
                  << protoClusters[clusterIDToMerge].size()
                  << " ). " << endmsg;
      }
      // Fill all cells into cluster, and assigned cells to new cluster
      alreadyUsedCells[neighbourID] = clusterIDToMerge;
      for (const auto& cell : protoClusters[aClusterID]) {
        alreadyUsedCells[cell.getCellID()] = clusterIDToMerge;
        // make sure that already assigned cells are not added
        for (const auto& cellToMerge : protoClusters[clusterIDToMerge]) {
          if (cellToMerge == cell) continue;
        }
        protoClusters[clusterIDToMerge].push_back(cell);
      }
      protoClusters.erase(aClusterID);
      // changed clusterId -> if more neighbours are found, correct assignment
      verbose() << "Cluster Id changed to " << clusterIDToMerge << endmsg;
      aClusterID = clusterIDToMerge;
      // found neighbour for next search
      additionalNeighbours.push_back(std::make_pair(neighbourID, aClusterID));
      // end loop to ensure correct cluster assignment
      break;
    }
  }

  return additionalNeighbours;
}


void CaloTopoClusterFCCee::createCache(
    const edm4hep::CalorimeterHitCollection* aCells) {
  // cache the minimum offset and noise per layer for faster lookups down the chain.

  std::unordered_map<int, std::vector<double>> offsets;
  std::unordered_map<int, std::vector<double>> noises;
  std::unordered_set<int> layers;

  // fill all noises and offsets values
  for (const auto& cell : *aCells) {
    int system = m_decoder->get(cell.getCellID(), m_index_system);
    if (system == 4) { //ECal barrel
      int layer = m_decoder_ecal->get(cell.getCellID(), m_index_layer_ecal);
      if (layers.find(layer) == layers.end()) {
        offsets[layer] = std::vector<double>{};
        noises[layer] = std::vector<double>{};
        layers.insert(layer);
      }
      offsets[layer].push_back(m_noiseTool->noiseOffset(cell.getCellID()));
      noises[layer].push_back(m_noiseTool->noiseRMS(cell.getCellID()));
    }
  }

  // then compute the minima
  int num_layers = *std::max_element(layers.begin(), layers.end());
  m_min_noise.resize(num_layers+1);
  m_min_offset.resize(num_layers+1);

  for (auto& off : offsets) {
    m_min_offset[off.first] = *std::min_element(off.second.begin(), off.second.end());
  }
  for (auto& n : noises) {
    m_min_noise[n.first] = *std::min_element(n.second.begin(), n.second.end());
  }
}
