#include "CaloTopoClusterFCCee.h"
#include "../../../RecCalorimeter/src/components/NoiseCaloCellsFromFileTool.h"

// std
#include <algorithm>
#include <memory>
#include <numeric>
#include <unordered_set>
#include <vector>

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"

// k4FWCore
#include "k4Interface/IGeoSvc.h"

// EDM4hep
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DD4hep/Readout.h"

DECLARE_COMPONENT(CaloTopoClusterFCCee)

CaloTopoClusterFCCee::CaloTopoClusterFCCee(const std::string& name, ISvcLocator* svcLoc) : Gaudi::Algorithm(name, svcLoc) {
  declareProperty("TopoClusterInput", m_inputTool, "Handle for input map of cells");
  declareProperty("noiseTool", m_noiseTool, "Handle for the cells noise tool");
  declareProperty("neigboursTool", m_neighboursTool, "Handle for tool to retrieve cell neighbours");
  declareProperty("positionsECalBarrelTool", m_cellPositionsECalBarrelTool,
                  "Handle for tool to retrieve cell positions in ECal Barrel");
  declareProperty("positionsECalEndcapTool", m_cellPositionsECalEndcapTool,
                  "Handle for tool to retrieve cell positions in ECal Endcap");
  declareProperty("positionsHCalBarrelTool", m_cellPositionsHCalBarrelTool=0,
                  "Handle for tool to retrieve cell positions in HCal Barrel");
  declareProperty("positionsHCalBarrelNoSegTool", m_cellPositionsHCalBarrelNoSegTool=0,
                  "Handle for tool to retrieve cell positions in HCal Barrel without DD4hep segmentation");
  declareProperty("positionsHCalExtBarrelTool", m_cellPositionsHCalExtBarrelTool=0,
                  "Handle for tool to retrieve cell positions in HCal ext Barrel");
  // declareProperty("positionsEMECTool", m_cellPositionsEMECTool, "Handle for tool to retrieve cell positions in EMEC");
  // declareProperty("positionsHECTool", m_cellPositionsHECTool, "Handle for tool to retrieve cell positions in HEC");
  // declareProperty("positionsEMFwdTool", m_cellPositionsEMFwdTool, "Handle for tool to retrieve cell positions EM Fwd");
  // declareProperty("positionsHFwdTool", m_cellPositionsHFwdTool, "Handle for tool to retrieve cell positions Had Fwd");
  declareProperty("clusters", m_clusterCollection, "Handle for calo clusters (output collection)");
  declareProperty("clusterCells", m_clusterCellsCollection, "Handle for clusters (output collection)");
}

StatusCode CaloTopoClusterFCCee::initialize() {
  if (Gaudi::Algorithm::initialize().isFailure()) return StatusCode::FAILURE;
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service. "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }
  if (!m_inputTool.retrieve()) {
    error() << "Unable to retrieve the topo cluster input tool!!!" << endmsg;
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
  // Check if cell position ECal Barrel tool available
  if (!m_cellPositionsECalBarrelTool.empty()) {
    if (!m_cellPositionsECalBarrelTool.retrieve()) {
      error() << "Unable to retrieve ECal Barrel cell positions tool!!!" << endmsg;
      return StatusCode::FAILURE;
    }
  }
  // Check if cell position HCal Barrel tool available (only if name is set so that we can also do ECAL-only topoclustering)
  if (!m_cellPositionsHCalBarrelTool.empty()) {
    if (!m_cellPositionsHCalBarrelTool.retrieve()) {
      error() << "Unable to retrieve HCal Barrel cell positions tool!!!" << endmsg;
      return StatusCode::FAILURE;
    }
  }
  // Check if cell position HCal Barrel no seg tool available (only if name is set so that we can also do ECAL-only topoclustering)
  if (!m_cellPositionsHCalBarrelNoSegTool.empty()) {
    if (!m_cellPositionsHCalBarrelNoSegTool.retrieve()) {
      error() << "Unable to retrieve HCal Barrel no segmentation cell positions tool!!!" << endmsg;
      return StatusCode::FAILURE;
    }
  }
  // Check if cell position HCal Endcap tool available (only if name is set so that we can also do ECAL-only topoclustering)
  if (!m_cellPositionsHCalExtBarrelTool.empty()) {
    if (!m_cellPositionsHCalExtBarrelTool.retrieve()) {
      error() << "Unable to retrieve HCal Ext Barrel cell positions tool!!!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  m_decoder_ecal = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder();
  m_index_layer_ecal = m_decoder_ecal->index("layer");

  // initialise the list of metadata for the clusters
  std::vector<std::string> shapeParameterNames = {"dR_over_E"};
  m_shapeParametersHandle.put(shapeParameterNames);

  return StatusCode::SUCCESS;
}

StatusCode CaloTopoClusterFCCee::execute(const EventContext&) const {
  // Create output collections
  auto edmClusters = m_clusterCollection.createAndPut();
  auto edmClusterCells = m_clusterCellsCollection.createAndPut();

  std::unordered_map<uint64_t, double> allCells;
  std::vector<std::pair<uint64_t, double>> firstSeeds;

  // Get input cell map from the input tool
  StatusCode sc_prepareCellMap = m_inputTool->cellIDMap(allCells);
  if (sc_prepareCellMap.isFailure()) {
    error() << "Unable to create cell map!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (allCells.empty()) {
    debug() << "No active cells, skipping event..." << endmsg;
    return StatusCode::SUCCESS;
  }
  debug() << "Number of active cells:                              : " << allCells.size() << endmsg;

  // On first event, create cache
  if (m_min_noise.empty()) {
    createCache(allCells);
  }

  // Find the seeds
  CaloTopoClusterFCCee::findingSeeds(allCells, m_seedSigma, firstSeeds);
  debug() << "Number of seeds found                                : " << firstSeeds.size() << endmsg;

  // Sort the seeds in decending order of their energy
  std::sort(firstSeeds.begin(), firstSeeds.end(),
            [](const std::pair<uint64_t, double>& lhs, const std::pair<uint64_t, double>& rhs) {
              return lhs.second < rhs.second;
            });

  // Build the clusters
  std::map<uint, std::vector<std::pair<uint64_t, int>>> preClusterCollection;
  StatusCode sc_buildProtoClusters = buildingProtoCluster(m_neighbourSigma,
                                                          m_lastNeighbourSigma,
                                                          firstSeeds, allCells,
                                                          preClusterCollection);
  if (sc_buildProtoClusters.isFailure()) {
    error() << "Unable to build the protoclusters!" << endmsg;
    return StatusCode::FAILURE;
  }

  // Calculate the cluster properties and save them in the edm
  debug() << "Building " << preClusterCollection.size() << " clusters" << endmsg;
  double checkTotEnergy = 0.;
  double checkTotEnergyAboveThreshold = 0.;
  int clusterWithMixedCells = 0;
  // loop over the clusters
  for (auto i : preClusterCollection) {
    //auto& clusterCore = cluster.core();
    double posX = 0.;
    double posY = 0.;
    double posZ = 0.;
    double energy = 0.;
    double deltaR = 0.;
    std::vector<double> posPhi (i.second.size());
    std::vector<double> posTheta (i.second.size());
    std::vector<double> vecEnergy (i.second.size());
    double sumPhi = 0.;
    double sumTheta = 0.;
    std::map<int,int> system;

    // calculate cluster energy and decide whether to keep it
    for (auto pair : i.second) {
      dd4hep::DDSegmentation::CellID cID = pair.first;
      energy += allCells[cID];
    }
    verbose() << "Cluster energy:     " << energy << endmsg;
    checkTotEnergy += energy;
  
    if (energy<m_minClusterEnergy) {
      continue;
    }

    edm4hep::MutableCluster cluster;
    cluster.setEnergy(energy);
    checkTotEnergyAboveThreshold += cluster.getEnergy();
    // loop over the cells attached to the cluster to calculate cluster barycenter and attach cells
    for (auto pair : i.second) {
      dd4hep::DDSegmentation::CellID cID = pair.first;
      // get CalorimeterHit by cellID
      auto newCell = edm4hep::MutableCalorimeterHit();
      newCell.setEnergy(allCells[cID]);
      newCell.setCellID(cID);
      newCell.setType(pair.second);

      // identify calo system, retrieve positioning tool and
      // get cell position by cellID
      auto systemId = m_decoder->get(cID, m_index_system);
      system[int(systemId)]++;
      dd4hep::Position posCell;
      if (systemId == 4)  // ECAL BARREL system id
        posCell = m_cellPositionsECalBarrelTool->xyzPosition(cID);
      else if (systemId == 5)
	posCell = m_cellPositionsECalEndcapTool->xyzPosition(cID);
      else if (systemId == 8){  // HCAL BARREL system id
	if (m_noSegmentationHCalUsed)
	  posCell = m_cellPositionsHCalBarrelNoSegTool->xyzPosition(cID);
	else
	  posCell = m_cellPositionsHCalBarrelTool->xyzPosition(cID);
      }
      else if (systemId == 9)  // HCAL ENDCAP system id
        posCell = m_cellPositionsHCalExtBarrelTool->xyzPosition(cID);
      //else if (systemId == 6)  // EMEC system id
      //  posCell = m_cellPositionsEMECTool->xyzPosition(cID);
      //else if (systemId == 7)  // HEC system id
      //  posCell = m_cellPositionsHECTool->xyzPosition(cID);
      //else if (systemId == 10)  // EMFWD system id
      //  posCell = m_cellPositionsEMFwdTool->xyzPosition(cID);
      //else if (systemId == 11)  // HFWD system id
      //  posCell = m_cellPositionsHFwdTool->xyzPosition(cID);
      else
	warning() << "No cell positions tool found for system id " << systemId << ". " << endmsg;
      newCell.setPosition(edm4hep::Vector3f{
	  static_cast<float>(posCell.X() / dd4hep::mm),
	  static_cast<float>(posCell.Y() / dd4hep::mm),
	  static_cast<float>(posCell.Z() / dd4hep::mm)});

      posX += posCell.X() * newCell.getEnergy();
      posY += posCell.Y() * newCell.getEnergy();
      posZ += posCell.Z() * newCell.getEnergy();

      posPhi.push_back(posCell.Phi());
      posTheta.push_back(posCell.Theta());
      vecEnergy.push_back(newCell.getEnergy());
      sumPhi += posCell.Phi() * newCell.getEnergy();
      sumTheta += posCell.Theta() * newCell.getEnergy();

      cluster.addToHits(newCell);
      edmClusterCells->push_back(newCell);
      auto er = allCells.erase(cID);
      if (er!=1)
	info() << "Problem in erasing cell ID from map." << endmsg;
      }
    cluster.setPosition(edm4hep::Vector3f{
	static_cast<float>((posX / energy) / dd4hep::mm),
	static_cast<float>((posY / energy) / dd4hep::mm),
	static_cast<float>((posZ / energy) / dd4hep::mm)});
    // store deltaR of cluster in time for the moment..
    sumPhi = sumPhi / energy;
    sumTheta = sumTheta / energy;
    int counter = 0;
    for (auto entryTheta : posTheta){
      deltaR += sqrt(pow(entryTheta-sumTheta,2) + pow(posPhi[counter]-sumPhi,2)) * vecEnergy[counter];
      counter++;
    }
    cluster.addToShapeParameters(deltaR / energy);
    edmClusters->push_back(cluster);

    if (system.size() > 1)
      clusterWithMixedCells++;

    posPhi.clear();
    posTheta.clear();
    vecEnergy.clear();
  }

  debug() << "Number of clusters                                   : " << preClusterCollection.size() << endmsg;
  debug() << "Total energy of clusters                             : " << checkTotEnergy << endmsg;
  debug() << "Number of clusters above threshold                   : " << edmClusters->size() << endmsg;
  debug() << "Total energy of clusters above threshold             : " << checkTotEnergyAboveThreshold << endmsg;
  debug() << "Clusters above threshold with cells in ECal and HCal : " << clusterWithMixedCells << endmsg;
  debug() << "Leftover cells (from clusters above threshold)       : " << allCells.size() << endmsg;
  return StatusCode::SUCCESS;
}

void CaloTopoClusterFCCee::findingSeeds(const std::unordered_map<uint64_t, double>& aCells,
                                        int aNumSigma,
                                        std::vector<std::pair<uint64_t, double>>& aSeeds) const {
  for (const auto& cell : aCells) {
    // retrieve the noise const and offset assigned to cell
    // first try to use the cache
    int system = m_decoder->get(cell.first, m_index_system);
    if (system == 4) { //ECal barrel
      int layer = m_decoder_ecal->get(cell.first, m_index_layer_ecal);

      double min_threshold = m_min_offset[layer] + m_min_noise[layer] * aNumSigma;

      verbose() << "m_min_offset[layer]   = " << m_min_offset[layer] << endmsg;
      verbose() << "m_min_noise[layer]   = " << m_min_noise[layer] << endmsg;
      verbose() << "aNumSigma   = " << aNumSigma << endmsg;
      verbose() << "min_threshold   = " << min_threshold << endmsg;
      verbose() << "abs(cell.second)   = " << abs(cell.second) << endmsg;

      if (abs(cell.second) < min_threshold) {
        // if we are below the minimum threshold for the full layer, no need to retrieve the exact value
        continue;
      }
    }

    // we are above the minimum threshold of the layer. Let's see if we are above the threshold for this cell.
    double threshold = m_noiseTool->noiseOffset(cell.first) + ( m_noiseTool->noiseRMS(cell.first) * aNumSigma);
    if (msgLevel() <= MSG::VERBOSE){
      debug() << "noise offset    = " << m_noiseTool->noiseOffset(cell.first) << "GeV " << endmsg;
      debug() << "noise rms       = " << m_noiseTool->noiseRMS(cell.first) << "GeV " << endmsg;
      debug() << "seed threshold  = " << threshold << "GeV " << endmsg;
      debug() << "======================================" << endmsg;
    }
    if (abs(cell.second) > threshold) {
      aSeeds.emplace_back(cell.first, cell.second);
    }
  }
}

StatusCode CaloTopoClusterFCCee::buildingProtoCluster(
    int aNumSigma,
    int aLastNumSigma,
    std::vector<std::pair<uint64_t, double>>& aSeeds,
    const std::unordered_map<uint64_t, double>& aCells,
    std::map<uint, std::vector< std::pair<uint64_t, int>>>& aPreClusterCollection) const {
  // Map of cellIDs to clusterIds
  std::map<uint64_t, uint> clusterOfCell;

  // Loop over every seed in Calo to create first cluster
  uint iSeeds = 0;
  verbose() << "seeds to loop over : " << aSeeds.size() << endmsg;
  for (auto& itSeed : aSeeds) {
    iSeeds++;
    verbose() << "Seed num: " << iSeeds << endmsg;
    auto seedId = itSeed.first;
    auto cellInCluster = clusterOfCell.find(seedId);
    if (cellInCluster != clusterOfCell.end()) {
      verbose() << "Seed is already assigned to another cluster!" << endmsg;
      continue;
    } else {
      // new cluster starts with seed
      // set cell Bits to 1 for seed cell
      aPreClusterCollection[iSeeds].push_back(std::make_pair(seedId, 1));
      uint clusterId = iSeeds;
      clusterOfCell[seedId] = clusterId;

      std::vector<std::vector<std::pair<uint64_t, uint>>> vecNextNeighbours(1);
      vecNextNeighbours.reserve(100);
      vecNextNeighbours[0] = CaloTopoClusterFCCee::searchForNeighbours(seedId, clusterId, aNumSigma, aCells, clusterOfCell,
                                                     aPreClusterCollection, true);
      // first loop over seeds neighbours
      verbose() << "Found " << vecNextNeighbours[0].size() << " neighbours.." << endmsg;
      int it = 0;
      while (vecNextNeighbours[it].size() > 0) {
        it++;
        vecNextNeighbours.emplace_back(std::vector<std::pair<uint64_t, uint>>{});
        for (auto& id : vecNextNeighbours[it - 1]) {
          if (id.first == 0){
            error() << "Building of cluster is stopped due to missing id in neighbours map." << endmsg;
            return StatusCode::FAILURE;
          }
          verbose() << "Next neighbours assigned to clusterId : " << clusterId << endmsg;
          auto vec = CaloTopoClusterFCCee::searchForNeighbours(id.first, clusterId, aNumSigma, aCells, clusterOfCell,
                                                               aPreClusterCollection, true);
          vecNextNeighbours[it].insert(vecNextNeighbours[it].end(), vec.begin(), vec.end());
        }
        verbose() << "Found " << vecNextNeighbours[it].size() << " more neighbours.." << endmsg;
      }
      // last try with different condition on neighbours
      if (vecNextNeighbours[it].size() == 0) {
        auto clusteredCells = aPreClusterCollection[clusterId];
        // loop over all clustered cells
        for (auto& id : clusteredCells) {
          if (id.second <= 2){
            verbose() << "Add neighbours of " << id.first << " in last round with thr = " << aLastNumSigma << " x sigma." << endmsg;
            auto lastNeighours = CaloTopoClusterFCCee::searchForNeighbours(id.first, clusterId, aLastNumSigma, aCells, clusterOfCell,
                                                                      aPreClusterCollection, false);
          }
        }
      }
    }
  }
  return StatusCode::SUCCESS;
}

std::vector<std::pair<uint64_t, uint> >
CaloTopoClusterFCCee::searchForNeighbours(const uint64_t aCellId,
                                     uint& aClusterID,
                                     int aNumSigma,
                                     const std::unordered_map<uint64_t, double>& aCells,
                                     std::map<uint64_t, uint>& aClusterOfCell,
                                     std::map<uint, std::vector<std::pair<uint64_t, int>>>& aPreClusterCollection,
                                     bool aAllowClusterMerge) const {
  // EWV debug
  unsigned modSeed = (aCellId >> 28) & 0x7ff;
  unsigned rhoSeed = (aCellId >> 39) & 0xff;
  unsigned zSeed = (aCellId >> 47) & 0xff; 
  // Fill vector to be returned, next cell ids and cluster id for which neighbours are found
  std::vector<std::pair<uint64_t, uint>> addedNeighbourIds;
  // Retrieve cellIDs of neighbours
  auto neighboursVec = m_neighboursTool->neighbours(aCellId);
  if (neighboursVec.size() == 0) {
    error() << "No neighbours for cellID found! " << endmsg;
    error() << "to cellID :  " << aCellId << endmsg;
    error() << "in system:   " << m_decoder->get(aCellId, m_index_system) << endmsg;
    addedNeighbourIds.resize(0);
    addedNeighbourIds.push_back(std::make_pair(0, 0));
  } else {

    verbose() << "For cluster: " << aClusterID << endmsg;
    // loop over neighbours
    for (auto& itr : neighboursVec) {
      auto neighbourID = itr;
      debug() << "neighbor ID is " << neighbourID << endmsg;
      unsigned modNeighbor = (neighbourID >> 28) & 0x7ff;
      unsigned rhoNeighbor = (neighbourID >> 39) & 0xff;
      unsigned zNeighbor = (neighbourID >> 47) & 0xff;
      if (rhoSeed != rhoNeighbor) {
	debug() << "checking neighbor with different rho than seed " << modSeed << " " << rhoSeed << " " << zSeed << " " << modNeighbor << " " << rhoNeighbor << " " << zNeighbor << " " << aCellId << " " << neighbourID << endmsg;
      }

      // Find the neighbour in the Calo cells list
      auto itAllCells = aCells.find(neighbourID);
      auto itAllUsedCells = aClusterOfCell.find(neighbourID);

      if (itAllCells == aCells.end()) {
	debug() << "Can't find this neighbour in the list of cells!" << endmsg;
      }
      // If cell is hit.. and is not assigned to a cluster
      if (itAllCells != aCells.end() && itAllUsedCells == aClusterOfCell.end()) {
        verbose() << "Found neighbour with CellID: " << neighbourID << endmsg;

        auto neighbouringCellEnergy = itAllCells->second;
        bool addNeighbour = false;
        int cellType = 2;
        // retrieve the cell noise level [GeV]
	//
	// first try to use the cache
	bool is_below = false;
	int system = m_decoder->get(neighbourID, m_index_system);
	if (system == 4) { //ECal barrel
	  int layer = m_decoder_ecal->get(neighbourID, m_index_layer_ecal);
	  double min_threshold = m_min_offset[layer] + m_min_noise[layer] * aNumSigma;
	  if (abs(neighbouringCellEnergy) < min_threshold) {
	    // if we are below the minimum threshold for the full layer, no need to retrieve the exact value
	    is_below = true;
	  }
	}
	
	if (is_below) {
	  addNeighbour = false;
	}
	else {
	  double thr = m_noiseTool->noiseOffset(neighbourID) + (aNumSigma * m_noiseTool->noiseRMS(neighbourID));
	  debug() << "Neighboring cell energy is " << neighbouringCellEnergy << endmsg;
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
          aPreClusterCollection[aClusterID].push_back(std::make_pair(neighbourID, cellType));
          aClusterOfCell[neighbourID] = aClusterID;
          addedNeighbourIds.push_back(std::make_pair(neighbourID, aClusterID));
        }
      }
      // If cell is hit.. but is assigned to another cluster
      else if (itAllUsedCells != aClusterOfCell.end() && itAllUsedCells->second != aClusterID && aAllowClusterMerge) {
        uint clusterIDToMerge = itAllUsedCells->second;
        if (msgLevel() <= MSG::VERBOSE){
          verbose() << "This neighbour was found in cluster " << clusterIDToMerge << ", cluster " << aClusterID
                    << " will be merged!" << endmsg;
          verbose() << "Assigning all cells ( " << aPreClusterCollection.find(aClusterID)->second.size() << " ) to Cluster "
                    << clusterIDToMerge << " with ( " << aPreClusterCollection.find(clusterIDToMerge)->second.size()
                    << " ). " << endmsg;
        }
        // Fill all cells into cluster, and assigned cells to new cluster
        aClusterOfCell[neighbourID] = clusterIDToMerge;
        for (auto& i : aPreClusterCollection.find(aClusterID)->second) {
          aClusterOfCell[i.first] = clusterIDToMerge;
          bool found = false;
          // make sure that already assigned cells are not added
          for (auto& j : aPreClusterCollection[clusterIDToMerge]) {
            if (j.first == i.first) found = true;
          }
          if (!found) {
            aPreClusterCollection[clusterIDToMerge].push_back(std::make_pair(i.first, i.second));
          }
        }
        aPreClusterCollection.erase(aClusterID);
        // changed clusterId -> if more neighbours are found, correct assignment
        verbose() << "Cluster Id changed to " << clusterIDToMerge << endmsg;
        aClusterID = clusterIDToMerge;
        // found neighbour for next search
        addedNeighbourIds.push_back(std::make_pair(neighbourID, aClusterID));
        // end loop to ensure correct cluster assignment
        break;
      }
    }
  }
  return addedNeighbourIds;
}

StatusCode CaloTopoClusterFCCee::finalize() { return Gaudi::Algorithm::finalize(); }


/**
 * \brief Cache the minimum offset and noise per layer for faster lookups down
 * the chain.
 */
void CaloTopoClusterFCCee::createCache(const std::unordered_map<uint64_t, double>& aCells) const {
  std::unordered_map<int, std::vector<double>> offsets;
  std::unordered_map<int, std::vector<double>> noises;
  std::unordered_set<int> layers;

  // Fill all noises and offsets values
  for (const auto& cell : aCells) {
    int system = m_decoder->get(cell.first, m_index_system);
    if (system == 4 || system == 5) { //ECal barrel or endcap
      int layer = m_decoder_ecal->get(cell.first, m_index_layer_ecal);
      if (layers.find(layer) == layers.end()) {
        offsets[layer] = std::vector<double>{};
        noises[layer] = std::vector<double>{};
        layers.insert(layer);
      }
      offsets[layer].push_back(m_noiseTool->noiseOffset(cell.first));
      noises[layer].push_back(m_noiseTool->noiseRMS(cell.first));
    }
  }

  // if ECal barrel cells are not included in the input then return, otherwise it will crash
  if(layers.empty()) return;

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
