#ifndef RECFCCEECALORIMETER_CALOTOPOCLUSTERFCCEE_H
#define RECFCCEECALORIMETER_CALOTOPOCLUSTERFCCEE_H

// std
#include <cstdint>
#include <map>
#include <sys/types.h>
#include <utility>
#include <vector>

// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"

// Key4HEP
#include "RecCaloCommon/ICaloReadNeighboursMap.h"
#include "RecCaloCommon/INoiseConstTool.h"
#include "k4FWCore/DataHandle.h"
#include "k4Interface/IGeoSvc.h"

// EDM4HEP
namespace edm4hep {
class CalorimeterHit;
class CalorimeterHitCollection;
class ClusterCollection;
} // namespace edm4hep

// DD4HEP
namespace dd4hep {
namespace DDSegmentation {
  class Segmentation;
  class BitFieldCoder;
} // namespace DDSegmentation
} // namespace dd4hep

/** @class CaloTopoClusterFCCee k4RecCalorimeter/RecFCCeeCalorimeter/src/components/CaloTopoClusterFCCee.h
 *
 *  Algorithm building the topological clusters for the energy reconstruction, following ATLAS note
 *  ATL-LARG-PUB-2008-002.
 *  1. Finds the seeds for the cluster, looking for cells exceeding the signal/noise ratio that is given by "seedSigma".
 *  2. The vector of seeds are sorted energy.
 *  3. Adds the neighbouring cells to the cluster in case their signal/noise ratio is larger than the "neighbourSigma".
 *  4. The found and added neighbours function as next seeds and their neighbours are added until no more cells exceed
 * the threshold.
 *  5. In the last step the neighbours that did not exceed the threshold the first time are tested on
 * "lastNeighbourSigma". In case that a neighbour is found that has already been assigned to another cluster, both
 * clusters are merged and assigned to the "older" clusterID, this is the one originating from a higher seed energy. The
 * iteration over neighburing cellIDs is continued.
 *
 *  @author Coralie Neubueser
 *  @author Giovanni Marchiori - algorithm rewritten for significant speed-up
 */

/// internal cell representation used for clustering, to avoid cloning EDM objects repeatedly
struct FastCell {
  uint64_t cellID;
  float energy;
  float x;
  float y;
  float z;
  uint8_t type; // 0=unused,1=seed,2=neighbour,3=lastNeighbour
  float SoverN;
};
using FastCluster = std::vector<FastCell>;
using FastClusterMap = std::map<uint32_t, FastCluster>; // TODO make it unordered or a vector

class CaloTopoClusterFCCee : public Gaudi::Algorithm {
public:
  CaloTopoClusterFCCee(const std::string& name, ISvcLocator* svcLoc);

  /**
   *
   */
  StatusCode initialize();

  /** Build clusters from the found seeds.
   * First the function initialises a cluster in the preClusterCollection for the seed cells,
   * then it calls the CaloTopoClusterFCCee::searchForNeighbours function to retrieve the vector of next cellIDs to add
   * and loop over to find neighbours. The iteration of search for neighbours is continued until no more neihgbours are
   * found. Then a last round of adding neighbouring cells to the cluster is run where the parameter lastNeighbourSigma
   * is applied.
   *   @param[in] seedCells, collection of seeding cells (vector of fastcells)
   *   @param[in] allCells, collection of all cells (map cellID -> fastcell)
   *   @param[in] clusters, collection of clusters to be filled by the algorithm (map of clusterID -> FastCluster)
   */
  StatusCode buildClusters(const std::vector<FastCell>& seedCells,
                           const std::unordered_map<uint64_t, FastCell>& allCells, FastClusterMap& clusters) const;
  /** Search for neighbours and add them to cluster collection
   *   @param[in] cellID, the cell ID for which to find the neighbours
   *   @param[in] clusterID, the current cluster ID
   *   @param[in] nSigma, the signal/noise ratio to be exceeded by the neighbouring cell to be added to cluster
   *   @param[in] allCells, map of all cells (CellID -> FastCell)
   *   @param[in] usedCells, map of used cells (CellID -> cluster ID)
   *   @param[in] clusters, map that is filled with clusterID pointing to the associated cells, in a pair of
   *              cluster index and cell collection
   *   @param[in] clusterMembers, map (cluster ID -> set of CellIDs) that is filled by the algorithm to keep track of
   * clustered cells
   *   @param[in] allowClusterMerge, bool to allow for clusters to be merged
   *   return vector of cellID of found neighbours
   */
  std::vector<uint64_t> searchForNeighbours(const uint64_t cellID, uint& clusterID, int nSigma,
                                            const std::unordered_map<uint64_t, FastCell>& allCells,
                                            std::unordered_map<uint64_t, uint32_t>& usedCells, FastClusterMap& clusters,
                                            std::unordered_map<uint32_t, std::unordered_set<uint64_t>>& clusterMembers,
                                            bool allowClusterMerge) const;

  StatusCode execute(const EventContext&) const;

  StatusCode finalize();

private:
  /// List of input cell collections
  Gaudi::Property<std::vector<std::string>> m_cellCollections{
      this, "cells", {}, "Names of CalorimeterHit collections to read"};
  /// the vector of input k4FWCore::DataHandles for the input cell collections
  std::vector<k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection>*> m_cellCollectionHandles;
  // Cluster collection (output)
  mutable k4FWCore::DataHandle<edm4hep::ClusterCollection> m_clusterCollection{"clusters", Gaudi::DataHandle::Writer,
                                                                               this};
  // Cluster cells in collection (output)
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_clusterCellsCollection{
      "clusterCells", Gaudi::DataHandle::Writer, this};

  Gaudi::Property<std::vector<int>> m_caloIDs{this, "calorimeterIDs", {}, "Corresponding list of calorimeter IDs"};

  /// Handle for the cells noise tool
  mutable ToolHandle<k4::recCalo::INoiseConstTool> m_noiseTool{"TopoCaloNoisyCells", this};
  /// Handle for neighbours tool
  mutable ToolHandle<k4::recCalo::ICaloReadNeighboursMap> m_neighboursTool{"TopoCaloNeighbours", this};
  // flag to use a pre-calculated neighbor map
  Gaudi::Property<bool> m_useNeighborMap{this, "useNeighborMap", true, "use pre-calculated neighbor map"};
  // use GeoSvc when the neighbor map is not present
  SmartIF<IGeoSvc> m_geoSvc;
  // name of the readout: only needed if useNeighborMap is set to false
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "",
                                             "name of the readout (needed if useNeighborMap=false)"};
  // pointer to the segmentation object
  dd4hep::DDSegmentation::Segmentation* m_segmentation = nullptr;

  /// Seed threshold in sigma
  Gaudi::Property<int> m_seedSigma{this, "seedSigma", 4, "number of sigma in noise threshold"};
  /// Neighbour threshold in sigma
  Gaudi::Property<int> m_neighbourSigma{this, "neighbourSigma", 2, "number of sigma in noise threshold"};
  /// Last neighbour threshold in sigma
  Gaudi::Property<int> m_lastNeighbourSigma{this, "lastNeighbourSigma", 0, "number of sigma in noise threshold"};
  /// Cluster energy threshold
  Gaudi::Property<float> m_minClusterEnergy{this, "minClusterEnergy", 0., "minimum cluster energy"};

  /// System encoding string
  Gaudi::Property<std::string> m_systemEncoding{this, "systemEncoding", "system:4", "System encoding string"};

  /// Flag if a new output cell collection of clustered cells should be created
  Gaudi::Property<bool> m_createClusterCellCollection{this, "createClusterCellCollection", false};
  /// General decoder to encode the calorimeter sub-system to determine which
  /// positions tool to use
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  int m_indexSystem;
};
#endif /* RECFCCEECALORIMETER_CALOTOPOCLUSTERFCCEE_H */
