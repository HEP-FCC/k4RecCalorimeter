#ifndef RECFCCEECALORIMETER_CALOTOPOCLUSTERFCCEE_H
#define RECFCCEECALORIMETER_CALOTOPOCLUSTERFCCEE_H

// std
#include <cstdint>
#include <sys/types.h>
#include <vector>
#include <map>

// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"

#include "edm4hep/Constants.h"

// Key4HEP
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/MetaDataHandle.h"
#include "k4Interface/ICaloReadCellNoiseMap.h"
#include "k4Interface/ICaloReadNeighboursMap.h"

// k4SimGeant4
class IGeoSvc;

// EDM4HEP
namespace edm4hep {
class CalorimeterHit;
class CalorimeterHitCollection;
class ClusterCollection;
}

// DD4HEP
namespace dd4hep {
namespace DDSegmentation {
class Segmentation;
class BitFieldCoder;
}
}

/** @class CaloTopoClusterFCCee k4RecCalorimeter/RecFCCeeCalorimeter/src/components/CaloTopoClusterFCCee.h
 *
 *  Algorithm building the topological clusters for the energy reconstruction, following ATLAS note
 *  ATL-LARG-PUB-2008-002.
 *  1. Finds the seeds for the cluster, looking for cells exceeding the signal/noise ratio that is given by "seedSigma".
 *  2. The vector of seeds are sorted energy.
 *  3. Adds the neighbouring cells to the cluster in case their signal/noise ratio is larger than the "neighbourSigma".
 *  4. The found and added neighbours function as next seeds and their neighbours are added until no more cells exceed the threshold.
 *  5. In the last step the neighbours that did not exceed the threshold the first time are tested on "lastNeighbourSigma".
 *  In case that a neighbour is found that has already been assigned to another cluster, both clusters are merged and assigned to the "older" clusterID, this is the one originating from a higher seed energy. The iteration over neighburing cellIDs is continued.
 *  @author Coralie Neubueser
 *  @author Giovanni Marchiori, based on code from Juraj Smiesko
 */

class CaloTopoClusterFCCee : public Gaudi::Algorithm {
public:
  CaloTopoClusterFCCee(const std::string& name, ISvcLocator* svcLoc);

  /**
   *
   */
  StatusCode initialize();

  /**  Find cells with a signal to noise ratio > m_seedSigma.
   *   @param[in] allCells, the map of all cells.
   *   @param[out] the collection of seed cells to build proto-clusters.
   */
  edm4hep::CalorimeterHitCollection findSeeds(const edm4hep::CalorimeterHitCollection* allCells) const;

  /** Build proto-clusters from the found seeds.
   * First the function initialises a cluster in the preClusterCollection for the seed cells,
   * then it calls the CaloTopoClusterFCCee::searchForNeighbours function to retrieve the vector of next cellIDs to add and loop over to find neighbours.
   * The iteration of search for neighbours is continued until no more neihgbours are found. Then a last round of adding neighbouring cells to the cluster is run where the parameter lastNeighbourSigma is applied.
   *   @param[in] seedCells, collection of seeding cells.
   *   @param[in] allCells, collection of all cells.
   *   @param[in] protoClusters, map that is filled with clusterID pointing to the associated cells, in a pair of clsuter index and cell collection
   */
  StatusCode buildProtoClusters(
    const edm4hep::CalorimeterHitCollection& seedCells,
    const edm4hep::CalorimeterHitCollection* allCells,
    std::map<uint32_t, edm4hep::CalorimeterHitCollection>& protoClusters) const;
  /** Search for neighbours and add them to preClusterCollection
   * The 
   *   @param[in] aCellId, the cell ID for which to find the neighbours.
   *   @param[in] aClusterID, the current cluster ID.
   *   @param[in] aNumSigma, the signal/noise ratio to be exceeded by the neighbouring cell to be added to cluster.
   *   @param[in] aCellsMap, map of all cells (CellID, cell pointer).
   *   @param[in] aClusterOfCell, map of cellID to clusterID.
   *   @param[in] protoClusters, map that is filled with clusterID pointing to the associated cells, in a pair of cluster index and cell collection.
   *   @param[in] allowClusterMerge, bool to allow for clusters to be merged, set to false in case of last iteration in CaloTopoClusterFCCee::buildingProtoCluster.
   *   return vector of pairs with cellID and energy of found neighbours.
   */
  std::vector<std::pair<uint64_t, uint32_t>>
  searchForNeighbours(
      const uint64_t aCellId,
      uint32_t& aClusterID,
      const int aNumSigma,
      std::map<uint64_t, const edm4hep::CalorimeterHit>& aCellsMap,
      std::map<uint64_t, uint32_t>& aClusterOfCell,
      std::map<uint32_t, edm4hep::CalorimeterHitCollection>& protoClusters,
      const bool aAllowClusterMerge) const;

  StatusCode execute(const EventContext&) const;

  StatusCode finalize();

private:
  /// List of input cell collections
  Gaudi::Property<std::vector<std::string>> m_cellCollections{this, "cells", {}, "Names of CalorimeterHit collections to read"};
  /// the vector of input DataHandles for the input cell collections
  /// std::vector<DataObjectHandleBase*> m_cellCollectionHandles;
  std::vector<DataHandle<edm4hep::CalorimeterHitCollection>*> m_cellCollectionHandles;
  // Cluster collection (output)
  mutable DataHandle<edm4hep::ClusterCollection> m_clusterCollection{"clusters", Gaudi::DataHandle::Writer, this};
  // Cluster cells in collection (output)
  mutable DataHandle<edm4hep::CalorimeterHitCollection> m_clusterCellsCollection{"clusterCells", Gaudi::DataHandle::Writer, this};
  /// Handle for the cluster shape metadata to write
  MetaDataHandle<std::vector<std::string>> m_shapeParametersHandle{
    m_clusterCollection,
    edm4hep::labels::ShapeParameterNames,
    Gaudi::DataHandle::Writer};
  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;
  /// Handle for the cells noise tool
  mutable ToolHandle<ICaloReadCellNoiseMap> m_noiseTool{"TopoCaloNoisyCells", this};
  /// Handle for neighbours tool
  mutable ToolHandle<ICaloReadNeighboursMap> m_neighboursTool{"TopoCaloNeighbours", this};

  /// Seed threshold in sigma
  Gaudi::Property<int> m_seedSigma{this, "seedSigma", 4, "number of sigma in noise threshold"};
  /// Neighbour threshold in sigma
  Gaudi::Property<int> m_neighbourSigma{this, "neighbourSigma", 2, "number of sigma in noise threshold"};
  /// Last neighbour threshold in sigma
  Gaudi::Property<int> m_lastNeighbourSigma{this, "lastNeighbourSigma", 0, "number of sigma in noise threshold"};
  /// Cluster energy threshold
  Gaudi::Property<float> m_minClusterEnergy{this, "minClusterEnergy", 0., "minimum cluster energy"};
  /// Name of the electromagnetic calorimeter readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "ECalBarrelModuleThetaMerged"};

  /// System encoding string
  Gaudi::Property<std::string> m_systemEncoding{
    this, "systemEncoding", "system:4", "System encoding string"};
  /// General decoder to encode the calorimeter sub-system to determine which
  /// positions tool to use
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  int m_indexSystem;

  /// Decoder for ECal layers
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder_ecal;
  int m_index_layer_ecal;

  // minimum noise and offset per barrel ECal layer
  // this serves as a very small cache for fast lookups
  // and avoid looking into the huge map for most of the cells.
  mutable std::vector<double> m_min_offset;
  mutable std::vector<double> m_min_noise;

  // Utility functions
  void createCache(const edm4hep::CalorimeterHitCollection* aCells) const;
  inline bool cellIdInColl(const uint64_t cellId,
    const edm4hep::CalorimeterHitCollection& coll) const;
};
#endif /* RECFCCEECALORIMETER_CALOTOPOCLUSTERFCCEE_H */
