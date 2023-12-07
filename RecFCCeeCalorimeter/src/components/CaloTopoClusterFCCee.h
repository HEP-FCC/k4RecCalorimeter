#ifndef RECFCCEECALORIMETER_CALOTOPOCLUSTERFCCEE_H
#define RECFCCEECALORIMETER_CALOTOPOCLUSTERFCCEE_H

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"

// FCCSW
#include "k4FWCore/DataHandle.h"
#include "k4Interface/ICaloReadCellNoiseMap.h"
#include "k4Interface/ICaloReadNeighboursMap.h"
// #include "k4Interface/ICalorimeterTool.h"
#include "k4Interface/ICellPositionsTool.h"
#include "k4Interface/ITopoClusterInputTool.h"

#include <Gaudi/Property.h>
#include <cstdint>
#include <sys/types.h>
#include <vector>
#include <map>

class IGeoSvc;

// datamodel
namespace edm4hep {
class CalorimeterHit;
class CalorimeterHitCollection;
class ClusterCollection;
}

namespace DD4hep {
namespace DDSegmentation {
class Segmentation;
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
 */

class CaloTopoClusterFCCee : public GaudiAlgorithm {
public:
  CaloTopoClusterFCCee(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize();

  /**  Find cells with a signal to noise ratio > seedSigma.
   *   @param[in] allCells, collection of all cells.
   *   return collection of seed cells.
   */
  edm4hep::CalorimeterHitCollection findSeeds(
      const edm4hep::CalorimeterHitCollection* allCells);

  /** Building proto-clusters from the found seeds.
   * First the function initialises a cluster in the preClusterCollection for the seed cells,
   * then it calls the CaloTopoClusterFCCee::searchForNeighbours function to retrieve the vector of next cellIDs to add and loop over to find neighbours.
   * The iteration of search for neighbours is continued until no more neihgbours are found. Then a last round of adding neighbouring cells to the cluster is run where the parameter lastNeighbourSigma is applied.
   *   @param[in] aSeeds, vector of seeding cells.
   *   @param[in] aCells, map of all cells.
   *   @param[in] aPreClusterCollection, map that is filled with clusterID pointing to the associated cells, in a pair of cellID and cellType.
   */
  StatusCode buildProtoClusters(
      const edm4hep::CalorimeterHitCollection& seedCells,
      const edm4hep::CalorimeterHitCollection* allCells,
      std::map<uint32_t, edm4hep::CalorimeterHitCollection>& protoClusters);

  /** Search for neighbours and add them to preClusterCollection
   * The 
   *   @param[in] aCellId, the cell ID for which to find the neighbours.
   *   @param[in] aClusterID, the current cluster ID.
   *   @param[in] aNumSigma, the signal/noise ratio to be exceeded by the neighbouring cell to be added to cluster.
   *   @param[in] aCells, map of all cells.
   *   @param[in] aClusterOfCell, map of cellID to clusterID.
   *   @param[in] aPreClusterCollection, map that is filled with clusterID pointing to the associated cells, in a pair of cellID and cellType.
   *   @param[in] aAllowClusterMerge, bool to allow for clusters to be merged, set to false in case of last iteration in CaloTopoClusterFCCee::buildingProtoCluster.
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
      const bool aAllowClusterMerge);

  StatusCode execute();

  StatusCode finalize();

private:
  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;

  /// Handle for electromagnetic barrel cells (input collection)
  DataHandle<edm4hep::CalorimeterHitCollection> m_cells{
      "cells", Gaudi::DataHandle::Reader, this};
  // Cluster collection (output)
  DataHandle<edm4hep::ClusterCollection> m_clusterCollection{
      "clusters", Gaudi::DataHandle::Writer, this};
  // Cluster cells in collection (output)
  DataHandle<edm4hep::CalorimeterHitCollection> m_clusterCellsCollection{
      "clusterCells", Gaudi::DataHandle::Writer, this};

  /// Handle for the cells noise tool
  ToolHandle<ICaloReadCellNoiseMap> m_noiseTool{"TopoCaloNoisyCells", this};
  /// Handle for neighbors tool
  ToolHandle<ICaloReadNeighboursMap> m_neighboursTool{"TopoCaloNeighbours", this};

  /// No segmentation used in HCal
  Gaudi::Property<bool> m_noSegmentationHCalUsed{this, "noSegmentationHCal", true, "HCal Barrel readout without DD4hep theta-module segmentation used."};
  /// Seed threshold in sigma
  Gaudi::Property<int> m_seedSigma{this, "seedSigma", 4, "number of sigma in noise threshold"};
  /// Neighbour threshold in sigma
  Gaudi::Property<int> m_neighbourSigma{this, "neighbourSigma", 2, "number of sigma in noise threshold"};
  /// Last neighbour threshold in sigma
  Gaudi::Property<int> m_lastNeighbourSigma{this, "lastNeighbourSigma", 0, "number of sigma in noise threshold"};
  /// Name of the electromagnetic calorimeter readout
  // Gaudi::Property<std::string> m_readoutName{this, "readoutName", "ECalBarrelModuleThetaMerged"};
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", ""};

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
  // this serves as a very small cache for fast lookups and avoid looking into the huge map for most of the cells.
  std::vector<double> m_min_offset;
  std::vector<double> m_min_noise;

  // Utility functions
  void createCache(const edm4hep::CalorimeterHitCollection* aCells);
  inline bool cellIdInColl(const uint64_t cellId,
                           const edm4hep::CalorimeterHitCollection& coll);

};
#endif /* RECFCCEECALORIMETER_CALOTOPOCLUSTERFCCEE_H */
