#ifndef RECCALORIMETER_CALOTOPOCLUSTER_H
#define RECCALORIMETER_CALOTOPOCLUSTER_H

// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"

#include "edm4hep/Constants.h"

// Key4HEP
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/MetaDataHandle.h"
#include "k4Interface/ICaloReadNeighboursMap.h"
#include "k4Interface/ICalorimeterTool.h"
#include "k4Interface/ICellPositionsTool.h"
#include "k4Interface/INoiseConstTool.h"
#include "k4Interface/ITopoClusterInputTool.h"

class IGeoSvc;

// EDM4HEP
namespace edm4hep {
class CalorimeterHit;
class CalorimeterHitCollection;
class ClusterCollection;
} // namespace edm4hep

// DD4HEP
namespace DD4hep {
namespace DDSegmentation {
  class Segmentation;
}
} // namespace DD4hep

/** @class CaloTopoClusterAlgorithm Reconstruction/RecCalorimeter/src/components/CaloTopoCluster.h
 * CombinedCaloTopoCluster.h
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
 *  @author Coralie Neubueser
 */

class CaloTopoCluster : public Gaudi::Algorithm {
public:
  CaloTopoCluster(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize();

  /**  Find cells with a signal to noise ratio > nSigma.
   *   @param[in] aCells, parse the map of all cells.
   *   @param[in] aNumSigma, the signal to noise ratio that the cell has to exceed to become seed.
   *   @param[in] aSeeds, the vector of seed cell ids anf their energy to build proto-clusters.
   */
  virtual void findingSeeds(const std::unordered_map<uint64_t, double>& aCells, int aNumSigma,
                            std::vector<std::pair<uint64_t, double>>& aSeeds) const;

  /** Building proto-clusters from the found seeds.
   * First the function initialises a cluster in the preClusterCollection for the seed cells,
   * then it calls the CaloTopoCluster::searchForNeighbours function to retrieve the vector of next cellIDs to add and
   * loop over to find neighbours. The iteration of search for neighbours is continued until no more neihgbours are
   * found. Then a last round of adding neighbouring cells to the cluster is run where the parameter lastNeighbourSigma
   * is applied.
   *   @param[in] aNumSigma, signal to noise ratio the neighbouring cell has to pass to be added to cluster.
   *   @param[in] aSeeds, vector of seeding cells.
   *   @param[in] aCells, map of all cells.
   *   @param[in] aPreClusterCollection, map that is filled with clusterID pointing to the associated cells, in a pair
   * of cellID and cellType.
   */
  StatusCode buildingProtoCluster(int aNumSigma, int aLastNumSigma, std::vector<std::pair<uint64_t, double>>& aSeeds,
                                  const std::unordered_map<uint64_t, double>& aCells,
                                  std::map<uint, std::vector<std::pair<uint64_t, int>>>& aPreClusterCollection) const;

  /** Search for neighbours and add them to preClusterCollection
   * The
   *   @param[in] aCellId, the cell ID for which to find the neighbours.
   *   @param[in] aClusterID, the current cluster ID.
   *   @param[in] aNumSigma, the signal/noise ratio to be exceeded by the neighbouring cell to be added to cluster.
   *   @param[in] aCells, map of all cells.
   *   @param[in] aClusterOfCell, map of cellID to clusterID.
   *   @param[in] aPreClusterCollection, map that is filled with clusterID pointing to the associated cells, in a pair
   * of cellID and cellType.
   *   @param[in] aAllowClusterMerge, bool to allow for clusters to be merged, set to false in case of last iteration in
   * CaloTopoCluster::buildingProtoCluster. return vector of pairs with cellID and energy of found neighbours.
   */
  std::vector<std::pair<uint64_t, uint>>
  searchForNeighbours(const uint64_t aCellId, uint& aClusterID, int aNumSigma,
                      const std::unordered_map<uint64_t, double>& aCells, std::map<uint64_t, uint>& aClusterOfCell,
                      std::map<uint, std::vector<std::pair<uint64_t, int>>>& aPreClusterCollection,
                      bool aAllowClusterMerge) const;

  StatusCode execute(const EventContext&) const;

  StatusCode finalize();

private:
  // Cluster collection
  mutable k4FWCore::DataHandle<edm4hep::ClusterCollection> m_clusterCollection{"calo/clusters",
                                                                               Gaudi::DataHandle::Writer, this};
  // Cluster cells in collection
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_clusterCellsCollection{
      "calo/clusterCells", Gaudi::DataHandle::Writer, this};
  /// Handle for the cluster shape metadata to write
  k4FWCore::MetaDataHandle<std::vector<std::string>> m_shapeParametersHandle{
      m_clusterCollection, edm4hep::labels::ShapeParameterNames, Gaudi::DataHandle::Writer};
  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;
  /// Handle for the input tool
  mutable ToolHandle<ITopoClusterInputTool> m_inputTool{"TopoClusterInput", this};
  /// Handle for the cells noise tool
  mutable ToolHandle<INoiseConstTool> m_noiseTool{"TopoCaloNoisyCells", this};
  /// Handle for neighbours tool
  mutable ToolHandle<ICaloReadNeighboursMap> m_neighboursTool{"TopoCaloNeighbours", this};
  /// Handle for tool to get positions in ECal Barrel
  mutable ToolHandle<ICellPositionsTool> m_cellPositionsECalBarrelTool{"CellPositionsECalBarrelTool", this};
  /// Handle for tool to get positions in HCal Barrel
  mutable ToolHandle<ICellPositionsTool> m_cellPositionsHCalBarrelNoSegTool{"CellPositionsHCalBarrelNoSegTool", this};
  /// Handle for tool to get positions in HCal Barrel
  mutable ToolHandle<ICellPositionsTool> m_cellPositionsHCalBarrelTool{"CellPositionsHCalBarrelTool", this};
  /// Handle for tool to get positions in HCal Barrel and Ext Barrel, no Segmentation
  mutable ToolHandle<ICellPositionsTool> m_cellPositionsHCalExtBarrelTool{"CellPositionsHCalBarrelNoSegTool", this};
  /// Handle for tool to get positions in Calo Discs
  mutable ToolHandle<ICellPositionsTool> m_cellPositionsEMECTool{"CellPositionsCaloDiscsTool", this};
  /// Handle for tool to get positions in Calo Discs
  mutable ToolHandle<ICellPositionsTool> m_cellPositionsHECTool{"CellPositionsCaloDiscsTool", this};
  /// Handle for tool to get positions in Calo Discs
  mutable ToolHandle<ICellPositionsTool> m_cellPositionsEMFwdTool{"CellPositionsCaloDiscsTool", this};
  /// Handle for tool to get positions in Calo Discs
  mutable ToolHandle<ICellPositionsTool> m_cellPositionsHFwdTool{"CellPositionsCaloDiscsTool", this};

  /// no segmentation used in HCal
  Gaudi::Property<bool> m_noSegmentationHCalUsed{this, "noSegmentationHCal", true,
                                                 "HCal Barrel readout without DD4hep eta-phi segmentation used."};
  /// Seed threshold in sigma
  Gaudi::Property<int> m_seedSigma{this, "seedSigma", 4, "number of sigma in noise threshold"};
  /// Neighbour threshold in sigma
  Gaudi::Property<int> m_neighbourSigma{this, "neighbourSigma", 2, "number of sigma in noise threshold"};
  /// Last neighbour threshold in sigma
  Gaudi::Property<int> m_lastNeighbourSigma{this, "lastNeighbourSigma", 0, "number of sigma in noise threshold"};
  /// General decoder to encode the calorimeter sub-system to determine which positions tool to use
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder = new dd4hep::DDSegmentation::BitFieldCoder("system:4");
};
#endif /* RECCALORIMETER_CALOTOPOCLUSTER_H */
