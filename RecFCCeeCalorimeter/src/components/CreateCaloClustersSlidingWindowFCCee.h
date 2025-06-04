#ifndef RECFCCEECALORIMETER_CREATECALOCLUSTERSSLIDINGWINDOWFCCEE_H
#define RECFCCEECALORIMETER_CREATECALOCLUSTERSSLIDINGWINDOWFCCEE_H

// GAUDI
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"

// k4FWCore
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/MetaDataHandle.h"
#include "k4Interface/ITowerToolThetaModule.h"

// edm4hep
namespace edm4hep {
class ClusterCollection;
}

/** @class CreateCaloClustersSlidingWindowFCCee
 *
 *  Algorithm for creating calorimeter clusters from cells.
 *
 *  Reconstruction is performed for the calorimeter with readout '\b readoutName' and volume ID retrieved using '\b
 *fieldNames' and '\b fieldValues'.
 *
 *  Sliding window algorithm:
 *  1. Create calorimeter towers.
 *     A tower contains all cells within certain theta and phi (tower size: '\b deltaThetaTower', '\b deltaPhiTower').
 *     Currently there is no support of cell splitting, so each cell should be completely inside the tower
 *  2. Find local maxima.
 *     Local maxima are found using the sliding window of a fixed size in phi x theta ('\b nThetaWindow' '\b nPhiWindow'
 *in units of tower size). If a local max is found and its energy is above threshold ('\b energyThreshold'), it is added
 *to the preclusters list. Each precluster contains the barycentre position and the transverse energy. Position is
 *recalculated using the window size in theta x phi ('\b nThetaPosition', '\b nPhiPosition') that may be smaller than
 *the sliding window to reduce the noise influence. Both windows are centred at the same tower. The energy of the
 *precluster is the energy calculated using the sliding window.
 *  3. Remove duplicates.
 *     If two pre-clusters are found next to each other (within window '\b nThetaDuplicates', '\b nPhiDuplicates'), the
 *pre-cluster with lower energy is removed.
 *     Currently there is no support on energy sharing between clusters, so if duplicate window is smaller than
 *sliding window, some towers may be taken twice (instead of the weighted energy).
 *  4. Build clusters.
 *     Clusters are created using the window of a fixed size in phi x theta ('\b nThetaFinal' '\b nPhiFinal' in units of
 *tower size) around the barycentre position.
 *     For each cluster the cell collection is searched and all those inside the cluster are attached.
 *
 *  Note: Sliding window performs well for electrons/gamma reconstruction. Topological clusters should be better for
 *jets.
 *
 *  @author Jana Faltova
 *  @author Anna Zaborowska
 *  @author Tong Li, for Theta-Module Merged readouts in FCCee
 *  @author Giovanni Marchiori
 */

class CreateCaloClustersSlidingWindowFCCee : public Gaudi::Algorithm {
public:
  CreateCaloClustersSlidingWindowFCCee(const std::string& name, ISvcLocator* svcLoc);
  /**  Initialize.
   *   @return status code
   */
  StatusCode initialize();
  /**  Execute.
   *   Perform the sliding window algorithm and build clusters.
   *   @return status code
   */
  StatusCode execute(const EventContext&) const;
  /**  Finalize.
   *   @return status code
   */
  StatusCode finalize();

private:
  // precluster
  struct precluster {
    float transEnergy;
    float theta;
    float phi;
    float X;
    float Y;
    float Z;
  };

  /**  Correct way to access the neighbour of the phi tower, taking into account the full coverage in phi.
   *   Full coverage means that first tower in phi, with ID = 0 is a direct neighbour
   *   of the last tower in phi with ID = m_nPhiTower - 1).
   *   @param[in] aIPhi requested ID of a phi tower, may be < 0 or >= m_nPhiTower
   *   @return  ID of a tower - shifted and corrected (in [0, m_nPhiTower) range)
   */
  unsigned int phiNeighbour(int aIPhi) const;
  /// Handle for calo clusters (output collection)
  mutable k4FWCore::DataHandle<edm4hep::ClusterCollection> m_clusters{"calo/clusters", Gaudi::DataHandle::Writer, this};
  /// Handle for calo cluster cells (output collection)
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_clusterCells{"calo/clusterCells",
                                                                                 Gaudi::DataHandle::Writer, this};
  /// Output collection metadata handles (saving a map of ID:collection does not work)
  k4FWCore::MetaDataHandle<std::vector<int>> m_caloIDsMetaData{m_clusters, "inputSystemIDs", Gaudi::DataHandle::Writer};
  k4FWCore::MetaDataHandle<std::vector<std::string>> m_cellsMetaData{m_clusters, "inputCellCollections",
                                                                     Gaudi::DataHandle::Writer};
  /// Handle for the tower building tool
  mutable ToolHandle<ITowerToolThetaModule> m_towerTool;
  /// Calorimeter towers
  mutable std::vector<std::vector<float>> m_towers;
  /// Vector of pre-clusters
  mutable std::vector<precluster> m_preClusters;
  /// number of towers in theta (calculated from m_deltaThetaTower and the theta size of the first layer)
  int m_nThetaTower;
  /// Number of towers in phi (calculated from m_deltaPhiTower)
  int m_nPhiTower;
  /// Size of the window in theta for pre-clusters (in units of tower size)
  Gaudi::Property<int> m_nThetaWindow{this, "nThetaWindow", 5};
  /// Size of the window in phi for pre-clusters (in units of tower size)
  Gaudi::Property<int> m_nPhiWindow{this, "nPhiWindow", 15};
  /// Size of the window in theta for cluster position calculation (in units of tower size)
  Gaudi::Property<int> m_nThetaPosition{this, "nThetaPosition", 3};
  /// Size of the window in phi for cluster position calculation (in units of tower size)
  Gaudi::Property<int> m_nPhiPosition{this, "nPhiPosition", 3};
  /// Size of the window in theta for the overlap removal (in units of tower size)
  Gaudi::Property<int> m_nThetaDuplicates{this, "nThetaDuplicates", 2};
  /// Size of the window in phi for the overlap removal (in units of tower size)
  Gaudi::Property<int> m_nPhiDuplicates{this, "nPhiDuplicates", 2};
  /// Energy threshold for cluster finding
  Gaudi::Property<float> m_energyThreshold{this, "energyThreshold", 10};
  /// Size of the window in theta for the final cluster building (in units of tower size)
  Gaudi::Property<int> m_nThetaFinal{this, "nThetaFinal", 5};
  /// Size of the window in phi for the final cluster building (in units of tower size)
  Gaudi::Property<int> m_nPhiFinal{this, "nPhiFinal", 15};
  /// Fraction of the energy threshold used in the position calculation
  Gaudi::Property<float> m_energyThresholdFraction{this, "energyThresholdFraction", 0.25};
  /// Flag for the application of the correction for energy sharing between clusters
  Gaudi::Property<bool> m_energySharingCorrection{this, "energySharingCorrection", false};
  /// Flag for the ellipse used in the final cluster instead of the rectangle
  Gaudi::Property<bool> m_ellipseFinalCluster{this, "ellipse", false};
  /// Flag if a new output cell collection of clustered cells should be created
  Gaudi::Property<bool> m_createClusterCellCollection{this, "createClusterCellCollection", false};
};

#endif /* RECFCCEECALORIMETER_CREATECALOCLUSTERSSLIDINGWINDOWFCCEE_H */
