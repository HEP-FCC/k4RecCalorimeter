#ifndef RECFCCEECALORIMETER_CALOTOWERTOOLFCCEE_H
#define RECFCCEECALORIMETER_CALOTOWERTOOLFCCEE_H

// from Gaudi
#include "GaudiKernel/AlgTool.h"

// k4FWCore
#include "k4FWCore/DataHandle.h"
#include "k4Interface/ITowerToolThetaModule.h"
class IGeoSvc;

namespace dd4hep {
namespace DDSegmentation {
  class Segmentation;
  class BitFieldCoder;
} // namespace DDSegmentation
} // namespace dd4hep

// edm4hep
namespace edm4hep {
class CalorimeterHitCollection;
class CalorimeterHit;
class Cluster;
} // namespace edm4hep

/** @class CaloTowerToolFCCee Reconstruction/RecFCCeeCalorimeter/src/components/CaloTowerToolFCCee.h
 *
 *  Tool building the calorimeter towers for the sliding window algorithm.
 *  This tool runs over all calorimeter systems (ECAL barrel, HCAL barrel + extended barrel, calorimeter endcaps,
 * forward calorimeters). If not all systems are available or not wanted to be used, create an empty collection using
 * CreateDummyCellsCollection algorithm.
 *  Towers are built of cells in theta-phi, summed over all radial layers.
 *  A tower contains all cells within certain theta and phi (tower size: '\b deltaThetaTower', '\b deltaPhiTower').
 *
 *  @author Anna Zaborowska
 *  @author Jana Faltova
 *  @Modified by Tong Li, for Theta-Module Merged readouts in FCCee
 */

class CaloTowerToolFCCee : public AlgTool, virtual public ITowerToolThetaModule {
public:
  CaloTowerToolFCCee(const std::string& type, const std::string& name, const IInterface* parent);
  virtual ~CaloTowerToolFCCee() = default;
  /**  Initialize.
   *   @return status code
   */
  virtual StatusCode initialize() final;
  /**  Finalize.
   *   @return status code
   */
  virtual StatusCode finalize() final;
  /**  Find number of calorimeter towers.
   *   Number of towers in phi is calculated from full azimuthal angle (2 pi) and the size of tower in phi ('\b
   * deltaPhiTower').
   *   Number of towers in theta is calculated from maximum detector theta ('\b thetaMax`) and the size of tower in
   * theta ('\b deltaThetaTower').
   *   @param[out] nTheta Number of towers in theta.
   *   @param[out] nPhi Number of towers in phi.
   */
  virtual void towersNumber(int& nTheta, int& nPhi) final;
  /**  Build calorimeter towers.
   *   Tower is defined by a segment in theta and phi, with the energy from all layers (no r segmentation).
   *   @param[out] aTowers Calorimeter towers.
   *   @param[in] fillTowersCells Whether to fill maps of cells into towers, for later use in attachCells
   *   @return Size of the cell collection.
   */
  virtual uint buildTowers(std::vector<std::vector<float>>& aTowers, bool fillTowersCells = true) final;
  /**  Get the map of cells contained within a tower.
   *   @return Map of cells in a tower
   */
  virtual std::map<std::pair<uint, uint>, std::vector<edm4hep::CalorimeterHit>> cellsInTowers() const final;
  /**  Get the tower IDs in theta.
   *   @param[in] aTheta Position of the calorimeter cell in theta
   *   @return ID (theta) of a tower
   */
  virtual uint idTheta(float aTheta) const final;
  /**  Get the tower IDs in phi.
   *   @param[in] aPhi Position of the calorimeter cell in phi
   *   @return ID (phi) of a tower
   */
  virtual uint idPhi(float aPhi) const final;
  /**  Get the theta position of the centre of the tower.
   *   @param[in] aIdTheta ID (theta) of a tower
   *   @return Position of the centre of the tower
   */
  virtual float theta(int aIdTheta) const final;
  /**  Get the phi position of the centre of the tower.
   *   @param[in] aIdPhi ID (phi) of a tower
   *   @return Position of the centre of the tower
   */
  virtual float phi(int aIdPhi) const final;
  /**  Find cells belonging to a cluster.
   *   @param[in] aTheta Position of the middle tower of a cluster in theta
   *   @param[in] aPhi Position of the middle tower of a cluster in phi
   *   @param[in] aHalfThetaFinal Half size of cluster in theta (in units of tower size). Cluster size is
   * 2*aHalfThetaFinal+1
   *   @param[in] aHalfPhiFinal Half size of cluster in phi (in units of tower size). Cluster size is 2*aHalfPhiFinal+1
   *   @param[out] aEdmCluster Cluster where cells are attached to
   */
  virtual void attachCells(float aTheta, float aPhi, uint aHalfThetaFinal, uint aHalfPhiFinal,
                           edm4hep::MutableCluster& aEdmCluster, edm4hep::CalorimeterHitCollection* aEdmClusterCells,
                           bool aEllipse = false) final;

private:
  /// Type of the segmentation
  enum class SegmentationType { kWrong, kModuleTheta, kMulti, kPhiTheta, kHCalPhiTheta, kHCalPhiRow, kEndcapTurbine };
  /**  Correct way to access the neighbour of the phi tower, taking into account
   * the full coverage in phi.
   *   Full coverage means that first tower in phi, with ID = 0 is a direct
   * neighbour of the last tower in phi with ID = m_nPhiTower - 1).
   *   @param[in] aIPhi requested ID of a phi tower, may be < 0 or >=m_nPhiTower
   *   @return ID of a tower - shifted and corrected (in [0, m_nPhiTower) range)
   */
  uint phiNeighbour(int aIPhi) const;
  /**  This is where the cell info is filled into towers
   *   @param[in] aTowers Calorimeter towers.
   *   @param[in] aCells Calorimeter cells collection.
   *   @param[in] fillTowerCells If true, make a list of the cells in each tower
   */
  void CellsIntoTowers(std::vector<std::vector<float>>& aTowers, const edm4hep::CalorimeterHitCollection* aCells,
                       bool fillTowersCells);
  /**  Find the maximum phi, theta covered by a readout
   *   @param[in] aReadoutName Readout name to be checked for maximum phi, theta
   *   @param[out] phiThetaPair  Values of the maximum phi and theta
   */
  StatusCode retrievePhiThetaExtrema(std::string aReadoutName, std::pair<double, double>& phiThetaPair);
  /**  Check if the readout name exists. If so, it returns the segmentation.
   *   @param[in] aReadoutName Readout name to be retrieved
   */
  std::pair<dd4hep::DDSegmentation::Segmentation*, SegmentationType> retrieveSegmentation(std::string aReadoutName);
  /// Handle for electromagnetic barrel cells (input collection)
  mutable DataHandle<edm4hep::CalorimeterHitCollection> m_ecalBarrelCells{"ecalBarrelCells", Gaudi::DataHandle::Reader,
                                                                          this};
  /// Handle for ecal endcap calorimeter cells (input collection)
  mutable DataHandle<edm4hep::CalorimeterHitCollection> m_ecalEndcapCells{"ecalEndcapCells", Gaudi::DataHandle::Reader,
                                                                          this};
  /// Handle for ecal forward calorimeter cells (input collection)
  mutable DataHandle<edm4hep::CalorimeterHitCollection> m_ecalFwdCells{"ecalFwdCells", Gaudi::DataHandle::Reader, this};
  /// Handle for hadronic barrel cells (input collection)
  mutable DataHandle<edm4hep::CalorimeterHitCollection> m_hcalBarrelCells{"hcalBarrelCells", Gaudi::DataHandle::Reader,
                                                                          this};
  /// Handle for hadronic extended barrel cells (input collection)
  mutable DataHandle<edm4hep::CalorimeterHitCollection> m_hcalExtBarrelCells{"hcalExtBarrelCells",
                                                                             Gaudi::DataHandle::Reader, this};
  /// Handle for hcal endcap calorimeter cells (input collection)
  mutable DataHandle<edm4hep::CalorimeterHitCollection> m_hcalEndcapCells{"hcalEndcapCells", Gaudi::DataHandle::Reader,
                                                                          this};
  /// Handle for hcal forward calorimeter cells (input collection)
  mutable DataHandle<edm4hep::CalorimeterHitCollection> m_hcalFwdCells{"hcalFwdCells", Gaudi::DataHandle::Reader, this};

  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;
  /// Name of the electromagnetic barrel readout
  Gaudi::Property<std::string> m_ecalBarrelReadoutName{this, "ecalBarrelReadoutName", "",
                                                       "name of the ecal barrel readout"};
  /// Name of the ecal endcap calorimeter readout
  Gaudi::Property<std::string> m_ecalEndcapReadoutName{this, "ecalEndcapReadoutName", "",
                                                       "name of the ecal endcap readout"};
  /// Name of the ecal forward calorimeter readout
  Gaudi::Property<std::string> m_ecalFwdReadoutName{this, "ecalFwdReadoutName", "", "name of the ecal fwd readout"};
  /// Name of the hadronic barrel readout
  Gaudi::Property<std::string> m_hcalBarrelReadoutName{this, "hcalBarrelReadoutName", "",
                                                       "name of the hcal barrel readout"};
  /// Name of the hadronic extended barrel readout
  Gaudi::Property<std::string> m_hcalExtBarrelReadoutName{this, "hcalExtBarrelReadoutName", "",
                                                          "name of the hcal extended barrel readout"};
  /// Name of the hcal endcap calorimeter readout
  Gaudi::Property<std::string> m_hcalEndcapReadoutName{this, "hcalEndcapReadoutName", "",
                                                       "name of the hcal endcap readout"};
  /// Name of the hcal forward calorimeter readout
  Gaudi::Property<std::string> m_hcalFwdReadoutName{this, "hcalFwdReadoutName", "", "name of the hcal fwd readout"};
  /// Type of segmentation of the electromagnetic barrel
  SegmentationType m_ecalBarrelSegmentationType;
  /// Type of segmentation of the ecal endcap calorimeter
  SegmentationType m_ecalEndcapSegmentationType;
  /// Type of segmentation of the ecal forward calorimeter
  SegmentationType m_ecalFwdSegmentationType;
  /// Type of segmentation of the hadronic barrel
  SegmentationType m_hcalBarrelSegmentationType;
  /// Type of segmentation of the hadronic extended barrel
  SegmentationType m_hcalExtBarrelSegmentationType;
  /// Type of segmentation of the hcal endcap calorimeter
  SegmentationType m_hcalEndcapSegmentationType;
  /// Type of segmentation of the hcal forward calorimeter
  SegmentationType m_hcalFwdSegmentationType;
  /// decoder: only for barrel
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

  /// Maximum theta of detector
  float m_thetaMax;
  /// Maximum phi of the detector
  float m_phiMax;
  /// Size of the tower in theta
  Gaudi::Property<float> m_deltaThetaTower{this, "deltaThetaTower", 0.01, "Size of the tower in theta"};
  /// Size of the tower in phi
  Gaudi::Property<float> m_deltaPhiTower{this, "deltaPhiTower", 0.01, "Size of the tower in phi"};
  /// number of towers in theta (calculated from m_deltaThetaTower and m_thetaMax)
  int m_nThetaTower;
  /// Number of towers in phi (calculated from m_deltaPhiTower)
  int m_nPhiTower;
  /// map to cells contained within a tower so they can be attached to a reconstructed cluster
  std::map<std::pair<uint, uint>, std::vector<edm4hep::CalorimeterHit>> m_cellsInTowers;
  /// Use only a part of the calorimeter (in depth)
  Gaudi::Property<bool> m_useHalfTower{this, "halfTower", false, "Use half tower"};
  /// Max layer
  Gaudi::Property<uint> m_max_layer{
      this, "max_layer", 12,
      "Specify which radial layer are used. The condition is 'if(cellLayer > m_max_layer) skip this cell'."};
};

#endif /* RECFCCEECALORIMETER_CALOTOWERTOOLFCCEE_H */
