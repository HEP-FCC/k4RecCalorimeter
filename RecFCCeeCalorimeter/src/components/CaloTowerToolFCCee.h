#ifndef RECFCCEECALORIMETER_CALOTOWERTOOLFCCEE_H
#define RECFCCEECALORIMETER_CALOTOWERTOOLFCCEE_H

// from Gaudi
#include "GaudiKernel/AlgTool.h"

// k4FWCore
#include "k4FWCore/DataHandle.h"
#include "k4Interface/ITowerToolThetaModule.h"

#include <cmath>

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
 *  @author Tong Li: implement theta-based grid
 *  @author Giovanni Marchiori: cleanup, generalise
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
  /**  Correct way to obtain the phi index of a tower, taking into account
   * the phi periodicity.
   *   @param[in] aIPhi requested ID of a phi tower, may be < 0 or >=m_nPhiTower
   *   @return ID of a tower - shifted and corrected (in [0, m_nPhiTower) range)
   */
  uint phiIndexTower(int aIPhi) const;
  /**  This is where the cell info is filled into towers
   *   @param[in] aTowers Calorimeter towers.
   *   @param[in] aCells Calorimeter cells collection.
   *   @param[in] fillTowerCells If true, make a list of the cells in each tower
   *   @return number of clustered cells
   */
  uint CellsIntoTowers(std::vector<std::vector<float>>& aTowers, const edm4hep::CalorimeterHitCollection* aCells,
                       bool fillTowersCells);

  /// List of input cell collections
  Gaudi::Property<std::vector<std::string>> m_cellCollections{
    this, "cells", {}, "Names of CalorimeterHit collections to read"};

    /// The vector of input DataHandles for the input cell collections
  std::vector<DataHandle<edm4hep::CalorimeterHitCollection>*> m_cellCollectionHandles;

  /// Maximum theta of towers
  /// Can be left to pi, it won't hurt
  // (there will just be towers beyond the detector acceptance with zero energy)
  Gaudi::Property<float> m_thetaMax{this, "thetaMax", M_PI, "Maximum theta of towers"};
  /// Maximum theta of towers
  Gaudi::Property<float> m_thetaMin{this, "thetaMin", 0., "Minimum theta of towers"};
  /// Maximum phi of towers (note: phi is in -pi..pi)
  Gaudi::Property<float> m_phiMax{this, "phiMax", M_PI, "Maximum phi of towers"};
  /// Minimum phi of towers
  Gaudi::Property<float> m_phiMin{this, "phiMin", -M_PI, "Minimum phi of towers"};
  /// Size of the tower in theta (default: pi/314 ~ 0.01 rad)
  Gaudi::Property<float> m_deltaThetaTower{this, "deltaThetaTower", M_PI/314, "Size of the tower in theta"};
  /// Size of the tower in phi (default: pi/314 ~ 0.01 rad)
  Gaudi::Property<float> m_deltaPhiTower{this, "deltaPhiTower", M_PI/314, "Size of the tower in phi"};

  /// Number of calorimeters (e.g. ecal/hcal/muon) for which to calculate the
  /// total subdetector energy of the cluster and save it in the subdetectorEnergies vector
  /// Should be equal to max(caloID) where caloID is defined as calculated in
  /// CreatePositionedCaloCells
  /// Barrel and endcap are considered as part of the same subdetector
  /// If the property is 0, subdetector energy information will not be computed
  Gaudi::Property<unsigned int> m_nSubDetectors{this, "nSubDetectors", 0, "Number of calorimeter systems"};

  /// number of towers in theta
  int m_nThetaTower;
  /// Number of towers in phi
  int m_nPhiTower;
  /// map to cells contained within a tower so they can be attached to a reconstructed cluster
  std::map<std::pair<uint, uint>, std::vector<edm4hep::CalorimeterHit>> m_cellsInTowers;
};

#endif /* RECFCCEECALORIMETER_CALOTOWERTOOLFCCEE_H */
