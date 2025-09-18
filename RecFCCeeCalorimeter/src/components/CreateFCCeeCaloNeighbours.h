#ifndef RECFCCEECALORIMETER_CREATEFCCEECALONEIGHBOURS_H
#define RECFCCEECALORIMETER_CREATEFCCEECALONEIGHBOURS_H

// Gaudi
#include "GaudiKernel/Service.h"

// k4FWCore
#include "k4FWCore/DataHandle.h"
#include "k4Interface/ICaloCreateMap.h"
#include "k4Interface/ICellPositionsTool.h"
#include "k4Interface/IGeoSvc.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"

// DD4hep
#include "DD4hep/Readout.h"
#include "DD4hep/Volumes.h"
#include "DDSegmentation/Segmentation.h"

// ROOT
#include "TGeoManager.h"

class IGeoSvc;
class TH1F;

/** @class CreateFCCeeCaloNeighbours
 *
 *  Service building a map of neighbours for all existing cells in the geometry.
 *  The volumes for which the neighbour map is created can be either segmented in theta-module (e.g. ECal inclined),
 *
 *  @author Giovanni Marchiori
 */

class CreateFCCeeCaloNeighbours : public extends1<Service, ICaloCreateMap> {
public:
  /// Standard constructor
  explicit CreateFCCeeCaloNeighbours(const std::string& aName, ISvcLocator* aSL);
  /// Standard destructor
  virtual ~CreateFCCeeCaloNeighbours();
  /**  Initialize the map creator service.
   *   @return status code
   */
  virtual StatusCode initialize() final;
  /**  Finalize the map creator service.
   *   @return status code
   */
  virtual StatusCode finalize() final;

private:
  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;

  /// Names of the detector readout for the volumes
  Gaudi::Property<std::vector<std::string>> m_readoutNamesSegmented{
    this, "readoutNames", {"ECalBarrelModuleThetaMerged", "ECalEndcapTurbine", "HCalBarrelReadout", "HCalEndcapReadout"}};
  /// Name of the fields describing the segmented volume
  Gaudi::Property<std::vector<std::string>> m_fieldNamesSegmented{this, "systemNames", {"system", "system", "system", "system"}};
  /// Values of the fields describing the segmented volume
  Gaudi::Property<std::vector<int>> m_fieldValuesSegmented{this, "systemValues", {4, 5, 8, 9}};
  /// Names of the active volume in geometry along radial axis (e.g. layer), the others are "module" or "phi", "theta"
  Gaudi::Property<std::vector<std::string>> m_activeFieldNamesSegmented{
    this, "activeFieldNames", {"layer", "layer", "layer", "layer"}};
  /// Number of layers in the segmented volumes
  Gaudi::Property<std::vector<unsigned int>> m_activeVolumesNumbersSegmented{
    this, "activeVolumesNumbers", {11, 98, 13, 37}};
  // Theta ranges of layers in the segmented volumes
  Gaudi::Property<std::vector<std::vector<double>>> m_activeVolumesTheta{this, "activeVolumesTheta"};
  /// Whether to consider diagonal cells as neighbours or not
  Gaudi::Property<bool> m_includeDiagonalCells{
      this, "includeDiagonalCells", false,
      "If True will consider also diagonal neighbours in volumes with theta-module segmentation"};
  /// Whether to consider diagonal cells as neighbours or not for HCal phi-theta segmentation
  Gaudi::Property<bool> m_includeDiagonalCellsHCal{
      this, "includeDiagonalCellsHCal", false,
      "If True will consider also diagonal neighbours in volumes with phi-theta segmentation"};

  // System ID of ECAL and HCAL barrels
  Gaudi::Property<uint> m_ecalBarrelSysId{this, "ecalBarrelSysId", 4};
  Gaudi::Property<uint> m_ecalEndcapSysId{this, "ecalEndcapSysId", 5};  
  Gaudi::Property<uint> m_hcalBarrelSysId{this, "hcalBarrelSysId", 8};
  // System ID of HCAL endcap
  Gaudi::Property<uint> m_hcalEndcapSysId{this, "hcalEndcapSysId", 9};

  // For combination of barrels: flag if ECal and HCal barrels should be merged
  Gaudi::Property<bool> m_connectBarrels{this, "connectBarrels", true};
  // For combination of HCal Barrel and Endcap: flag if HCal Barrel and Endcap should be merged
  Gaudi::Property<bool> m_connectHCal{this, "connectHCal", true};
  // For combination of ECal Barrel and Endcap: flag if ECal Barrel and Endcap should be merged
  Gaudi::Property<bool> m_connectECal{this, "connectECal", true};
  // For combination of barrels: size of HCal cell along z axis
  Gaudi::Property<double> m_hCalZSize{this, "hCalZsize", 18};
  // For combination of barrels: offset of HCal detector in z (lower edge)
  Gaudi::Property<double> m_hCalZOffset{this, "hCalZoffset", -4590};
  // For combination of barrels: HCal inner radius for calculation of eta from z ???
  Gaudi::Property<double> m_hCalRinner{this, "hCalRinner", 2850};
  // For combination of barrels: offset of HCal modules in phi (lower edge)
  Gaudi::Property<double> m_hCalPhiOffset{this, "hCalPhiOffset"};

  /// Name of output file
  std::string m_outputFileName;

  // Lookup tables/histograms for cell positions in the barrel and endcap EM
  // calorimeters. The vectors of histograms function as lookups for the module ID  // vs phi, since each histogram has bins in phi with the content set to the
  // module ID associated with the phi range of that bin.  Each histogram in the
  // vector represents one layer of the barrel or one rho bin of the endcap.
  // The vectors of Float_t, on the other hand, store the central value of
  // theta or phi for a given layer (barrel) or rho bin (endcap).
  // By filling these vectors at the start of processing, the code can perform
  // the spatial matching between the barrel and endcap cells without
  // re-calculating positions inside the double loop over cells that occurs when
  // looking for neighbors.
  
  std::vector< TH1F* > m_EMB_h_module_vs_phi;
  std::vector< std::vector<Float_t>> m_EMB_phi_lookup;
  std::vector<Float_t> m_EMB_theta_lookup;

  std::vector< TH1F* > m_EMEC_h_module_vs_phi_pos;
  std::vector< TH1F* > m_EMEC_h_module_vs_phi_neg;
  std::vector< std::vector<Float_t>> m_EMEC_pos_phi_lookup;
  std::vector< std::vector<Float_t>> m_EMEC_neg_phi_lookup;
  std::vector<Float_t> m_EMEC_theta_lookup;

  StatusCode initialize_lookups();
  
};

#endif /* RECFCCEECALORIMETER_CREATEFCCEECALONEIGHBOURS_H */
