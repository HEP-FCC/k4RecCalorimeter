#ifndef RECFCCEECALORIMETER_CREATEFCCEECALOXTALKNEIGHBOURS_H
#define RECFCCEECALORIMETER_CREATEFCCEECALOXTALKNEIGHBOURS_H

// Gaudi
#include "GaudiKernel/Service.h"

// k4FWCore
#include "k4FWCore/DataHandle.h"
#include "k4Interface/ICaloCreateMap.h"
#include "k4Interface/ICellPositionsTool.h"
#include "k4Interface/IGeoSvc.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"
#include "detectorCommon/xtalk_neighbors_moduleThetaMergedSegmentation.h"

// DD4hep
#include "DD4hep/Readout.h"
#include "DD4hep/Volumes.h"
#include "DDSegmentation/Segmentation.h"

// ROOT
#include "TGeoManager.h"

class IGeoSvc;

/** @class CreateFCCeeCaloXTalkNeighbours
 *
 *  Service building a map of crosstalk neighbours for all existing cells in the geometry.
 *  Only applicable to the ALLEGRO ECAL barrel with the theta-merged segmentation
 *
 *  @author Zhibo Wu
 */

class CreateFCCeeCaloXTalkNeighbours : public extends1<Service, ICaloCreateMap> {
public:
  /// Standard constructor
  explicit CreateFCCeeCaloXTalkNeighbours(const std::string& aName, ISvcLocator* aSL);
  /// Standard destructor
  virtual ~CreateFCCeeCaloXTalkNeighbours();
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
      this, "readoutNames", {"ECalBarrelModuleThetaMerged", "BarHCal_Readout_phitheta"}};
  /// Name of the fields describing the segmented volume
  Gaudi::Property<std::vector<std::string>> m_fieldNamesSegmented{this, "systemNames", {"system", "system"}};
  /// Values of the fields describing the segmented volume
  Gaudi::Property<std::vector<int>> m_fieldValuesSegmented{this, "systemValues", {4, 8}};
  /// Names of the active volume in geometry along radial axis (e.g. layer), the others are "module" or "phi", "theta"
  Gaudi::Property<std::vector<std::string>> m_activeFieldNamesSegmented{this, "activeFieldNames", {"layer", "layer"}};
  /// Number of layers in the segmented volumes
  Gaudi::Property<std::vector<unsigned int>> m_activeVolumesNumbersSegmented{this, "activeVolumesNumbers", {12, 13}};
  // Theta ranges of layers in the segmented volumes
  Gaudi::Property<std::vector<std::vector<double>>> m_activeVolumesTheta{this, "activeVolumesTheta"};

  // System ID of ECAL and HCAL barrels
  Gaudi::Property<uint> m_ecalBarrelSysId{this, "ecalBarrelSysId", 4};
  Gaudi::Property<uint> m_hcalBarrelSysId{this, "hcalBarrelSysId", 8};

  // For combination of barrels: flag if ECal and HCal barrels should be merged
  Gaudi::Property<bool> m_connectBarrels{this, "connectBarrels", true};
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
};

#endif /* RECFCCEECALORIMETER_CREATEFCCEECALOXTALKNEIGHBOURS_H */
