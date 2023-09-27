#ifndef RECCALORIMETER_CREATEFCCEECALONEIGHBOURS_H
#define RECCALORIMETER_CREATEFCCEECALONEIGHBOURS_H

// Gaudi
#include "GaudiKernel/Service.h"
#include "k4Interface/ICaloCreateMap.h"

#include "DetCommon/DetUtils.h"
#include "k4Interface/IGeoSvc.h"
#include "DetSegmentation/FCCSWGridPhiEta.h"
#include "DetSegmentation/FCCSWGridPhiTheta.h"
#include "DetSegmentation/FCCSWGridModuleThetaMerged.h"
#include "k4FWCore/DataHandle.h"
#include "k4Interface/ICellPositionsTool.h"

#include "DD4hep/Readout.h"
#include "DD4hep/Volumes.h"
#include "DDSegmentation/Segmentation.h"
#include "TGeoManager.h"

class IGeoSvc;

/** @class CreateFCCeeCaloNeighbours
 *
 *  Service building a map of neighbours for all existing cells in the geometry.
 *  The volumes for which the neighbour map is created can be either segmented in theta-module (e.g. ECal inclined),
 *  or can contain nested volumes (e.g. HCal barrel).
 *
 *  @author Anna Zaborowska
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

  /// Names of the detector readout for volumes with theta-module segmentation
  Gaudi::Property<std::vector<std::string>> m_readoutNamesSegmented{this, "readoutNamesModuleTheta", {"ECalBarrelModuleThetaMerged"}};
  /// Name of the fields describing the segmented volume
  Gaudi::Property<std::vector<std::string>> m_fieldNamesSegmented{this, "systemNamesModuleTheta", {"system"}};
  /// Values of the fields describing the segmented volume
  Gaudi::Property<std::vector<int>> m_fieldValuesSegmented{this, "systemValuesModuleTheta", {5}};
  /// Names of the active volume in geometry along radial axis (e.g. layer), the others are "module", "theta"
  Gaudi::Property<std::vector<std::string>> m_activeFieldNamesSegmented{this, "activeFieldNamesModuleTheta", {"layer"}};
  /// Number of layers in the segmented volume
  Gaudi::Property<std::vector<unsigned int>> m_activeVolumesNumbersSegmented{this, "activeVolumesNumbers", {12}};
  // Radii of layers in the segmented volume
  Gaudi::Property<std::vector<double>> m_activeVolumesTheta{this, "activeVolumesTheta"};

  /// Names of the detector readout for volumes with nested volume structure and no segmentation
  Gaudi::Property<std::vector<std::string>> m_readoutNamesNested{this, "readoutNamesVolumes"};
  /// Name of the field describing the nested volume
  Gaudi::Property<std::string> m_fieldNameNested{this, "systemNameNested"};
  /// Values of the fields describing the nested volume
  Gaudi::Property<std::vector<int>> m_fieldValuesNested{this, "systemValuesNested"};
  /// Names of the active volume in geometry: along radial axis, azimuthal angle, and along z axis
  Gaudi::Property<std::vector<std::string>> m_activeFieldNamesNested{
      this, "activeFieldNamesNested"};
  /// Names of the nested volumes - to retrieve the number of active volumes, need to correspond to m_activeFieldNamesNested
  Gaudi::Property<std::vector<std::string>> m_activeVolumeNamesNested{
      this,
      "activeVolumeNamesNested",
      {"layerVolume", "moduleVolume", "wedgeVolume"}};  // to find out number of volumes
  /// Name of output file
  std::string m_outputFileName;

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
};

#endif /* RECALORIMETER_CREATEFCCHHCALONEIGHBOURS_H */
