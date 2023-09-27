#ifndef RECCALORIMETER_CREATEFCCEECALONOISELEVELMAP_H
#define RECCALORIMETER_CREATEFCCEECALONOISELEVELMAP_H

// Gaudi
#include "GaudiKernel/Service.h"
#include "k4Interface/ICaloCreateMap.h"
#include "k4Interface/INoiseConstTool.h"
#include "k4Interface/ICellPositionsTool.h"

class IGeoSvc;

/** @class CreateFCCeeCaloNoiseLevelMap
 *
 *  Service building a map from cellIds to noise level per cell.
 *  The volumes for which the neighbour map is created can be either segmented in Module-Theta (e.g. ECal inclined),
 *  or can contain nested volumes (e.g. HCal barrel).
 *
 *  @author Coralie Neubueser
 */

class CreateFCCeeCaloNoiseLevelMap : public extends1<Service, ICaloCreateMap> {
public:
  /// Standard constructor
  explicit CreateFCCeeCaloNoiseLevelMap(const std::string& aName, ISvcLocator* aSL);
  /// Standard destructor
  virtual ~CreateFCCeeCaloNoiseLevelMap();
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

  /// Handle for the cells noise tool in ECal
  ToolHandle<INoiseConstTool> m_ecalBarrelNoiseTool{"ReadNoiseFromFileTool", this};
  Gaudi::Property<uint> m_ecalBarrelSysId{this, "ecalBarrelSysId", 5};
  /// Handle for the cells noise tool in HCal
  ToolHandle<INoiseConstTool> m_hcalBarrelNoiseTool{"ReadNoiseFromFileTool", this};
  Gaudi::Property<uint> m_hcalBarrelSysId{this, "hcalBarrelSysId", 8};

  /// Names of the detector readout for volumes with Module-Theta segmentation
  Gaudi::Property<std::vector<std::string>> m_readoutNamesSegmented{this, "readoutNamesModuleTheta", {"ECalBarrelModuleThetaMerged"}};
  /// Name of the fields describing the segmented volume
  Gaudi::Property<std::vector<std::string>> m_fieldNamesSegmented{this, "systemNamesModuleTheta", {"system"}};
  /// Values of the fields describing the segmented volume
  Gaudi::Property<std::vector<int>> m_fieldValuesSegmented{this, "systemValuesModuleTheta", {5}};
  /// Names of the active volume in geometry along radial axis (e.g. layer), the others are "module", "theta"
  Gaudi::Property<std::vector<std::string>> m_activeFieldNamesSegmented{this, "activeFieldNamesModuleTheta", {"layer"}};
  /// Number of layers in the segmented volume
  Gaudi::Property<std::vector<unsigned int>> m_activeVolumesNumbersSegmented{this, "activeVolumesNumbers", {8}};
  // Radii of layers in the segmented volume
  Gaudi::Property<std::vector<double>> m_activeVolumesTheta{this, "activeVolumesTheta"};

  /// Names of the detector readout for volumes with nested volume structure and no segmentation
  Gaudi::Property<std::vector<std::string>> m_readoutNamesNested{this, "readoutNamesVolumes", {"HCalBarrelReadout"}};
  /// Name of the field describing the nested volume
  Gaudi::Property<std::string> m_fieldNameNested{this, "systemNameNested", "system"};
  /// Values of the fields describing the nested volume
  Gaudi::Property<std::vector<int>> m_fieldValuesNested{this, "systemValuesNested", {8}};
  /// Names of the active volume in geometry: along radial axis, azimuthal angle, and along z axis
  Gaudi::Property<std::vector<std::string>> m_activeFieldNamesNested{
      this, "activeFieldNamesNested", {"layer", "module", "row"}};
  /// Names of the nested volumes - to retrieve the number of active volumes, need to correspond to
  /// m_activeFieldNamesNested
  Gaudi::Property<std::vector<std::string>> m_activeVolumeNamesNested{
      this,
      "activeVolumeNamesNested",
      {"layerVolume", "moduleVolume", "wedgeVolume"}};  // to find out number of volumes

  /// Name of output file
  std::string m_outputFileName;
};

#endif /* RECALORIMETER_CREATEFCCEECALONOISELEVELMAP_H */
