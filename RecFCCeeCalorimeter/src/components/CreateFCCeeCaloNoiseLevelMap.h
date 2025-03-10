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
 *  or phi-theta (e.g. HCal barrel).
 *
 *  @author Coralie Neubueser
 *  @author Giovanni Marchiori
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

  /// ideally could replace with vector of handles (so that we can also have ecal endcap and so on)
  /// Handle for the cells noise tool in ECal barrel
  ToolHandle<INoiseConstTool> m_ecalBarrelNoiseTool{"ReadNoiseFromFileTool", this};
  /// Handle for the cells noise tool in HCal
  ToolHandle<INoiseConstTool> m_ecalEndcapNoiseTool{"ReadNoiseFromFileTool", this};
  /// Handle for the cells noise tool in HCal
  ToolHandle<INoiseConstTool> m_hcalBarrelNoiseTool{"ReadNoiseFromFileTool", this};
  /// Handle for the cells noise tool in HCal
  ToolHandle<INoiseConstTool> m_hcalEndcapNoiseTool{"ReadNoiseFromFileTool", this};

  // System ID of ECAL and HCAL detectors
  Gaudi::Property<uint> m_ecalBarrelSysId{this, "ecalBarrelSysId", 4};
  Gaudi::Property<uint> m_ecalEndcapSysId{this, "ecalEndcapSysId", 5};
  Gaudi::Property<uint> m_hcalBarrelSysId{this, "hcalBarrelSysId", 8};
  Gaudi::Property<uint> m_hcalEndcapSysId{this, "hcalEndcapSysId", 9};

  /// Names of the detector readout for the volumes
  Gaudi::Property<std::vector<std::string>> m_readoutNamesSegmented{this, "readoutNames", {"ECalBarrelModuleThetaMerged", "ECalEndcapTurbine", "HCalBarrelReadout","HCalEndcapReadout"}};
  /// Name of the fields describing the segmented volume
  Gaudi::Property<std::vector<std::string>> m_fieldNamesSegmented{this, "systemNames", {"system", "system", "system", "system"}};
  /// Values of the fields describing the segmented volume
  Gaudi::Property<std::vector<int>> m_fieldValuesSegmented{this, "systemValues", {4, 5, 8, 9}};
  /// Names of the active volume in geometry along radial axis (e.g. layer), the others are "module" or "phi", "theta"
  Gaudi::Property<std::vector<std::string>> m_activeFieldNamesSegmented{this, "activeFieldNames", {"layer", "layer", "layer", "layer"}};
  /// Number of layers in the segmented volumes
  Gaudi::Property<std::vector<unsigned int>> m_activeVolumesNumbersSegmented{this, "activeVolumesNumbers", {11, 98, 13, 37}};
  // Theta ranges of layers in the segmented volumes (if needed)
  Gaudi::Property<std::vector<std::vector<double>>> m_activeVolumesTheta{this, "activeVolumesTheta"};

  /// Name of output file
  std::string m_outputFileName;
};

#endif /* RECALORIMETER_CREATEFCCEECALONOISELEVELMAP_H */
