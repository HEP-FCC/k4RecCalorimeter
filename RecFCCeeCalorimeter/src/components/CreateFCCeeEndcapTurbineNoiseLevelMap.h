#ifndef RECCALORIMETER_CREATEFCCEEENDCAPTURBINENOISELEVELMAP_H
#define RECCALORIMETER_CREATEFCCEEENDCAPTURBINENOISELEVELMAP_H

// Gaudi
#include "GaudiKernel/Service.h"
#include "k4Interface/ICaloCreateMap.h"
#include "k4Interface/INoiseConstTool.h"
#include "k4Interface/ICellPositionsTool.h"

class IGeoSvc;

/** @class CreateFCCeeEndcapTurbineNoiseLevelMap
 *
 *  Service building a map from cellIds to noise level per cell for the
 *  endcap turbine calorimeter.  Values are just dummies at this point.
 *
 *  @author Erich Varnes
 */

class CreateFCCeeEndcapTurbineNoiseLevelMap : public extends1<Service, ICaloCreateMap> {
public:
  /// Standard constructor
  explicit CreateFCCeeEndcapTurbineNoiseLevelMap(const std::string& aName, ISvcLocator* aSL);
  /// Standard destructor
  virtual ~CreateFCCeeEndcapTurbineNoiseLevelMap();
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

 
  /// Names of the detector readout for volumes with eta-phi segmentation
  Gaudi::Property<std::vector<std::string>> m_readoutNamesSegmented{this, "readoutNamesSegmented", {"ECalEndcapTurbine"}};
  /// Name of the fields describing the segmented volume
  Gaudi::Property<std::vector<std::string>> m_fieldNamesSegmented{this, "systemNamesSegmented", {"system"}};
  /// Values of the fields describing the segmented volume
  Gaudi::Property<std::vector<int>> m_fieldValuesSegmented{this, "systemValuesSegmented", {5}};
  /// Names of the active volume in geometry 
  Gaudi::Property<std::vector<std::string>> m_activeFieldNamesSegmented{this, "activeFieldNamesSegmented", {"layer"}};
  /// Number of layers in the segmented volume
  Gaudi::Property<std::vector<unsigned int>> m_activeVolumesNumbersSegmented{this, "activeVolumesNumbers", {8}};

  /// Name of output file
  std::string m_outputFileName;
};

#endif /* RECALORIMETER_CREATEFCCEEENDCAPTURBINENOISELEVELMAP_H */
