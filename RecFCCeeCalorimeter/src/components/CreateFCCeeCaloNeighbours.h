#ifndef RECFCCEECALORIMETER_CREATEFCCEECALONEIGHBOURS_H
#define RECFCCEECALORIMETER_CREATEFCCEECALONEIGHBOURS_H

// Gaudi
#include "GaudiKernel/Service.h"

// k4FWCore
#include "k4Interface/ICaloCreateMap.h"
#include "k4Interface/IGeoSvc.h"
#include "k4FWCore/DataHandle.h"
#include "k4Interface/ICellPositionsTool.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"

// DD4hep
#include "DD4hep/Readout.h"
#include "DD4hep/Volumes.h"
#include "DDSegmentation/Segmentation.h"

// ROOT
#include "TGeoManager.h"

class IGeoSvc;

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

  /// Names of the detector readout for volumes with theta-module segmentation
  Gaudi::Property<std::vector<std::string>> m_readoutNamesSegmented{this, "readoutNames", {"ECalBarrelModuleThetaMerged"}};
  /// Name of the fields describing the segmented volume
  Gaudi::Property<std::vector<std::string>> m_fieldNamesSegmented{this, "systemNames", {"system"}};
  /// Values of the fields describing the segmented volume
  Gaudi::Property<std::vector<int>> m_fieldValuesSegmented{this, "systemValues", {5}};
  /// Names of the active volume in geometry along radial axis (e.g. layer), the others are "module"/"phi", "theta"
  Gaudi::Property<std::vector<std::string>> m_activeFieldNamesSegmented{this, "activeFieldNames", {"layer"}};
  /// Number of layers in the segmented volume
  Gaudi::Property<std::vector<unsigned int>> m_activeVolumesNumbersSegmented{this, "activeVolumesNumbers", {12}};
  // Radii of layers in the segmented volume
  Gaudi::Property<std::vector<double>> m_activeVolumesTheta{this, "activeVolumesTheta"};
  /// Whether to create the geant4 geometry or not
  Gaudi::Property<bool> m_includeDiagonalCells{this, "includeDiagonalCells", false, "If True will consider also diagonal neighbours in volumes with theta-module segmentation"};
  
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

#endif /* RECFCCEECALORIMETER_CREATEFCCEECALONEIGHBOURS_H */
