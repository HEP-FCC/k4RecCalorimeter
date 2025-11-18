#ifndef RECFCCEECALORIMETER_CREATEFCCEECALOSIGNALSHAPES_H
#define RECFCCEECALORIMETER_CREATEFCCEECALOSIGNALSHAPES_H

// Gaudi
#include "GaudiKernel/Service.h"

// k4FWCore
#include "k4Interface/ICaloCreateMap.h"
#include "k4Interface/IGeoSvc.h"
#include "k4FWCore/DataHandle.h"
#include "k4Interface/ICellPositionsTool.h"

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

/** @class CreateFCCeeCaloSignalShapes
 *
 *  Service building a map of signal shapes for the ECAL detector.
 *  Only applicable to the ALLEGRO ECAL barrel with the theta-merged segmentation
 *
 *  @author Sahibjeet Singh
 */

class CreateFCCeeCaloSignalShapes : public extends1<Service, ICaloCreateMap> {
public:
  /// Standard constructor
  explicit CreateFCCeeCaloSignalShapes(const std::string& aName, ISvcLocator* aSL);
  /// Standard destructor
  virtual ~CreateFCCeeCaloSignalShapes();
  /**  Initialize the map creator service.
   *   @return status code
   */
  virtual StatusCode initialize() final;
  /**  Finalize the map creator service.
   *   @return status code
   */
  virtual StatusCode finalize() final;

  /**  Create a map of signal shapes for the ECAL detector.
   *   @return status code
   */
  std::vector<float> CreateSignalShape_BasicDirac(int nSamples, int DiracIdx, float DiracAmplitude);

private:
  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;

  /// Names of the detector readout for the volumes
  Gaudi::Property<std::vector<std::string>> m_readoutNamesSegmented{this, "readoutNames", {"ECalBarrelModuleThetaMerged"}};
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

  /// Name of output file
  std::string m_outputFileName;

  Gaudi::Property<int> m_len {this, "pulseSampling", 10, "Number of samples in pulse" };
  Gaudi::Property<int> m_idx {this, "pulseIndex", 2, "Index of Dirac pulse" };
  Gaudi::Property<float> m_amp {this, "pulseAmplitude", 2.0, "Amplitude of Dirac pulse" };
};

#endif /* RECFCCEECALORIMETER_CREATEFCCEECALOSIGNALSHAPES_H */
