#ifndef RECCALORIMETER_CALIBRATEBENCHMARKMETHOD_H
#define RECCALORIMETER_CALIBRATEBENCHMARKMETHOD_H

// GAUDI
#include "GaudiAlg/GaudiAlgorithm.h"

#include "TObject.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <iostream>
#include "Math/IFunction.h"

// Key4HEP
#include "k4FWCore/DataHandle.h"
class IGeoSvc;

// EDM4HEP
namespace edm4hep {
  class CalorimeterHitCollection;
  class MCParticleCollection;
}

namespace DD4hep {
namespace DDSegmentation {
class Segmentation;
}
}

class TH1F;
class ITHistSvc;

/** @class CalibrateBenchmarkMethod CalibrateBenchmarkMethod.h
 *
 *  Sums energy deposited in every calorimeter layer separately, sums also energy deposited in the dead material of the
 *  calorimeter (cryostat).
 *
 *  In order to calculate upstream or downstream energy correction one needs to mark cryostat as sensitive in the
 *  detector XML files. Additionally, for the downstream correction, the thickness of the back cryostat needs to be
 *  enlarged to be unrealistically large (at least one meter) to capture all energy deposited behind the calorimeter.
 *
 *  Based on work done by Anna Zaborowska, Jana Faltova and Juraj Smiesko
 *
 *  @author Michaela Mlynarikova
 */

class CalibrateBenchmarkMethod : public GaudiAlgorithm {
public:
  explicit CalibrateBenchmarkMethod(const std::string&, ISvcLocator*);
  virtual ~CalibrateBenchmarkMethod();
  /**  Initialize.
   *   @return status code
   */
  virtual StatusCode initialize() final;
  /**  Fills the histograms.
   *   @return status code
   */
  virtual StatusCode execute() final;
  /**  Finalize.
   *   @return status code
   */
  virtual StatusCode finalize() final;

private:
  /// Pointer to the interface of histogram service
  ServiceHandle<ITHistSvc> m_histSvc;

  /// Handle for electromagnetic barrel cells (input collection)
  DataHandle<edm4hep::CalorimeterHitCollection> m_ecalBarrelCells{"ecalBarrelCells", Gaudi::DataHandle::Reader, this};
  /// Handle for hadronic barrel cells (input collection)
  DataHandle<edm4hep::CalorimeterHitCollection> m_hcalBarrelCells{"hcalBarrelCells", Gaudi::DataHandle::Reader, this};

  // Histogram of total deposited energy in the calorimeters
  TH1F* m_totalEnergyECal;
  TH1F* m_totalEnergyHCal;
  TH1F* m_totalEnergyBoth;
  TH1F* h_parameter; 

  // Maximum energy for the x-axis range
  Gaudi::Property<double> m_energy{this, "energy", 100, "Generated energy"};

  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;

  /// Number of ECal and HCal layers 
  Gaudi::Property<size_t> m_numLayersECal{this, "numLayersECal", 12, "Number of ECal layers"};
  Gaudi::Property<size_t> m_numLayersHCal{this, "numLayersHCal", 13, "Number of HCal layers"};


  /// Id of the first HCal layer
  Gaudi::Property<uint> m_firstLayerHCal{this, "firstLayerHCal", 0, "ID of first HCal layer"};

  /// Names of the detector readouts, corresponding to system IDs
  Gaudi::Property<std::vector<std::string>> m_readoutNames {this, "readoutNames", {"ECalBarrelReadout","HCalBarrelReadout"},"Names of the detector readout, corresponding to systemId"};
  Gaudi::Property<std::vector<uint>> m_SystemID {this, "systemID", {4,8},"systemId"};
  Gaudi::Property<uint> m_ECalSystemID{this, "ECalSystemID", 4, "ID of ECal system"};
  Gaudi::Property<uint> m_HCalSystemID{this, "HCalSystemID", 8, "ID of the HCal system"};

};
#endif /* RECCALORIMETER_CALIBRATEBENCHMARKMETHOD_H */