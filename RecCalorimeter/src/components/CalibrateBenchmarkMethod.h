#ifndef RECCALORIMETER_CALIBRATEBENCHMARKMETHOD_H
#define RECCALORIMETER_CALIBRATEBENCHMARKMETHOD_H

// GAUDI
#include "GaudiAlg/GaudiAlgorithm.h"

// Key4HEP
#include "k4FWCore/DataHandle.h"
class IGeoSvc;
class TH1F;

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
 * Benchmark calibration should be used for the combined simulation of ECal and HCal when using charged pions. 
 * As an input it expect ECal to be calibrated to EM scale and HCal to be calibrated to HAD scale. 
 * The aim of the benchmark calibration is to bring ECal to HAD scale and also to take into account
 * the energy loss between the ECal and HCal (e.g. in cryostat) - for this, the energy from the last ECal layer and the first HCal layer is used. 
 * The output parameters from the fit are stored in a histogram m_parameters and these parameters are supposed to be used by CorrectClustersBenchmarkMethod 
 * to apply the benchmark method output to clusters and correct output cluster energy.
 * To obtain the actual parameters run RecCalorimeter/tests/options/fcc_ee_caloBenchmarkCalibration.py
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
  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;

  /// Pointer to the interface of histogram service
  ServiceHandle<ITHistSvc> m_histSvc;

  double chiSquareFitBarrel(const Double_t *par) const;

  /// Handle for electromagnetic barrel cells (input collection)
  DataHandle<edm4hep::CalorimeterHitCollection> m_ecalBarrelCells{"ecalBarrelCells", Gaudi::DataHandle::Reader, this};
  /// Handle for hadronic barrel cells (input collection)
  DataHandle<edm4hep::CalorimeterHitCollection> m_hcalBarrelCells{"hcalBarrelCells", Gaudi::DataHandle::Reader, this};

  /// Histogram of total deposited energy in the calorimeters 
  TH1F* m_totalEnergyECal;
  TH1F* m_totalEnergyHCal;
  TH1F* m_totalEnergyBoth;

  /// An output histogram, the values of fit parameters are stored in the bins
  TH1F* m_parameters; 

  /// vectors to store the energy in each ECal/HCal layer 
  std::vector<double> m_energyInECalLayer;
  std::vector<double> m_energyInHCalLayer;


  /// Maximum energy for the x-axis range
  Gaudi::Property<double> m_energy{this, "energy", 100, "Generated energy"};

  /// Number of ECal and HCal layers 
  Gaudi::Property<size_t> m_numLayersECal{this, "numLayersECal", 12, "Number of ECal layers"};
  Gaudi::Property<size_t> m_numLayersHCal{this, "numLayersHCal", 13, "Number of HCal layers"};

  /// ID of the first HCal layer
  Gaudi::Property<uint> m_firstLayerHCal{this, "firstLayerHCal", 0, "ID of first HCal layer"};

  /// Names of the detector readouts, corresponding to system IDs
  Gaudi::Property<std::vector<std::string>> m_readoutNames {this, "readoutNames", {"ECalBarrelReadout","HCalBarrelReadout"},"Names of the detector readout, corresponding to systemId"};
  Gaudi::Property<std::vector<uint>> m_SystemID {this, "systemID", {4,8},"systemId"};
  Gaudi::Property<uint> m_ECalSystemID{this, "ECalSystemID", 4, "ID of ECal system"};
  Gaudi::Property<uint> m_HCalSystemID{this, "HCalSystemID", 8, "ID of the HCal system"};

  /// vectors containing the energy deposits to be used for minimization
  std::vector<double> m_vecEgenerated;
  std::vector<double> m_vecEinECaltotal;
  std::vector<double> m_vecEinHCaltotal;
  std::vector<double> m_vecEinHCalfirstLayer;
  std::vector<double> m_vecEinECallastLayer;

};
#endif /* RECCALORIMETER_CALIBRATEBENCHMARKMETHOD_H */