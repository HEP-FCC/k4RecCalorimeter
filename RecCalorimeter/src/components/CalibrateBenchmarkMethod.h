#ifndef RECCALORIMETER_CALIBRATEBENCHMARKMETHOD_H
#define RECCALORIMETER_CALIBRATEBENCHMARKMETHOD_H

// GAUDI
#include "Gaudi/Algorithm.h"

// Key4HEP
#include "k4FWCore/DataHandle.h"
class IGeoSvc;
class TH1F;

// EDM4HEP
namespace edm4hep {
class CalorimeterHitCollection;
class MCParticleCollection;
} // namespace edm4hep

namespace DD4hep {
namespace DDSegmentation {
  class Segmentation;
}
} // namespace DD4hep

class TH1F;
class ITHistSvc;

/** @class CalibrateBenchmarkMethod CalibrateBenchmarkMethod.h
 *
 * Benchmark calibration should be used for the combined simulation of ECal and HCal when using charged pions.
 * As an input it expect ECal to be calibrated to EM scale and HCal to be calibrated to HAD scale.
 * The aim of the benchmark calibration is to bring ECal to HAD scale and also to take into account
 * the energy loss between the ECal and HCal (e.g. in cryostat) - for this, the energy from the last ECal layer and the
 * first HCal layer is used. The output parameters from the fit are stored in a histogram m_parameters and these
 * parameters are supposed to be used by CorrectClustersBenchmarkMethod to apply the benchmark method output to clusters
 * and correct output cluster energy. To obtain the actual parameters run
 * RecCalorimeter/tests/options/fcc_ee_caloBenchmarkCalibration.py
 *
 *  Based on work done by Anna Zaborowska, Jana Faltova and Juraj Smiesko
 *
 *  @author Michaela Mlynarikova
 */

class CalibrateBenchmarkMethod : public Gaudi::Algorithm {
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
  virtual StatusCode execute(const EventContext&) const final;
  /**  Finalize.
   *   @return status code
   */
  virtual StatusCode finalize() final;

private:
  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;

  /// Pointer to the interface of histogram service
  ServiceHandle<ITHistSvc> m_histSvc;

  double chiSquareFitBarrel(const double* par) const;
  void registerHistogram(const std::string& path, TH1F*& histogramName);
  void runMinimization(int n_param, const std::vector<double>& variable, const std::vector<double>& steps,
                       const std::vector<int>& fixedParameters) const;

  /// Handle for electromagnetic barrel cells (input collection)
  mutable DataHandle<edm4hep::CalorimeterHitCollection> m_ecalBarrelCells{"ecalBarrelCells", Gaudi::DataHandle::Reader,
                                                                          this};
  /// Handle for hadronic barrel cells (input collection)
  mutable DataHandle<edm4hep::CalorimeterHitCollection> m_hcalBarrelCells{"hcalBarrelCells", Gaudi::DataHandle::Reader,
                                                                          this};

  /// Histogram of total deposited energy in the calorimeters
  TH1F* m_totalEnergyECal;
  TH1F* m_totalEnergyHCal;
  TH1F* m_totalEnergyBoth;

  /// An output histogram, the values of fit parameters are stored in the bins
  TH1F* m_parameters;

  /// vectors to store the energy in each ECal/HCal layer
  mutable std::vector<double> m_energyInLayerECal;
  mutable std::vector<double> m_energyInLayerHCal;

  /// Maximum energy for the x-axis range
  Gaudi::Property<double> m_energy{this, "energy", 100, "Generated energy"};

  /// Number of ECal and HCal layers
  Gaudi::Property<size_t> m_numLayersECal{this, "numLayersECal", 12, "Number of ECal layers"};
  Gaudi::Property<size_t> m_numLayersHCal{this, "numLayersHCal", 13, "Number of HCal layers"};

  /// ID of the first ECal/HCal layer
  Gaudi::Property<uint> m_firstLayerECal{this, "firstLayerECal", 0, "ID of first ECal layer"};
  Gaudi::Property<uint> m_firstLayerHCal{this, "firstLayerHCal", 0, "ID of first HCal layer"};

  /// Names of the detector readouts, corresponding to system IDs
  Gaudi::Property<std::vector<std::string>> m_readoutNames{this,
                                                           "readoutNames",
                                                           {"ECalBarrelReadout", "HCalBarrelReadout"},
                                                           "Names of the detector readout, corresponding to systemId"};
  Gaudi::Property<std::vector<uint>> m_systemID{this, "systemID", {4, 8}, "systemId"};
  Gaudi::Property<uint> m_systemIDECal{this, "ECalSystemID", 4, "ID of ECal system"};
  Gaudi::Property<uint> m_systemIDHCal{this, "HCalSystemID", 8, "ID of the HCal system"};

  /// vectors containing the energy deposits to be used for minimization
  mutable std::vector<double> m_vecGeneratedEnergy;
  mutable std::vector<double> m_vecTotalEnergyinECal;
  mutable std::vector<double> m_vecTotalEnergyinHCal;
  mutable std::vector<double> m_vecEnergyInFirstLayerECal;
  mutable std::vector<double> m_vecEnergyInLastLayerECal;
  mutable std::vector<double> m_vecEnergyInFirstLayerHCal;

  // benchmark parameters which should be fixed
  // p[1] because HCal is already calibrated to HAD scale
  // p[5] is the constant term and was not bringing large improvement, hence not minimized for the moment (might be
  // tried again in the future)
  Gaudi::Property<std::vector<int>> m_fixedParameters = {
      this, "fixedParameters", {1, 5}, "Fixed parameters that will not be minimized"};
};
#endif /* RECCALORIMETER_CALIBRATEBENCHMARKMETHOD_H */
