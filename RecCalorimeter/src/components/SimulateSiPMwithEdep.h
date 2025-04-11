#ifndef SimulateSiPMwithEdep_h
#define SimulateSiPMwithEdep_h 1

#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/TimeSeriesCollection.h"
#include "edm4hep/EventHeaderCollection.h"

#include "k4FWCore/DataHandle.h"

#include "Gaudi/Algorithm.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/ToolHandle.h"

// Check for SiPMSensor header location (similar to how DigiSiPM handles this)
#if __has_include("SiPMSensor.h")
#include "SiPMSensor.h"
#elif __has_include("sipm/SiPMSensor.h")
#include "sipm/SiPMSensor.h"
#endif

/** @class SimulateSiPMwithEdep
 *
 *  Algorithm for digitizing the SiPM response using SIM energy deposits,
 *  taking SimCalorimeterHitCollection (energy deposited in scintillators) as input
 *  and produces CalorimeterHitCollection (ADC count from SiPM) as output.
 *
 *  The algorithm simulates SiPM response including effects like:
 *  - Dark count rate (DCR)
 *  - Crosstalk (Xt)
 *  - Afterpulsing
 *  - Cell recovery time
 *  - Signal rise and fall times
 *
 *  @author Sanghyun Ko
 *  @date   2025-03-27
 */

class SimulateSiPMwithEdep : public Gaudi::Algorithm {
public:
  SimulateSiPMwithEdep(const std::string& name, ISvcLocator* svcLoc);
  virtual ~SimulateSiPMwithEdep() {};

  StatusCode initialize();
  StatusCode execute(const EventContext&) const;
  StatusCode finalize();

private:
  std::vector<double> integral(const std::vector<double>& wavelen, const std::vector<double>& yval) const;

  // Random Number Service
  SmartIF<IRndmGenSvc> m_randSvc;
  Rndm::Numbers m_rndmUniform;
  Rndm::Numbers m_rndmExp;

  // input collection names
  Gaudi::Property<std::string> m_hitColl{this, "inputHitCollection", "DRcaloSiPMreadout_scint", "input calo collection name"};
  Gaudi::Property<std::string> m_outColl{this, "outputHitCollection", "DRcaloSiPMreadoutDigiHit_scint", "output calo collection name"};
  Gaudi::Property<std::string> m_outTimeColl{this, "outputTimeStructCollection", "DRcaloSiPMreadoutDigiWaveform_scint", "output waveform collection name"};

  mutable DataHandle<edm4hep::SimCalorimeterHitCollection> m_scintHits{m_hitColl, Gaudi::DataHandle::Reader, this};
  mutable DataHandle<edm4hep::CalorimeterHitCollection> m_digiHits{m_outColl, Gaudi::DataHandle::Writer, this};
  mutable DataHandle<edm4hep::TimeSeriesCollection> m_waveforms{m_outTimeColl, Gaudi::DataHandle::Writer, this};

  std::unique_ptr<sipm::SiPMSensor> m_sensor;

  // Hamamatsu S14160-1310PS
  // Signal properties
  Gaudi::Property<double> m_sigLength{this, "signalLength", 200., "signal length in ns"};
  Gaudi::Property<double> m_sampling{this, "sampling", 0.1, "SiPM sampling rate in ns"};
  Gaudi::Property<double> m_risetime{this, "risetime", 1., "signal rise time in ns"};
  Gaudi::Property<double> m_falltimeFast{this, "falltimeFast", 6.5, "signal fast component decay time in ns"};

  // SiPM physical properties
  Gaudi::Property<double> m_sipmSize{this, "SiPMsize", 1.3, "width of photosensitive area in mm"};
  Gaudi::Property<double> m_cellPitch{this, "cellpitch", 10., "SiPM cell size in um"};
  Gaudi::Property<double> m_recovery{this, "recovery", 10., "SiPM cell recovery time in ns"}; // https://arxiv.org/abs/2001.10322

  // Noise parameters
  Gaudi::Property<double> m_Dcr{this, "DCR", 120e3, "SiPM DCR"}; // dark count rate
  Gaudi::Property<double> m_Xt{this, "Xtalk", 0.01, "SiPM crosstalk"};
  Gaudi::Property<double> m_afterpulse{this, "afterpulse", 0.03, "afterpulse probability"};
  Gaudi::Property<double> m_snr{this, "SNR", 20., "SNR value in dB"}; // signal-to-noise ratio

  // integration parameters
  Gaudi::Property<double> m_gateStart{this, "gateStart", 5., "Integration gate starting time in ns"};
  Gaudi::Property<double> m_gateL{this, "gateLength", 95., "Integration gate length in ns"};  // Should be approx 5 times fallTimeFast (see above)
  Gaudi::Property<double> m_thres{this, "threshold", 1.5, "Integration threshold"};  // Threshold in p.e. (1.5 to suppress DCR)

  // SiPM efficiency, filter efficiency
  Gaudi::Property<std::vector<double>> m_wavelen{this, "wavelength", {1000., 100.}, "wavelength vector in nm (decreasing order)"};
  Gaudi::Property<std::vector<double>> m_sipmEff{this, "sipmEfficiency", {0.1, 0.1}, "SiPM efficiency vs wavelength"};
  Gaudi::Property<std::vector<double>> m_filterEff{this, "filterEfficiency", {0.1, 0.1}, "optical filter efficiency vs wavelength"};

  // scintillator optical properties
  Gaudi::Property<double> m_refractiveIndex{this, "refractiveIndex", 1.6, "scintillator refractive index"};
  Gaudi::Property<std::vector<double>> m_absLen{this, "absorptionLength", {9999., 9999.}, "scintillator absorption length in meter"};

  // Scintillation spectrum and decay time
  Gaudi::Property<std::vector<double>> m_scintSpectrum{this, "scintSpectrum", {1., 1.}, "scintillation spectrum vs wavelength"};
  Gaudi::Property<double> m_scintDecaytime{this, "scintDecaytime", 2.8, "scintillation decay time in ns"};

  // TB-based empirical values after considering attenuation, NIM time window, etc.
  // set yield to 10^-6 (i.e. keV/GeV) if the stored Edep is equal to the number of photon
  Gaudi::Property<double> m_scintYield{this, "scintYield", 2.5, "scintillation yield in /keV"};

  // scale ADC to energy
  Gaudi::Property<double> m_scaleADC{this, "scaleADC", 1., "calibration factor for scaling ADC to energy"};

  // integral table - initialized at SimulateSiPMwithEdep::initialize()
  std::vector<double> m_integral;
  double m_efficiency;
};

#endif
