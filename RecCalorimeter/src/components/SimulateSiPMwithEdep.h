#ifndef SimulateSiPMwithEdep_h
#define SimulateSiPMwithEdep_h 1

#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/EventHeaderCollection.h"

#include "k4FWCore/DataHandle.h"

#include "Gaudi/Algorithm.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/ToolHandle.h"

#include "sipm/SiPMSensor.h"

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
  Gaudi::Property<std::string> m_outColl{this, "outputHitCollection", "DigiCalorimeterHits_scint", "output calo collection name"};

  mutable DataHandle<edm4hep::SimCalorimeterHitCollection> m_scintHits{m_hitColl, Gaudi::DataHandle::Reader, this};
  mutable DataHandle<edm4hep::CalorimeterHitCollection> m_digiHits{m_outColl, Gaudi::DataHandle::Writer, this};

  std::unique_ptr<sipm::SiPMSensor> m_sensor;

  // Hamamatsu S14160-1310PS
  Gaudi::Property<double> m_sigLength{this, "signalLength", 200., "signal length in ns"};
  Gaudi::Property<double> m_sipmSize{this, "SiPMsize", 1.3, "width of photosensitive area in mm"};
  Gaudi::Property<double> m_Dcr{this, "DCR", 120e3, "SiPM DCR"}; // dark count rate
  Gaudi::Property<double> m_Xt{this, "Xtalk", 0.01, "SiPM crosstalk"};
  Gaudi::Property<double> m_sampling{this, "sampling", 0.1, "SiPM sampling rate in ns"};
  Gaudi::Property<double> m_recovery{this, "recovery", 10., "SiPM cell recovery time in ns"}; // https://arxiv.org/abs/2001.10322
  Gaudi::Property<double> m_cellPitch{this, "cellpitch", 10., "SiPM cell size in um"};
  Gaudi::Property<double> m_afterpulse{this, "afterpulse", 0.03, "afterpulse probability"};
  Gaudi::Property<double> m_falltimeFast{this, "falltimeFast", 6.5, "signal fast component decay time in ns"};
  Gaudi::Property<double> m_risetime{this, "risetime", 1., "signal rise time in ns"};
  Gaudi::Property<double> m_snr{this, "SNR", 20., "SNR value in dB"}; // signal-to-noise ratio

  // integration parameters
  Gaudi::Property<double> m_gateStart{this, "gateStart", 5., "Integration gate starting time in ns"};
  Gaudi::Property<double> m_gateL{this, "gateLength", 95., "Integration gate length in ns"};  // Should be approx 5 times fallTimeFast (see above)
  Gaudi::Property<double> m_thres{this, "threshold", 1.5, "Integration threshold"};  // Threshold in p.e. (1.5 to suppress DCR)

  Gaudi::Property<std::vector<double>> m_wavelen{this, "wavelength", {1000., 100.}, "wavelength vector in nm"};
  Gaudi::Property<std::vector<double>> m_sipmEff{this, "sipmEfficiency", {0.1, 0.1}, "SiPM efficiency vs wavelength"};
  Gaudi::Property<std::vector<double>> m_scintSpectrum{this, "scintSpectrum", {1., 1.}, "scintillation spectrum vs wavelength"};
  Gaudi::Property<std::vector<double>> m_filterEff{this, "filterEfficiency", {0.1, 0.1}, "optical filter efficiency vs wavelength"};
  Gaudi::Property<double> m_scintDecaytime{this, "scintDecaytime", 2.8, "scintillation decay time in ns"};

  // TB-based empirical values after considering attenuation, NIM time window, etc.
  Gaudi::Property<double> m_scintYield{this, "scintYield", 2.5, "scintillation yield in /keV"};
  Gaudi::Property<double> m_timeInterval{this, "timeInterval", 15., "time to the first signal in ns"};

  // integral table - initialized at SimulateSiPMwithEdep::initialize()
  std::vector<double> m_integral;
  double m_efficiency;
};

#endif
