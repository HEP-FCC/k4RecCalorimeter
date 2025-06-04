#ifndef RECCALORIMETER_SimulateSiPMwithOpticalPhoton_H
#define RECCALORIMETER_SimulateSiPMwithOpticalPhoton_H

// EDM4HEP includes
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/RawTimeSeriesCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/TimeSeriesCollection.h"

// k4FWCore includes
#include "k4FWCore/DataHandle.h"

// Gaudi includes
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

/** @class SimulateSiPMwithOpticalPhoton
 *
 *  Algorithm for digitizing the SiPM response
 *  This digitizer takes RawTimeSeriesCollection and RawCalorimeterHitCollection as input
 *  and produces TimeSeriesCollection and CalorimeterHitCollection as output.
 *
 *  The algorithm simulates SiPM response including effects like:
 *  - Dark count rate (DCR)
 *  - Crosstalk (Xt)
 *  - Afterpulsing
 *  - Cell recovery time
 *  - Signal rise and fall times
 *
 *  @author Sanghyun Ko
 *  @author Sungwon Kim
 *  @date   2025-03-18
 */

class SimulateSiPMwithOpticalPhoton : public Gaudi::Algorithm {
public:
  SimulateSiPMwithOpticalPhoton(const std::string& name, ISvcLocator* svcLoc);
  virtual ~SimulateSiPMwithOpticalPhoton() {};

  StatusCode initialize() override;
  StatusCode execute(const EventContext&) const override;
  StatusCode finalize() override;

private:
  // integral function for the wavelength random generation
  std::vector<double> integral(const std::vector<double>& wavelen, const std::vector<double>& yval) const;

  // Random Number Service
  SmartIF<IRndmGenSvc> m_randSvc;
  Rndm::Numbers m_rndmUniform;

  // input collection names
  Gaudi::Property<std::string> m_hitColl{this, "inputHitCollection", "DRcaloSiPMreadoutSimHit",
                                         "input calo collection name"};
  Gaudi::Property<std::string> m_outColl{this, "outputHitCollection", "DRcaloSiPMreadoutDigiHit",
                                         "output calo collection name"};

  Gaudi::Property<std::string> m_inTimeColl{this, "inputTimeStructCollection", "DRcaloSiPMreadoutTimeStruct",
                                            "input time structure collection name"};
  Gaudi::Property<std::string> m_inWavlenColl{this, "inputWavlenCollection", "DRcaloSiPMreadoutWaveLen",
                                              "input wavelength collection name"};
  Gaudi::Property<std::string> m_outTimeColl{this, "outputTimeStructCollection", "DRcaloSiPMreadoutDigiWaveform",
                                             "output waveform collection name"};

  // Input collections
  mutable k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection> m_simHits{m_hitColl, Gaudi::DataHandle::Reader,
                                                                               this};
  mutable k4FWCore::DataHandle<edm4hep::RawTimeSeriesCollection> m_timeStruct{m_inTimeColl, Gaudi::DataHandle::Reader,
                                                                              this};
  mutable k4FWCore::DataHandle<edm4hep::RawTimeSeriesCollection> m_wavelenStruct{m_inWavlenColl,
                                                                                 Gaudi::DataHandle::Reader, this};

  // Output collections
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_digiHits{m_outColl, Gaudi::DataHandle::Writer,
                                                                             this};
  mutable k4FWCore::DataHandle<edm4hep::TimeSeriesCollection> m_waveforms{m_outTimeColl, Gaudi::DataHandle::Writer,
                                                                          this};

  // SiPM sensor model
  std::unique_ptr<sipm::SiPMSensor> m_sensor;

  // SiPM properties (defaults set based on Hamamatsu S14160-1310PS)
  // Signal properties
  Gaudi::Property<double> m_sigLength{this, "signalLength", 200., "Signal length in ns"};
  Gaudi::Property<double> m_sampling{this, "sampling", 0.1, "SiPM sampling rate in ns"};
  Gaudi::Property<double> m_risetime{this, "risetime", 1., "Signal rise time in ns"};
  Gaudi::Property<double> m_falltimeFast{this, "falltimeFast", 6.5, "Signal fast component decay time in ns"};

  // SiPM physical properties
  Gaudi::Property<double> m_sipmSize{this, "SiPMsize", 1.3, "Width of photosensitive area in mm"};
  Gaudi::Property<double> m_cellPitch{this, "cellpitch", 10., "SiPM cell size in um"};
  Gaudi::Property<double> m_recovery{this, "recovery", 10., "SiPM cell recovery time in ns"};

  // Noise parameters
  Gaudi::Property<double> m_Dcr{this, "DCR", 120e3, "SiPM dark count rate in Hz"};
  Gaudi::Property<double> m_Xt{this, "Xtalk", 0.01, "SiPM optical crosstalk probability"};
  Gaudi::Property<double> m_afterpulse{this, "afterpulse", 0.03, "Afterpulse probability"};
  Gaudi::Property<double> m_snr{this, "SNR", 20., "Signal-to-noise ratio in dB"};

  // Integration parameters
  Gaudi::Property<double> m_gateStart{this, "gateStart", 5., "Integration gate starting time in ns"};
  Gaudi::Property<double> m_gateL{this, "gateLength", 95., "Integration gate length in ns"};
  Gaudi::Property<double> m_thres{this, "threshold", 1.5, "Integration threshold in photoelectrons"};

  // SiPM efficiency
  Gaudi::Property<std::vector<double>> m_wavelen{
      this, "wavelength", {1000., 100.}, "wavelength vector in nm (decreasing order)"};
  Gaudi::Property<std::vector<double>> m_sipmEff{this, "sipmEfficiency", {0.1, 0.1}, "SiPM efficiency vs wavelength"};

  // scale ADC to energy
  Gaudi::Property<double> m_scaleADC{this, "scaleADC", 1., "calibration factor for scaling ADC to energy"};
};

#endif // RECCALORIMETER_SimulateSiPMwithOpticalPhoton_H
