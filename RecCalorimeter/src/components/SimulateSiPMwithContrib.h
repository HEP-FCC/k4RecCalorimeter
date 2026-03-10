#ifndef SimulateSiPMwithContrib_h
#define SimulateSiPMwithContrib_h 1

#include "edm4hep/CaloHitSimCaloHitLinkCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/TimeSeriesCollection.h"

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

/** @class SimulateSiPMwithContrib
 *
 *  This algorithm was created to apply the SiPM response over the
 *  hit created with IDEA_o2 simulation.
 *  Contributions from such hits containt the photo-electrons (fired cells)
 *  and their time of arrivals at SiPMs.
 *  This algorithm takes all the input collections form the IDEA_o2
 *  hadronic calorimeter and returns corresponding digitized hits.
 *
 *  The algorithm simulates SiPM response including effects like:
 *  - Dark count rate (DCR)
 *  - Crosstalk (Xt)
 *  - Afterpulsing
 *  - Cell recovery time
 *  - Signal rise and fall times
 *
 *  @author Lorenzo Pezzotti
 *  @date   2026-02-17
 */

class SimulateSiPMwithContrib : public Gaudi::Algorithm {
public:
  SimulateSiPMwithContrib(const std::string& name, ISvcLocator* svcLoc);
  virtual ~SimulateSiPMwithContrib() {};

  StatusCode initialize();
  StatusCode execute(const EventContext&) const;
  StatusCode finalize();

private:
  // Random Number Service
  SmartIF<IRndmGenSvc> m_randSvc;
  Rndm::Numbers m_rndmExp;

  // readout name and segmentation (of specific type)
  // Not used for the moment, might be used in future
  // Gaudi::Property<std::string> m_readoutName{this, "readoutName", "", "name of the readout"};

  mutable k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection> m_simHits{"", Gaudi::DataHandle::Reader, this};
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_digiHits{"", Gaudi::DataHandle::Writer, this};
  mutable k4FWCore::DataHandle<edm4hep::CaloHitSimCaloHitLinkCollection> m_hitLinks{"", Gaudi::DataHandle::Writer,
                                                                                    this};

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
  Gaudi::Property<double> m_recovery{this, "recovery", 10.,
                                     "SiPM cell recovery time in ns"}; // https://arxiv.org/abs/2001.10322

  // Noise parameters
  Gaudi::Property<double> m_Dcr{this, "DCR", 120e3, "SiPM DCR"}; // dark count rate
  Gaudi::Property<double> m_Xt{this, "Xtalk", 0.01, "SiPM crosstalk"};
  Gaudi::Property<double> m_afterpulse{this, "afterpulse", 0.03, "afterpulse probability"};
  Gaudi::Property<double> m_snr{this, "SNR", 20., "SNR value in dB"}; // signal-to-noise ratio

  // integration parameters
  Gaudi::Property<double> m_gateStart{this, "gateStart", 5., "Integration gate starting time in ns"};
  Gaudi::Property<double> m_gateL{this, "gateLength", 95.,
                                  "Integration gate length in ns"}; // Should be approx 5 times fallTimeFast (see above)
  Gaudi::Property<double> m_thres{this, "threshold", 1.5,
                                  "Integration threshold"}; // Threshold in p.e. (1.5 to suppress DCR)

  // other parameters (attention, will override above parameters if set)
  Gaudi::Property<std::map<std::string, double>> m_params{this, "params", {}, "optional parameters"};

  // switch to select if you are processing scintillation or Cherenkov hits
  Gaudi::Property<bool> m_isCherenkov{this, "isCherenkov", false,
                                      "switch to select if you are processing scintillation or Cherenkov hits"};
  // SiPM efficiency
  Gaudi::Property<std::vector<double>> m_wavelen{
      this, "wavelength", {1000., 100.}, "wavelength vector in nm (decreasing order)"};
  Gaudi::Property<std::vector<double>> m_sipmEff{this, "sipmEfficiency", {1.0, 1.0}, "SiPM efficiency vs wavelength"};

  Gaudi::Property<double> m_scintDecaytime{this, "scintDecaytime", 2.8, "scintillation decay time in ns"};

  // option to store full waveform (for debugging)
  // Not used for the moment, might be used in future for debugging
  // Gaudi::Property<bool> m_storeFullWaveform{this, "storeFullWaveform", false, "Store full waveform for debugging"};
};

#endif
