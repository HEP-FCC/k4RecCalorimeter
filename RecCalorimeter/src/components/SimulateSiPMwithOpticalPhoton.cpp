#include "SimulateSiPMwithOpticalPhoton.h"

#include <cmath>
#include <stdexcept>

// DECLARE_COMPONENT macro connects the algorithm to the framework
DECLARE_COMPONENT(SimulateSiPMwithOpticalPhoton)

SimulateSiPMwithOpticalPhoton::SimulateSiPMwithOpticalPhoton(const std::string& aName, ISvcLocator* aSvcLoc)
  : Gaudi::Algorithm(aName, aSvcLoc) {
  // Constructor doesn't need to do anything special
}

StatusCode SimulateSiPMwithOpticalPhoton::initialize() {
  // First call the base class initialization
  StatusCode sc = Gaudi::Algorithm::initialize();
  if (sc.isFailure()) return sc;

  if (m_wavelen.size() < 2) {
    error() << "SimulateSiPMwithOpticalPhoton: "
            << "The wavelength vector size must be greater or equal than 2" << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_wavelen.size()!=m_sipmEff.size()) {
    error() << "SimulateSiPMwithOpticalPhoton: "
            << "The SiPM efficiency vector size should be equal to the wavelength vector size" << endmsg;
    return StatusCode::FAILURE;
  }

  // Create and initialize SiPM sensor with the properties from job options
  sipm::SiPMProperties properties;
  properties.setSignalLength(m_sigLength);
  properties.setSize(m_sipmSize);
  properties.setDcr(m_Dcr);
  properties.setXt(m_Xt);
  properties.setSampling(m_sampling);
  properties.setRecoveryTime(m_recovery);
  properties.setPitch(m_cellPitch);
  properties.setAp(m_afterpulse);
  properties.setFallTimeFast(m_falltimeFast);
  properties.setRiseTime(m_risetime);
  properties.setSnr(m_snr);
  // Set the PDE type to spectrum PDE
  properties.setPdeType(sipm::SiPMProperties::PdeType::kSpectrumPde);
  // Set the PDE spectrum
  properties.setPdeSpectrum(m_wavelen,m_sipmEff);

  // Create the SiPM sensor model
  m_sensor = std::make_unique<sipm::SiPMSensor>(properties);

  info() << "SimulateSiPMwithOpticalPhoton initialized" << endmsg;
  info() << properties << endmsg; // sipm::SiPMProperties has std::ostream& operator<<

  return StatusCode::SUCCESS;
}

StatusCode SimulateSiPMwithOpticalPhoton::execute(const EventContext&) const {
  // Get input collections
  const edm4hep::RawTimeSeriesCollection* timeStructs = m_timeStruct.get();
  const edm4hep::RawTimeSeriesCollection* waveLenStructs = m_waveLen.get();
  const edm4hep::RawCalorimeterHitCollection* rawHits = m_rawHits.get();

  // Create output collections
  edm4hep::TimeSeriesCollection* waveforms = m_waveforms.createAndPut();
  edm4hep::CalorimeterHitCollection* digiHits = m_digiHits.createAndPut();

  // Process each hit
  for (unsigned int idx = 0; idx < timeStructs->size(); idx++) {
    const auto& timeStruct = timeStructs->at(idx);
    const auto& waveLength = waveLenStructs->at(idx);
    const auto& rawhit = rawHits->at(idx);

    // Validate that both collections have matching CellIDs
    if (timeStruct.getCellID() != rawhit.getCellID()) {
      error() << "CellIDs of RawCalorimeterHits & RawTimeSeries are different! "
              << "This should never happen." << endmsg;
      return StatusCode::FAILURE;
    }

    // Extract photon arrival times from the time structure
    std::vector<double> vecTimes;
    std::vector<double> vecWavelengths;
    vecTimes.reserve(rawhit.getAmplitude());
    vecWavelengths.reserve(rawhit.getAmplitude());

    for (unsigned int bin = 0; bin < timeStruct.adcCounts_size(); bin++) {
      int counts = static_cast<int>(timeStruct.getAdcCounts(bin));
      double timeBin = timeStruct.getTime() + timeStruct.getInterval() * (static_cast<float>(bin) + 0.5);

      // Add a photon arrival time for each count in this bin
      for (int num = 0; num < counts; num++) {
        vecTimes.emplace_back(timeBin);
      }
    }

    // Extract photon wavelengths from the wave length structure
    for (unsigned int bin = 0; bin < waveLength.adcCounts_size(); bin++) {
      int counts = static_cast<int>(waveLength.getAdcCounts(bin));
      double waveLenBin = waveLength.getTime() + waveLength.getInterval() * (static_cast<float>(bin) + 0.5);

      // Add a photon arrival time for each count in this bin
      for (int num = 0; num < counts; num++) {
        vecWavelengths.emplace_back(waveLenBin);
      }
    }

    // Reset the SiPM state and run the simulation
    m_sensor->resetState();
    m_sensor->addPhotons(vecTimes, vecWavelengths); // Sets photon arrival times (in ns) & wavelengths (in nm)
    m_sensor->runEvent();        // Runs the SiPM simulation

    auto digiHit = digiHits->create();
    auto waveform = waveforms->create();

    // Get the analog signal from the SiPM
    const sipm::SiPMAnalogSignal anaSignal = m_sensor->signal();

    // Compute integral (energy) and time of arrival
    const double integral = anaSignal.integral(m_gateStart, m_gateL, m_thres); // (intStart, intGate, threshold)
    const double toa = anaSignal.toa(m_gateStart, m_gateL, m_thres);           // (intStart, intGate, threshold)
    const double gateEnd = m_gateStart.value() + m_gateL.value();

    // Set digitized hit properties
    digiHit.setEnergy(integral);
    digiHit.setCellID(rawhit.getCellID());
    digiHit.setTime(toa + m_gateStart); // Toa and m_gateStart are in ns

    // Set waveform properties
    waveform.setInterval(m_sampling);
    waveform.setTime(timeStruct.getTime());
    waveform.setCellID(timeStruct.getCellID());

    // Fill the waveform with amplitude values
    // The sipm::SiPMAnalogSignal can be iterated as an std::vector<double>
    for (unsigned bin = 0; bin < anaSignal.size(); bin++) {
      double amp = anaSignal[bin];

      double tStart = static_cast<double>(bin) * m_sampling;
      double tEnd = static_cast<double>(bin + 1) * m_sampling;
      double center = (tStart + tEnd) / 2.;

      // Only include samples within our time window of interest
      if (center < timeStruct.getTime())
        continue;

      if (center > gateEnd)
        continue;

      waveform.addToAmplitude(amp);
    }
  }

  return StatusCode::SUCCESS;
}

StatusCode SimulateSiPMwithOpticalPhoton::finalize() {
  // Just delegate to the base class
  return Gaudi::Algorithm::finalize();
}
