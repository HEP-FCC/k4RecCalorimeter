#include "SimulateSiPMwithOpticalPhoton.h"

// DECLARE_COMPONENT macro connects the algorithm to the framework
DECLARE_COMPONENT(SimulateSiPMwithOpticalPhoton)

SimulateSiPMwithOpticalPhoton::SimulateSiPMwithOpticalPhoton(const std::string& aName, ISvcLocator* aSvcLoc)
    : Gaudi::Algorithm(aName, aSvcLoc) {
  // Constructor doesn't need to do anything special
}

StatusCode SimulateSiPMwithOpticalPhoton::initialize() {
  // First call the base class initialization
  StatusCode sc = Gaudi::Algorithm::initialize();

  if (sc.isFailure())
    return sc;

  // Initialize random services
  m_randSvc = service("RndmGenSvc", false);

  if (!m_randSvc) {
    error() << "Couldn't get RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_rndmUniform.initialize(m_randSvc, Rndm::Flat(0., 1.)).isFailure()) {
    error() << "Couldn't initialize RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_wavelen.size() < 2) {
    error() << "SimulateSiPMwithOpticalPhoton: "
            << "The wavelength vector size must be greater or equal than 2" << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_wavelen.size() != m_sipmEff.size()) {
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
  properties.setPdeSpectrum(m_wavelen.value(), m_sipmEff.value());

  // set other parameters if provided
  for (const auto& [key, value] : m_params.value())
    properties.setProperty(key, value);

  // Create the SiPM sensor model
  m_sensor = std::make_unique<sipm::SiPMSensor>(properties);

  info() << "SimulateSiPMwithOpticalPhoton initialized" << endmsg;
  info() << properties << endmsg; // sipm::SiPMProperties has std::ostream& operator<<

  return StatusCode::SUCCESS;
}

std::vector<double> SimulateSiPMwithOpticalPhoton::integral(const std::vector<double>& wavelen,
                                                            const std::vector<double>& yval) const {
  double val = 0.;
  double prevXval = wavelen.front();
  std::vector<double> result = {0.};
  result.reserve(wavelen.size());

  for (unsigned idx = 1; idx < wavelen.size(); idx++) {
    double intervalX = std::abs(prevXval - wavelen.at(idx));
    double avgY = (yval.at(idx) + yval.at(idx - 1)) / 2.;
    val += intervalX * avgY;
    result.push_back(val);
  }

  return result;
}

StatusCode SimulateSiPMwithOpticalPhoton::execute(const EventContext&) const {
  // Get input collections
  const edm4hep::RawTimeSeriesCollection* timeStructs = m_timeStruct.get();
  const edm4hep::RawTimeSeriesCollection* waveLenStructs = m_wavelenStruct.get();
  const edm4hep::SimCalorimeterHitCollection* simHits = m_simHits.get();

  // Create output collections
  edm4hep::TimeSeriesCollection* waveforms = m_waveforms.createAndPut();
  edm4hep::CalorimeterHitCollection* digiHits = m_digiHits.createAndPut();

  // Process each hit
  for (unsigned int idx = 0; idx < timeStructs->size(); idx++) {
    const auto& timeStruct = timeStructs->at(idx);
    const auto& waveLength = waveLenStructs->at(idx);
    const auto& simHit = simHits->at(idx);

    // Validate that both collections have matching CellIDs
    if (timeStruct.getCellID() != simHit.getCellID()) {
      error() << "CellIDs of RawCalorimeterHits & RawTimeSeries are different! "
              << "This should never happen." << endmsg;
      return StatusCode::FAILURE;
    }

    // Validate that both collections have matching CellIDs
    if (waveLength.getCellID() != simHit.getCellID()) {
      error() << "CellIDs of RawCalorimeterHits & RawTimeSeries are different! "
              << "This should never happen." << endmsg;
      return StatusCode::FAILURE;
    }

    // extract wavelength randomly according to the wavelength distribution
    // where we used edm4hep::TimeSeries to store the distribution in the SD
    std::vector<double> vecWavelenEdm;
    std::vector<double> vecSpectrum;
    vecWavelenEdm.reserve(waveLength.adcCounts_size());
    vecSpectrum.reserve(waveLength.adcCounts_size());
    // be aware that we loop in the decreasing order!
    for (unsigned bin = waveLength.adcCounts_size(); bin != 1; bin--) {
      double wavlenBinCenter = waveLength.getTime() + waveLength.getInterval() * (static_cast<float>(bin - 1) + 0.5);
      double spectrum = static_cast<double>(waveLength.getAdcCounts(bin - 1));

      vecWavelenEdm.push_back(wavlenBinCenter);
      vecSpectrum.push_back(spectrum);
    }

    std::vector<double> integralSpectrum = integral(vecWavelenEdm, vecSpectrum);

    // extract photon arrival times from the time structure
    // and generate wavelength according to the pdf
    std::vector<double> vecTimes;
    std::vector<double> vecWavelengths;
    vecTimes.reserve(static_cast<unsigned>(simHit.getEnergy()));
    vecWavelengths.reserve(static_cast<unsigned>(simHit.getEnergy()));

    // SimSiPM ignores negative time photons, so translate the whole time structure if needed
    double minTime = 0.;

    for (unsigned int bin = 0; bin < timeStruct.adcCounts_size(); bin++) {
      int counts = static_cast<int>(timeStruct.getAdcCounts(bin));

      if (counts == 0)
        continue;

      double timeBin = timeStruct.getTime() + timeStruct.getInterval() * static_cast<float>(bin);

      if (timeBin < minTime)
        minTime = timeBin;
    }

    // now fill the photon vectors
    for (unsigned int bin = 0; bin < timeStruct.adcCounts_size(); bin++) {
      int counts = static_cast<int>(timeStruct.getAdcCounts(bin));
      double timeBin = timeStruct.getTime() + timeStruct.getInterval() * (static_cast<float>(bin) + 0.5);

      // Add a photon arrival time for each count in this bin
      for (int num = 0; num < counts; num++) {
        vecTimes.emplace_back(timeBin - minTime); // shift to non-negative time

        // generate wavelength
        const double randval = integralSpectrum.back() * m_rndmUniform.shoot();
        unsigned xhigh = 1;

        for (; xhigh < integralSpectrum.size() - 1; xhigh++) {
          if (randval < integralSpectrum.at(xhigh))
            break;
        }

        const unsigned xlow = xhigh - 1;
        const double xdiff = vecWavelenEdm.at(xhigh) - vecWavelenEdm.at(xlow);
        const double ydiff = integralSpectrum.at(xhigh) - integralSpectrum.at(xlow);
        const double ydiffInv = (ydiff == 0.) ? 0. : 1. / ydiff;
        double valWav = vecWavelenEdm.at(xlow) + xdiff * (randval - integralSpectrum.at(xlow)) * ydiffInv;

        vecWavelengths.emplace_back(valWav);
      }
    }

    // Reset the SiPM state and run the simulation
    m_sensor->resetState();
    m_sensor->addPhotons(vecTimes, vecWavelengths); // Sets photon arrival times (in ns) & wavelengths (in nm)
    m_sensor->runEvent();                           // Runs the SiPM simulation

    auto digiHit = digiHits->create();
    auto waveform = waveforms->create();

    // Get the analog signal from the SiPM
    const sipm::SiPMAnalogSignal anaSignal = m_sensor->signal();

    // Compute integral (energy) and time of arrival
    // if the signal never exceeds the threshold, it will return -1.
    const double integral =
        std::max(0., anaSignal.integral(m_gateStart - minTime, m_gateL, m_thres));           // (intStart, intGate, threshold)
    const double toa = std::max(0., anaSignal.toa(m_gateStart - minTime, m_gateL, m_thres)); // (intStart, intGate, threshold)

    // Set digitized hit properties
    digiHit.setEnergy(integral * m_scaleADC.value());
    digiHit.setEnergyError(m_scaleADC.value() * std::sqrt(integral));
    digiHit.setPosition(simHit.getPosition());
    digiHit.setCellID(simHit.getCellID());
    digiHit.setTime(toa + m_gateStart); // Toa and m_gateStart are in ns

    // Set waveform properties
    waveform.setInterval(m_sampling);
    waveform.setTime(m_storeFullWaveform.value() ? minTime : toa + m_gateStart);
    waveform.setCellID(timeStruct.getCellID());

    // Fill the waveform with amplitude values
    // The sipm::SiPMAnalogSignal can be iterated as an std::vector<double>
    const double gateEnd = m_gateStart.value() + m_gateL.value();

    if (integral > 0.) {
      for (unsigned bin = 0; bin < anaSignal.size(); bin++) {
        float amp = anaSignal[bin];

        double tStart = static_cast<double>(bin) * m_sampling + minTime;
        double tEnd = static_cast<double>(bin + 1) * m_sampling + minTime;
        double center = (tStart + tEnd) / 2.;

        // Only include samples within our time window of interest
        if (!m_storeFullWaveform.value()) {
          if (center < toa + m_gateStart)
            continue;

          if (center > gateEnd)
            continue;

          if (amp < m_thres)
            continue;
        }

        waveform.addToAmplitude(amp);
      }
    }
  }

  return StatusCode::SUCCESS;
}

StatusCode SimulateSiPMwithOpticalPhoton::finalize() {
  // Just delegate to the base class
  return Gaudi::Algorithm::finalize();
}
