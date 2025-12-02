#include "SimulateSiPMwithEdep.h"
#include "DD4hep/DD4hepUnits.h"
#include <cmath>

DECLARE_COMPONENT(SimulateSiPMwithEdep)

SimulateSiPMwithEdep::SimulateSiPMwithEdep(const std::string& aName, ISvcLocator* aSvcLoc)
    : Gaudi::Algorithm(aName, aSvcLoc) {}

StatusCode SimulateSiPMwithEdep::initialize() {
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

  if (m_rndmExp.initialize(m_randSvc, Rndm::Exponential(m_scintDecaytime.value())).isFailure()) {
    error() << "Couldn't initialize RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_wavelen.size() < 2) {
    error() << "SimulateSiPMwithEdep: "
            << "The wavelength vector size must be greater or equal than 2" << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_wavelen.size() != m_sipmEff.size()) {
    error() << "SimulateSiPMwithEdep: "
            << "The SiPM efficiency vector size should be equal to the wavelength vector size" << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_wavelen.size() != m_scintSpectrum.size()) {
    error() << "SimulateSiPMwithEdep: "
               "The scintillation spectrum vector size should be equal to the wavelength vector size"
            << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_wavelen.size() != m_filterEff.size()) {
    error() << "SimulateSiPMwithEdep: "
               "The filter efficiency vector size should be equal to the wavelength vector size"
            << endmsg;
    return StatusCode::FAILURE;
  }

  m_geoSvc = service("GeoSvc");

  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service. "
            << "Make sure you have GeoSvc in the configuration." << endmsg;

    return StatusCode::FAILURE;
  }

  if (m_geoSvc->getDetector()->readouts().find(m_readoutName) ==
      m_geoSvc->getDetector()->readouts().end()) {
    error() << "Readout <<" << m_readoutName << ">> does not exist." << endmsg;

    return StatusCode::FAILURE;
  }

  // get segmentation (cast to type GridDRcalo_k4geo)
  m_segmentation = dynamic_cast<dd4hep::DDSegmentation::GridDRcalo_k4geo*>(m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());

  // initialize SiPM properties
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
  properties.setPdeType(sipm::SiPMProperties::PdeType::kSpectrumPde);
  properties.setPdeSpectrum(m_wavelen, m_sipmEff);

  // set other parameters if provided
  for (const auto& [key, value] : m_params.value())
    properties.setProperty(key, value);

  m_sensor = std::make_unique<sipm::SiPMSensor>(properties); // must be constructed from SiPMProperties

  info() << "SimulateSiPMwithEdep initialized" << endmsg;
  info() << properties << endmsg; // sipm::SiPMProperties has std::ostream& operator<<

  // calculate the total filtering efficiency
  std::vector<double> specTimesEff = m_filterEff; // copy

  for (unsigned ibin = 0; ibin < specTimesEff.size(); ibin++)
    specTimesEff.at(ibin) *= m_scintSpectrum.value().at(ibin);

  // integral scintillation spectrum times efficiency
  // and retrieve average efficiency
  auto integralSpectrum = integral(m_wavelen.value(), m_scintSpectrum.value());
  m_integral = integral(m_wavelen.value(), specTimesEff);
  m_efficiency = m_integral.back() / integralSpectrum.back();

  return StatusCode::SUCCESS;
}

std::vector<double> SimulateSiPMwithEdep::integral(const std::vector<double>& wavelen,
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

StatusCode SimulateSiPMwithEdep::execute(const EventContext&) const {
  const edm4hep::SimCalorimeterHitCollection* scintHits = m_scintHits.get();
  edm4hep::CalorimeterHitCollection* digiHits = m_digiHits.createAndPut();
  edm4hep::TimeSeriesCollection* waveforms = m_waveforms.createAndPut();
  edm4hep::CaloHitSimCaloHitLinkCollection* hitLinks = m_hitLinks.createAndPut();

  const double yield = m_scintYield.value() / dd4hep::keV;

  for (unsigned int idx = 0; idx < scintHits->size(); idx++) {
    const auto& scintHit = scintHits->at(idx);
    const auto sipmPos = m_segmentation->sipmPosition(scintHit.getCellID()); // in dd4hep unit

    std::vector<double> vecTimes;
    std::vector<double> vecWavelens;

    // SimSiPM ignores negative time photons, so translate the whole time structure if needed
    double minTime = 0.;

    for (auto contrib = scintHit.contributions_begin(); contrib != scintHit.contributions_end(); ++contrib) {
      const double edep = contrib->getEnergy() * dd4hep::GeV;
      double avgNphoton = edep * yield * m_efficiency;

      // generate the number of p.e. (npe)
      unsigned npe = 0;

      if (avgNphoton < 10.) {
        Rndm::Numbers rnd(m_randSvc, Rndm::Poisson(avgNphoton));
        npe = static_cast<unsigned>(std::floor(rnd.shoot() + 0.5));
      } else {
        Rndm::Numbers rnd(m_randSvc, Rndm::Gauss(avgNphoton, std::sqrt(avgNphoton)));
        // prevent underflow since Gaussian can shoot negative in rare case
        double val = std::floor(rnd.shoot() + 0.5);
        npe = static_cast<unsigned>(std::max(val, 0.));
      }

      vecTimes.reserve(vecTimes.size() + npe);
      vecWavelens.reserve(vecTimes.size() + npe);

      // calculate photon arrival time at SiPM
      // using the distance from the step to the SiPM
      const double initialTime = contrib->getTime(); // in ns
      const auto stepPos = contrib->getStepPosition(); // in edm4hep unit

      const float relX = sipmPos.x() - stepPos.x * dd4hep::mm;
      const float relY = sipmPos.y() - stepPos.y * dd4hep::mm;
      const float relZ = sipmPos.z() - stepPos.z * dd4hep::mm;
      const float dist = std::sqrt(relX * relX + relY * relY + relZ * relZ); // in dd4hep unit

      const double effSpeedOfLight = dd4hep::c_light / m_refractiveIndex.value();
      const double distTime = dist / effSpeedOfLight;
      const double arrivalTime =
          m_switchTime.value() ? initialTime : initialTime + distTime / dd4hep::nanosecond; // in ns

      for (unsigned ipho = 0; ipho < npe; ipho++) {
        // get photon wavelength
        // similar to
        // https://gitlab.cern.ch/geant4/geant4/-/blob/master/source/processes/electromagnetic/xrays/src/G4Scintillation.cc
        const double randval = m_integral.back() * m_rndmUniform.shoot();
        unsigned xhigh = 1;

        for (xhigh = 1; xhigh < m_integral.size() - 1; xhigh++) {
          if (randval < m_integral.at(xhigh))
            break;
        }

        // linear interpolation
        const unsigned xlow = xhigh - 1;
        const double xdiff = m_wavelen.value().at(xhigh) - m_wavelen.value().at(xlow);
        const double ydiff = m_integral.at(xhigh) - m_integral.at(xlow);
        const double ydiffInv = (ydiff == 0.) ? 0. : 1. / ydiff;
        double valWav = m_wavelen.value().at(xlow) + xdiff * (randval - m_integral.at(xlow)) * ydiffInv;

        // get absorption length
        unsigned iwave = 1;

        for (iwave = 1; iwave < m_wavelen.value().size() - 1; iwave++) { // decreasing order
          if (valWav > m_wavelen.value().at(iwave))
            break;
        }

        // absorption length is ill-defined if scintHit.getPosition() is not the sensor position
        if (!m_switchTime.value()) {
          // linear interpolation
          const unsigned waveLow = iwave - 1;
          const double waveDiff = m_wavelen.value().at(iwave) - m_wavelen.value().at(waveLow);
          const double waveDiffInv = (waveDiff == 0.) ? 0. : 1. / waveDiff;
          const double abslenDiff = m_absLen.value().at(iwave) - m_absLen.value().at(waveLow);
          double absLen =
              m_absLen.value().at(waveLow) + abslenDiff * waveDiffInv * (valWav - m_wavelen.value().at(waveLow));

          // check absorption
          // similar to https://gitlab.cern.ch/geant4/geant4/-/blob/master/source/processes/management/src/G4VProcess.cc
          const double nInteractionLengthLeft = -std::log(m_rndmUniform.shoot());
          const double nInteractionLength = dist / (absLen * dd4hep::meter);

          // absorb photons
          if (nInteractionLength > nInteractionLengthLeft)
            continue;
        }

        // get scintillation time
        double scintTime = arrivalTime + m_rndmExp.shoot();

        if (scintTime < minTime)
          minTime = scintTime;

        vecTimes.push_back(scintTime);
        vecWavelens.push_back(valWav);
      } // ipho
    } // contrib

    // shift times to non-negative
    std::vector<double> vecTimesShifted(vecTimes.size());
    std::transform(vecTimes.begin(), vecTimes.end(), vecTimesShifted.begin(),
                   [minTime](double t) { return t - minTime; });

    m_sensor->resetState();
    m_sensor->addPhotons(vecTimesShifted, vecWavelens); // Sets photon times & wavelengths
    m_sensor->runEvent();                        // Runs the simulation

    auto digiHit = digiHits->create();
    auto waveform = waveforms->create();
    auto hitLink = hitLinks->create();
    hitLink.setFrom(digiHit);
    hitLink.setTo(scintHit);

    // Using only analog signal (ADC conversion is still experimental)
    const sipm::SiPMAnalogSignal anaSignal = m_sensor->signal();

    // if the signal never exceeds the threshold, it will return -1.
    const double integral =
        std::max(0., anaSignal.integral(m_gateStart - minTime, m_gateL, m_thres));           // (intStart, intGate, threshold)
    const double toa = std::max(0., anaSignal.toa(m_gateStart - minTime, m_gateL, m_thres)); // (intStart, intGate, threshold)

    digiHit.setEnergy(integral * m_scaleADC.value());
    digiHit.setEnergyError(m_scaleADC.value() * std::sqrt(integral));
    digiHit.setPosition(scintHit.getPosition());
    digiHit.setCellID(scintHit.getCellID());
    // Toa and m_gateStart are in ns
    digiHit.setTime(toa + m_gateStart);

    // Set waveform properties
    waveform.setInterval(m_sampling);
    waveform.setTime(m_storeFullWaveform.value() ? minTime : toa + m_gateStart);
    waveform.setCellID(scintHit.getCellID());

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

StatusCode SimulateSiPMwithEdep::finalize() { return Gaudi::Algorithm::finalize(); }
