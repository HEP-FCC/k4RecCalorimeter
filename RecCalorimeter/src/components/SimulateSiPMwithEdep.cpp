#include "SimulateSiPMwithEdep.h"

#include "DDG4/Geant4Random.h"
#include "DD4hep/DD4hepUnits.h"

#include <cmath>
#include <stdexcept>

DECLARE_COMPONENT(SimulateSiPMwithEdep)

SimulateSiPMwithEdep::SimulateSiPMwithEdep(const std::string& aName, ISvcLocator* aSvcLoc) : Gaudi::Algorithm(aName, aSvcLoc) {}

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

  if (m_rndmUniform.initialize(m_randSvc, Rndm::Flat(0.,1.)).isFailure()) {
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

  if (m_wavelen.size()!=m_sipmEff.size()) {
    error() << "SimulateSiPMwithEdep: "
            << "The SiPM efficiency vector size should be equal to the wavelength vector size" << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_wavelen.size()!=m_scintSpectrum.size()) {
    error() << "SimulateSiPMwithEdep: "
               "The scintillation spectrum vector size should be equal to the wavelength vector size" << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_wavelen.size()!=m_filterEff.size()) {
    error() << "SimulateSiPMwithEdep: "
               "The filter efficiency vector size should be equal to the wavelength vector size" << endmsg;
    return StatusCode::FAILURE;
  }

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
  properties.setPdeSpectrum(m_wavelen,m_sipmEff);

  m_sensor = std::make_unique<sipm::SiPMSensor>(properties); // must be constructed from SiPMProperties

  info() << "SimulateSiPMwithEdep initialized" << endmsg;
  info() << properties << endmsg; // sipm::SiPMProperties has std::ostream& operator<<

  // calculate the total filtering efficiency
  std::vector<double> specTimesEff = m_filterEff; // copy

  for (unsigned ibin = 0; ibin < specTimesEff.size(); ibin++)
    specTimesEff.at(ibin) *= m_scintSpectrum.value().at(ibin);

  // integral scintillation spectrum times efficiency
  // and retrieve average efficiency
  auto integralSpectrum = integral(m_wavelen.value(),m_scintSpectrum.value());
  m_integral = integral(m_wavelen.value(),specTimesEff);
  m_efficiency = m_integral.back()/integralSpectrum.back();

  return StatusCode::SUCCESS;
}

std::vector<double> SimulateSiPMwithEdep::integral(const std::vector<double>& wavelen, const std::vector<double>& yval) const {
  double val = 0.;
  double prevXval = wavelen.front();
  std::vector<double> result = {0.};
  result.reserve(wavelen.size());

  for (unsigned idx = 1; idx < wavelen.size(); idx++) {
    double intervalX = std::abs(prevXval - wavelen.at(idx));
    double avgY = (yval.at(idx) + yval.at(idx-1))/2.;
    val += intervalX*avgY;
    result.push_back(val);
  }

  return result;
}

StatusCode SimulateSiPMwithEdep::execute(const EventContext&) const {
  const edm4hep::SimCalorimeterHitCollection* scintHits = m_scintHits.get();
  edm4hep::CalorimeterHitCollection* digiHits = m_digiHits.createAndPut();

  for (unsigned int idx = 0; idx < scintHits->size(); idx++) {
    const auto& scintHit = scintHits->at(idx);
    const double edep = scintHit.getEnergy()*dd4hep::GeV;
    const double yield = m_scintYield.value()/dd4hep::keV;

    double avgNphoton = edep*yield*m_efficiency;

    // generate the number of p.e. (npe)
    unsigned npe = 0;

    if (avgNphoton < 10.) {
      Rndm::Numbers rnd(m_randSvc, Rndm::Poisson(avgNphoton));
      npe = static_cast<unsigned>(std::floor(rnd.shoot() + 0.5));
    } else {
      Rndm::Numbers rnd(m_randSvc, Rndm::Gauss(avgNphoton,std::sqrt(avgNphoton)));
      npe = static_cast<unsigned>(std::floor(rnd.shoot() + 0.5));
    }

    std::vector<double> vecTimes;
    std::vector<double> vecWavelens;
    vecTimes.reserve(npe);
    vecWavelens.reserve(npe);

    for (unsigned ipho = 0; ipho < npe; ipho++) {
      // get photon wavelength
      // similar to https://gitlab.cern.ch/geant4/geant4/-/blob/master/source/processes/electromagnetic/xrays/src/G4Scintillation.cc
      const double randval = m_integral.back()*m_rndmUniform.shoot();
      unsigned xhigh = 1;

      for (xhigh = 1; xhigh < m_integral.size(); xhigh++) {
        if (randval < m_integral.at(xhigh))
          break;
      }

      const unsigned xlow = xhigh - 1;
      const double xdiff = m_wavelen.value().at(xhigh) - m_wavelen.value().at(xlow);
      const double ydiff = m_integral.at(xhigh) - m_integral.at(xlow);
      const double ydiffInv = (ydiff==0.) ? 0. : 1./ydiff;
      double valWav = m_wavelen.value().at(xlow) + xdiff*(randval-m_integral.at(xlow))*ydiffInv;

      double scintTime = m_timeInterval.value() + m_rndmExp.shoot();

      vecTimes.push_back(scintTime);
      vecWavelens.push_back(valWav);
    }

    m_sensor->resetState();
    m_sensor->addPhotons(vecTimes,vecWavelens); // Sets photon times & wavelengths
    m_sensor->runEvent(); // Runs the simulation

    auto digiHit = digiHits->create();

    // Using only analog signal (ADC conversion is still experimental)
    const sipm::SiPMAnalogSignal anaSignal = m_sensor->signal();

    const double integral = anaSignal.integral(m_gateStart,m_gateL,m_thres); // (intStart, intGate, threshold)
    const double toa = anaSignal.toa(m_gateStart,m_gateL,m_thres);           // (intStart, intGate, threshold)

    digiHit.setEnergy( integral );
    digiHit.setCellID( scintHit.getCellID() );
    // Toa and m_gateStart are in ns
    digiHit.setTime( toa+m_gateStart );
  }

  return StatusCode::SUCCESS;
}

StatusCode SimulateSiPMwithEdep::finalize() {
  return Gaudi::Algorithm::finalize();
}
