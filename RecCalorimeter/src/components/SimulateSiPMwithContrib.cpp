#include "SimulateSiPMwithContrib.h"
#include "DD4hep/DD4hepUnits.h"
#include <cmath>

DECLARE_COMPONENT(SimulateSiPMwithContrib)

SimulateSiPMwithContrib::SimulateSiPMwithContrib(const std::string& aName, ISvcLocator* aSvcLoc)
    : Gaudi::Algorithm(aName, aSvcLoc) {
    
  // register DataHandles with hit collection names
  declareProperty("inputHitCollection", m_simHits, "input calo collection name");
  declareProperty("outputHitCollection", m_digiHits, "output calo collection name");
}

StatusCode SimulateSiPMwithContrib::initialize() {
  StatusCode sc = Gaudi::Algorithm::initialize();

  if (sc.isFailure())
    return sc;

  // Initialize random services
  m_randSvc = service("RndmGenSvc", false);

  if (!m_randSvc) {
    error() << "Couldn't get RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_rndmExp.initialize(m_randSvc, Rndm::Exponential(m_scintDecaytime.value())).isFailure()) {
    error() << "Couldn't initialize RndmGenSvc!" << endmsg;
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
  properties.setPdeType(sipm::SiPMProperties::PdeType::kNoPde);
  properties.setPdeSpectrum(m_wavelen, m_sipmEff);

  // set other parameters if provided
  for (const auto& [key, value] : m_params.value())
    properties.setProperty(key, value);

  m_sensor = std::make_unique<sipm::SiPMSensor>(properties); // must be constructed from SiPMProperties

  info() << "SimulateSiPMwithEdep initialized" << endmsg;
  info() << properties << endmsg; // sipm::SiPMProperties has std::ostream& operator<<

  return StatusCode::SUCCESS;
}

StatusCode SimulateSiPMwithContrib::execute(const EventContext&) const {
  const edm4hep::SimCalorimeterHitCollection* scintHits = m_simHits.get();
  edm4hep::CalorimeterHitCollection* digiHits = m_digiHits.createAndPut();

  // loop through the hits (each hit corresponds to a fiber)
  for (unsigned int idx = 0; idx < scintHits->size(); idx++) {
    const auto& scintHit = scintHits->at(idx);

    if (!scintHit.isAvailable()) {
    std::cout << "ERROR: Hit not available!" << std::endl;
    continue;
    }

    std::vector<double> vecTimes;
    std::vector<double> vecWavelens;
    
    // loop through the hit contributions
    for (auto contrib = scintHit.contributions_begin(); contrib != scintHit.contributions_end(); ++contrib) {
	
      if (!contrib->isAvailable()) {
        std::cerr << "ERROR: Hit contribution not available!" << std::endl;
        continue;
      }
      const double npe = contrib->getEnergy(); // in IDEA_o2 this is the number of photo-electrons

      vecTimes.reserve(vecTimes.size() + npe); // resize vector to store new p.e.
      vecWavelens.reserve(vecTimes.size() + npe);

      // get the photon arrival time at SiPM
      const double time = contrib->getTime(); // in IDEA_o2 this it the toa at SiPM in ns
      double thisTime = time;
      // if it is a scintillation hit, add scintillation decay time
      if (!m_isCherenkov){
        thisTime += m_rndmExp.shoot();
      }

      for (unsigned int pe=0; pe<npe; pe++){
        vecTimes.push_back(thisTime);
        vecWavelens.push_back(1.);
      }
    } // contrib

    // sort time of arrivals from shortest to longest
    std::sort(vecTimes.begin(), vecTimes.end());

    m_sensor->resetState();
    m_sensor->addPhotons(vecTimes, vecWavelens); // Sets photon times & wavelengths
    m_sensor->runEvent(); // Runs the simulation

    auto digiHit = digiHits->create();

    // Using only analog signal (ADC conversion is still experimental)
    const sipm::SiPMAnalogSignal anaSignal = m_sensor->signal();

    // if the signal never exceeds the threshold, it will return -1.
    const double integral =
        std::max(0., anaSignal.integral(m_gateStart, m_gateL, m_thres));           // (intStart, intGate, threshold)
    const double toa = std::max(0., anaSignal.toa(m_gateStart, m_gateL, m_thres)); // (intStart, intGate, threshold)

    digiHit.setEnergy(integral/* * m_scaleADC.value()*/); // convert ADC to GeV
    digiHit.setEnergyError(/*m_scaleADC.value() **/ std::sqrt(integral));
    digiHit.setPosition(scintHit.getPosition());
    digiHit.setCellID(scintHit.getCellID());
    // Toa and m_gateStart are in ns
    digiHit.setTime(toa + m_gateStart);
  }

  return StatusCode::SUCCESS;
}

StatusCode SimulateSiPMwithContrib::finalize() { return Gaudi::Algorithm::finalize(); }
