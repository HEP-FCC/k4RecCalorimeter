#include "CaloDIGI.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "Randomize.h"
#include "DD4hep/Detector.h"
#include "DD4hep/CellIDPositionConverter.h"
#include <cmath>

DECLARE_COMPONENT(RECAL::CaloDIGI)

namespace RECAL {

CaloDIGI::CaloDIGI(const std::string& name, ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc) {
  declareProperty("InputSimHits", m_inputSimHits);
  declareProperty("OutputHits", m_outputHits);
  declareProperty("ADCgain", m_adcGain);
  declareProperty("Pedestal", m_pedestal);
  declareProperty("NoiseRMS", m_noiseRMS);
  declareProperty("TimeSmear", m_timeSmear);
  declareProperty("EnergyThreshold", m_energyThreshold);
  declareProperty("ADCbits", m_adcBits);
  declareProperty("Verbose", m_verbose);
}

StatusCode CaloDIGI::initialize() {
  if (service("GeoSvc", m_geoSvc).isFailure()) {
    error() << "Unable to locate Geometry Service." << endmsg;
    return StatusCode::FAILURE;
  }

  info() << "Initialized CaloDIGI with settings:"
         << " gain=" << m_adcGain
         << ", pedestal=" << m_pedestal
         << ", noiseRMS=" << m_noiseRMS
         << ", timeSmear=" << m_timeSmear
         << ", threshold=" << m_energyThreshold
         << ", adcBits=" << m_adcBits
         << endmsg;
  return StatusCode::SUCCESS;
}

double CaloDIGI::digitizeToMeV(double energyGeV) const {
  double adc = energyGeV * m_adcGain + m_pedestal;
  if (m_noiseRMS > 0)
    adc += gRandom->Gaus(0, m_noiseRMS);

  adc = std::max(0.0, std::min(adc, static_cast<double>((1 << m_adcBits) - 1)));

  return (adc / m_adcGain) * 1000.0;  // MeV
}

double CaloDIGI::calculateTOF(double radius) const {
  constexpr double speedOfLight = 299.792458;  // mm/ns
  return radius / speedOfLight;
}

double CaloDIGI::smearTime(double time) const {
  return time + gRandom->Gaus(0, m_timeSmear);
}

StatusCode CaloDIGI::execute() {
  const auto* simHits = get<edm4hep::SimCalorimeterHitCollection>(m_inputSimHits);
  auto* caloHits = create<edm4hep::CalorimeterHitCollection>(m_outputHits);

  dd4hep::CellIDPositionConverter posTool(m_geoSvc->detector());

  for (const auto& simHit : *simHits) {
    double energyMeV = digitizeToMeV(simHit.getEnergy());
    if (energyMeV < m_energyThreshold) continue;

    auto pos = posTool.position(simHit.getCellID());
    double tof = calculateTOF(pos.r());
    double time = smearTime(simHit.getTime() - tof);

    auto hit = caloHits->create();
    hit.setCellID(simHit.getCellID());
    hit.setEnergy(energyMeV);
    hit.setTime(time);
    hit.setPosition({(float)pos.x(), (float)pos.y(), (float)pos.z()});

    if (m_verbose) {
      debug() << "Hit cellID=" << simHit.getCellID()
              << " rawEnergy=" << simHit.getEnergy()
              << " MeV=" << energyMeV
              << " time=" << time
              << " pos=(" << pos.x() << "," << pos.y() << "," << pos.z() << ")"
              << endmsg;
    }
  }

  return StatusCode::SUCCESS;
}

StatusCode CaloDIGI::finalize() {
  info() << "CaloDIGI finalized." << endmsg;
  return StatusCode::SUCCESS;
}

} // namespace RECAL
