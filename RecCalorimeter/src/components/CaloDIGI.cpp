#include "CaloDIGI.h"

#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/RegistryEntry.h"

#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"

#include "DD4hep/Detector.h"
#include "DD4hep/IDDescriptor.h"

#include <random>
#include <cmath>

DECLARE_COMPONENT(RECAL::CaloDIGI)

namespace RECAL {

CaloDIGI::CaloDIGI(const std::string& name, ISvcLocator* svcLoc)
    : Gaudi::Algorithm(name, svcLoc) {}

StatusCode CaloDIGI::initialize() {
  if (Gaudi::Algorithm::initialize().isFailure()) return StatusCode::FAILURE;

  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service." << endmsg;
    return StatusCode::FAILURE;
  }

  info() << "Digitization tool initialized" << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode CaloDIGI::execute(const EventContext&) const {
  const auto* simHits = m_inputSimHits.get();
  if (!simHits) {
    error() << "SimCalorimeterHitCollection not found." << endmsg;
    return StatusCode::FAILURE;
  }
  auto* digiHits = m_outputHits.createAndPut();

  std::default_random_engine generator;
  std::normal_distribution<double> noise_dist(0.0, m_noiseRMS.value());
  std::normal_distribution<double> time_dist(0.0, m_timeSmear.value());

  for (const auto& hit : *simHits) {
    double energyGeV = hit.getEnergy();
    double energyMeV = digitizeToMeV(energyGeV);

    if (energyMeV < m_energyThreshold.value()) continue;

    double adc = energyMeV * m_adcGain.value() + m_pedestal.value();
    adc += noise_dist(generator);

    // Clip to ADC bit range
    double adcMax = std::pow(2, m_adcBits.value()) - 1;
    adc = std::max(0.0, std::min(adc, adcMax));

    // double time = hit.getTime() + time_dist(generator);

    auto digiHit = digiHits->create();
    digiHit.setCellID(hit.getCellID());
    digiHit.setEnergy(adc / m_adcGain.value()); // Convert back to MeV
    digiHit.setPosition(hit.getPosition());
    // digiHit.setTime(time);

    // digiHit.setEnergy(integral * m_scaleADC.value());
    // digiHit.setEnergyError(m_scaleADC.value() * std::sqrt(integral));
    // digiHit.setPosition(scintHit.getPosition());
    // Toa and m_gateStart are in ns
    // digiHit.setTime(toa + m_gateStart);

    if (m_verbose) {
      info() << "SimHit E=" << energyGeV << " GeV -> "
             << "ADC=" << adc << ", time=" << time << " ns" << endmsg;
    }
  }

  return StatusCode::SUCCESS;
}

StatusCode CaloDIGI::finalize() {
  return Gaudi::Algorithm::finalize();
}

double CaloDIGI::digitizeToMeV(double energyGeV) const {
  return energyGeV * 1e3;
}

double CaloDIGI::smearTime(double time) const {
  static std::default_random_engine generator;
  std::normal_distribution<double> dist(0.0, m_timeSmear.value());
  return time + dist(generator);
}

double CaloDIGI::calculateTOF(double radius) const {
  constexpr double speedOfLight = Gaudi::Units::c_light; // mm/ns
  return radius / speedOfLight;
}

} // namespace RECAL
