#include "DualCrysCalDigi.h"
#include "GaudiKernel/Service.h"
#include "Randomize.h"
#include "DD4hep/CellIDPositionConverter.h"
#include <cmath>

DECLARE_COMPONENT(DUAL::DualCrysCalDigi)

namespace DUAL {

DualCrysCalDigi::DualCrysCalDigi(const std::string& name, ISvcLocator* svcLoc)
 : GaudiAlgorithm(name, svcLoc) {
  declareProperty("SimHitsA", m_simHitsA);
  declareProperty("SimHitsB", m_simHitsB);
  declareProperty("OutputHits", m_outputHits);
  declareProperty("GainA", m_gainA);
  declareProperty("GainB", m_gainB);
  declareProperty("Pedestal", m_pedestal);
  declareProperty("Noise", m_noise);
  declareProperty("AdcBits", m_adcBits);
}

StatusCode DualCrysCalDigi::initialize() {
  if (service("GeoSvc", m_geoSvc).isFailure()) {
    error() << "Unable to retrieve GeoSvc" << endmsg;
    return StatusCode::FAILURE;
  }
  info() << "DualCrysCalDigi: gainA=" << m_gainA
         << ", gainB=" << m_gainB
         << ", noise=" << m_noise
         << ", adcBits=" << m_adcBits << endmsg;
  return StatusCode::SUCCESS;
}

float DualCrysCalDigi::digitizeEnergy(double energy, double gain) const {
  double adc = energy * gain + m_pedestal;
  if (m_noise > 0)
    adc += gRandom->Gaus(0, m_noise);

  double maxADC = double((1u<<m_adcBits) - 1);
  adc = std::clamp(adc, 0.0, maxADC);
  return static_cast<float>(adc / gain);
}

StatusCode DualCrysCalDigi::execute() {
  const auto* simA = get<edm4hep::SimCalorimeterHitCollection>(m_simHitsA);
  const auto* simB = get<edm4hep::SimCalorimeterHitCollection>(m_simHitsB);
  auto* out = create<edm4hep::CalorimeterHitCollection>(m_outputHits);

  dd4hep::CellIDPositionConverter posTool(m_geoSvc->detector());

  // Process hits in crystal A
  for (const auto& h : *simA) {
    float e = digitizeEnergy(h.getEnergy(), m_gainA);
    if (e <= 0) continue;
    auto hit = out->create();
    hit.setCellID(h.getCellID());
    hit.setEnergy(e);
    hit.setTime(h.getTime());
    auto p = posTool.position(h.getCellID());
    hit.setPosition({float(p.x()), float(p.y()), float(p.z())});
  }

  // Process hits in crystal B
  for (const auto& h : *simB) {
    float e = digitizeEnergy(h.getEnergy(), m_gainB);
    if (e <= 0) continue;
    auto hit = out->create();
    hit.setCellID(h.getCellID());
    hit.setEnergy(e);
    hit.setTime(h.getTime());
    auto p = posTool.position(h.getCellID());
    hit.setPosition({float(p.x()), float(p.y()), float(p.z())});
  }

  return StatusCode::SUCCESS;
}

StatusCode DualCrysCalDigi::finalize() {
  info() << "DualCrysCalDigi finished." << endmsg;
  return StatusCode::SUCCESS;
}

} // namespace DUAL
