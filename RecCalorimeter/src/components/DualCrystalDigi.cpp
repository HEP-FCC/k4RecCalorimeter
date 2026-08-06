#include "DualCrystalDigi.h"
#include "GaudiKernel/Service.h"
#include "DDRec/CellIDPositionConverter.h"
#include "CLHEP/Random/Randomize.h"
#include <cmath>
#include <random>

namespace DUAL {

DualCrystalDigi::DualCrystalDigi(const std::string& name, ISvcLocator* svcLoc)
 : Gaudi::Algorithm(name, svcLoc) {
  declareProperty("SimHitsA", m_simHitsA);
  declareProperty("SimHitsB", m_simHitsB);
  declareProperty("OutputHits", m_outputHits);
  declareProperty("GainA", m_gainA);
  declareProperty("GainB", m_gainB);
  declareProperty("Pedestal", m_pedestal);
  declareProperty("Noise", m_noise);
  declareProperty("AdcBits", m_adcBits);
}

StatusCode DualCrystalDigi::initialize() {
  if (!service("GeoSvc", m_geoSvc)) {
    error() << "Unable to locate Geometry Service." << endmsg;
    return StatusCode::FAILURE;
  }
  info() << "DualCrysCalDigi: gainA=" << m_gainA
         << ", gainB=" << m_gainB
         << ", noise=" << m_noise
         << ", adcBits=" << m_adcBits << endmsg;

  return StatusCode::SUCCESS;
}

float DualCrystalDigi::digitizeEnergy(double energy, double gain) const {
  double adc = energy * gain + m_pedestal;

  if (m_noise > 0) {
    static thread_local std::default_random_engine engine(std::random_device{}());
    std::normal_distribution<float> gauss_dist(0.0f, m_noise);
    adc += gauss_dist(engine);
  }

  double maxADC = double((1u << m_adcBits) - 1);
  adc = std::clamp(adc, 0.0, maxADC);
  return static_cast<float>(adc / gain);
}

StatusCode DualCrystalDigi::execute(const EventContext& /* ctx */ ) const {
  const auto* simA = m_simHitsA.get();
  const auto* simB = m_simHitsB.get();
  auto out = m_outputHits.createAndPut();

  dd4hep::rec::CellIDPositionConverter posTool(*m_geoSvc->getDetector());

  processHits(simA, m_gainA, out, posTool);
  processHits(simB, m_gainB, out, posTool);

  return StatusCode::SUCCESS;
}

void DualCrystalDigi::processHits(const edm4hep::SimCalorimeterHitCollection* simHits,
                                  double gain,
                                  edm4hep::CalorimeterHitCollection* out,
                                  dd4hep::rec::CellIDPositionConverter& posTool) const {
  for (const auto& h : *simHits) {
    float e = digitizeEnergy(h.getEnergy(), gain);
    if (e <= 0) continue;

    auto hit = out->create();
    hit.setCellID(h.getCellID());
    hit.setEnergy(e);

    // Try to set time if supported, else comment out or set to zero
    // Uncomment the following line if your edm4hep version supports time()
    // hit.setTime(h.time());

    auto p = posTool.position(h.getCellID());
    hit.setPosition(edm4hep::Vector3f(float(p.x()), float(p.y()), float(p.z())));
  }
}

StatusCode DualCrystalDigi::finalize() {
  info() << "DualCrystalDigi finished." << endmsg;
  return StatusCode::SUCCESS;
}

} // namespace DUAL

DECLARE_COMPONENT(DUAL::DualCrystalDigi)
