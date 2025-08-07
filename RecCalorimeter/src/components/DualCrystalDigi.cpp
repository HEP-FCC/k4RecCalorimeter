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

  ///if (!service("CellIDPositionConverter", m_cellIDPositionConverter)) {
  ///  error() << "Unable to locate CellID Position Converter Service." << endmsg;
  ///  return StatusCode::FAILURE;
  ///}  
  
  return StatusCode::SUCCESS;
}

float DualCrystalDigi::digitizeEnergy(double energy, double gain) const {
  // 1. Random engine
  std::random_device rd;                         // Non-deterministic seed
  std::default_random_engine generator(rd());    // Mersenne Twister engine

  // 2. Normal distribution: mean = 0.0, stddev = 1.0
  std::normal_distribution<double> gauss(0.0, 1.0);

  // 3. Generate 5 Gaussian-distributed values
  for (int i = 0; i < 5; ++i) {
    double sample = gauss(generator);
    std::cout << "Sample[" << i << "] = " << sample << '\n';
  }

  double adc = energy * gain + m_pedestal;
  /// if (m_noise > 0)
  ///   adc += gRandom->Gaus(0, m_noise);
  if (m_noise > 0) {
    // Static to avoid re-seeding every call — thread-safe alternative would use thread_local
    static thread_local std::default_random_engine engine(std::random_device{}());
    std::normal_distribution<float> gauss_dist(0.0f, m_noise);
    adc += gauss_dist(engine);
  }

  double maxADC = double((1u<<m_adcBits) - 1);
  adc = std::clamp(adc, 0.0, maxADC);
  return static_cast<float>(adc / gain);
}

StatusCode DualCrystalDigi::execute(const EventContext& /* ctx */ ) const {
  const auto* simA = m_simHitsA.get();
  const auto* simB = m_simHitsB.get();
  ///auto* out = create<edm4hep::CalorimeterHitCollection>(m_outputHits);
  auto out = m_outputHits.createAndPut();

  dd4hep::rec::CellIDPositionConverter posTool(*m_geoSvc->getDetector());

  // Process hits in crystal A
  for (const auto& h : *simA) {
    float e = digitizeEnergy(h.getEnergy(), m_gainA);
    if (e <= 0) continue;
    auto hit = out->create();
    hit.setCellID(h.getCellID());
    hit.setEnergy(e);
    /// FIXME: hit.setTime(h.getTime());
    auto p = posTool.position(h.getCellID());
    auto pos = edm4hep::Vector3f(float(p.x()), float(p.y()), float(p.z()));
    hit.setPosition(pos);
  }

  // Process hits in crystal B
  for (const auto& h : *simB) {
    float e = digitizeEnergy(h.getEnergy(), m_gainB);
    if (e <= 0) continue;
    auto hit = out->create();
    hit.setCellID(h.getCellID());
    hit.setEnergy(e);
    /// FIXME: hit.setTime(h.getTime());
    auto p = posTool.position(h.getCellID());
    auto pos = edm4hep::Vector3f(float(p.x()), float(p.y()), float(p.z()));
    hit.setPosition(pos);
  }

  return StatusCode::SUCCESS;
}

StatusCode DualCrystalDigi::finalize() {
  info() << "DualCrystalDigi finished." << endmsg;
  return StatusCode::SUCCESS;
}

} // namespace DUAL

DECLARE_COMPONENT(DUAL::DualCrystalDigi)

