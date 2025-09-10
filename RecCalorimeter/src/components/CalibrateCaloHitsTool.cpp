#include "CalibrateCaloHitsTool.h"
#include "RecCaloCommon/k4RecCalorimeter_check.h"

DECLARE_COMPONENT(CalibrateCaloHitsTool)

StatusCode CalibrateCaloHitsTool::initialize() {
  K4RECCALORIMETER_CHECK( AlgTool::initialize() );

  info() << "Calibration constant: 1/sampling fraction=" << m_invSamplingFraction << endmsg;
  return StatusCode::SUCCESS;
}

void CalibrateCaloHitsTool::calibrate(std::unordered_map<uint64_t, double>& aHits) const {
  // Loop through energy deposits, multiply energy to get cell energy at electromagnetic scale
  for (auto& p : aHits) {
    p.second *= m_invSamplingFraction;
  }
}

void CalibrateCaloHitsTool::calibrate(std::vector<std::pair<uint64_t, double> >& aHits) const {
  // Loop through energy deposits, multiply energy to get cell energy at electromagnetic scale
  for (auto& p : aHits) {
    p.second *= m_invSamplingFraction;
  }
}

