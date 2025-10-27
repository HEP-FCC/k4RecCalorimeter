#include "CalibrateInLayersTool.h"

// k4FWCore
#include "k4Interface/IGeoSvc.h"
#include "RecCaloCommon/k4RecCalorimeter_check.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DD4hep/Readout.h"

DECLARE_COMPONENT(CalibrateInLayersTool)

StatusCode CalibrateInLayersTool::initialize() {
  K4RECCALORIMETER_CHECK( AlgTool::initialize() );
  // check if readout exists
  if (m_geoSvc->getDetector()->readouts().find(m_readoutName) == m_geoSvc->getDetector()->readouts().end()) {
    error() << "Readout <<" << m_readoutName << ">> does not exist." << endmsg;
    return StatusCode::FAILURE;
  }
  m_decoder = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder();
  return StatusCode::SUCCESS;
}

void CalibrateInLayersTool::calibrate(std::unordered_map<uint64_t, double>& aHits) const {
  // Loop through energy deposits, multiply energy to get cell energy at electromagnetic scale
  for (auto& p : aHits) {
    calibrateCell (p.first, p.second);
  }
}

void CalibrateInLayersTool::calibrate(std::vector<std::pair<uint64_t, double> >& aHits) const {
  for (auto& p : aHits) {
    calibrateCell (p.first, p.second);
  }
}

void CalibrateInLayersTool::calibrateCell(uint64_t cID, double& energy) const
{
  // shift layer id if the numbering does not start at 0
  uint layer = m_decoder->get(cID, m_layerFieldName) - m_firstLayerId;
  if (layer < m_samplingFraction.size()) {
    energy /= m_samplingFraction[layer];
  } else {
    energy /= m_samplingFraction[m_samplingFraction.size() - 1];
    warning() << "Size of sampling fraction values is smaller than the number of existing layers."
              << " Taking the sampling fraction for last layer."
              << " Layer ID: " << layer << endmsg;
  }
}

