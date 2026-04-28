#include "ConstNoiseTool.h"
#include "RecCaloCommon/k4RecCalorimeter_check.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"

// k4FWCore
#include "k4Interface/IGeoSvc.h"

// DD4hep
#include "DD4hep/Detector.h"

DECLARE_COMPONENT(ConstNoiseTool)

StatusCode ConstNoiseTool::initialize() {
  // initialize decoder for retrieving system ID
  m_decoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(m_systemEncoding);
  m_systemIndex = m_decoder->index("system");

  // check that sizes of m_detectors and m_detectorNoiseRMS and m_detectorNoiseOffset are the same
  if (m_detectorsNoiseRMS.size() != m_detectors.size()) {
    error() << "Sizes of the detectors vector and detectorsNoiseRMS vector do not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_setNoiseOffset) {
    if (m_detectorsNoiseOffset.size() != m_detectors.size()) {
      error() << "Sizes of the detectors vector and detectorsNoiseOffset vector do not match, exiting!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  K4RECCALORIMETER_CHECK(m_geoSvc.retrieve());

  // loop over the detectors
  for (size_t iDet = 0; iDet < m_detectors.size(); iDet++) {
    std::string detName = "DetID_" + m_detectors[iDet];
    try {
      int detID = m_geoSvc->getDetector()->constant<int>(detName);
      if (detID < 0) {
        error() << "Bad detector ID " << detID << endmsg;
        return StatusCode::FAILURE;
      }
      if (static_cast<int>(m_noise.size()) <= detID)
        m_noise.resize(detID + 1);
      m_noise[detID].first = m_detectorsNoiseRMS[iDet];
      debug() << "Set noise RMS for detector " << detName << " (ID=" << detID << ") to: " << m_detectorsNoiseRMS[iDet]
              << endmsg;
      if (m_setNoiseOffset) {
        m_noise[detID].second = m_detectorsNoiseOffset[iDet];
        debug() << "Set noise offset for detector " << detName << " (ID=" << detID
                << ") to: " << m_detectorsNoiseOffset[iDet] << endmsg;
      }
    } catch (const std::exception& e) {
      error() << "Error: detector with name " << detName << " not found, exception raised: " << e.what() << endmsg;
      return StatusCode::FAILURE;
    }
  }

  K4RECCALORIMETER_CHECK(AlgTool::initialize());

  return StatusCode::SUCCESS;
}

double ConstNoiseTool::getNoiseRMSPerCell(CellID aCellId) const {

  // Determine sub detector containing the cell and retrieve corresponding noise
  dd4hep::DDSegmentation::CellID cID = aCellId;
  unsigned cellSystem = m_decoder->get(cID, m_systemIndex);
  if (cellSystem >= m_noise.size())
    return 0;
  return m_noise[cellSystem].first;
}

double ConstNoiseTool::getNoiseOffsetPerCell(CellID aCellId) const {

  if (!m_setNoiseOffset) {
    return 0.;
  }

  // Get cells global coordinate "system"
  dd4hep::DDSegmentation::CellID cID = aCellId;
  unsigned cellSystem = m_decoder->get(cID, m_systemIndex);
  if (cellSystem >= m_noise.size())
    return 0;
  return m_noise[cellSystem].second;
}

std::pair<double, double> ConstNoiseTool::getNoisePerCell(CellID aCellId) const {
  // Get cells global coordinate "system"
  dd4hep::DDSegmentation::CellID cID = aCellId;
  unsigned cellSystem = m_decoder->get(cID, m_systemIndex);
  if (cellSystem >= m_noise.size())
    return std::make_pair(0., 0.);
  return m_noise[cellSystem];
}
