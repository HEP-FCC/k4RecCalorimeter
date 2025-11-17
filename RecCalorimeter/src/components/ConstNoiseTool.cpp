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

  K4RECCALORIMETER_CHECK( m_geoSvc.retrieve() );

  // loop over the detectors
  for (size_t iDet = 0; iDet < m_detectors.size(); iDet++) {
    std::string detName = "DetID_" + m_detectors[iDet];
    try {
      int detID = m_geoSvc->getDetector()->constant<int>(detName);
      m_systemNoiseRMSMap.emplace(detID, m_detectorsNoiseRMS[iDet]);
      debug() << "Set noise RMS for detector " << detName << " (ID=" << detID << ") to: " << m_detectorsNoiseRMS[iDet]
              << endmsg;
      if (m_setNoiseOffset) {
        m_systemNoiseOffsetMap.emplace(detID, m_detectorsNoiseOffset[iDet]);
        debug() << "Set noise offset for detector " << detName << " (ID=" << detID
                << ") to: " << m_detectorsNoiseOffset[iDet] << endmsg;
      }
    } catch (const std::exception& e) {
      error() << "Error: detector with name " << detName << " not found, exception raised: " << e.what() << endmsg;
      return StatusCode::FAILURE;
    }
  }

  K4RECCALORIMETER_CHECK( AlgTool::initialize() );

  return StatusCode::SUCCESS;
}

double ConstNoiseTool::getNoiseRMSPerCell(uint64_t aCellId) const {

  // Get cells global coordinate "system"
  dd4hep::DDSegmentation::CellID cID = aCellId;
  unsigned cellSystem = m_decoder->get(cID, "system");
  // cell noise in system
  auto it = m_systemNoiseRMSMap.find(cellSystem);
  if (it == m_systemNoiseRMSMap.end()) {
    warning() << "No noise RMS set for subsystem " << cellSystem << "! Noise RMS of cell set to 0. " << endmsg;
    return 0;
  }
  return it->second;
}

double ConstNoiseTool::getNoiseOffsetPerCell(uint64_t aCellId) const {

  if (!m_setNoiseOffset) {
    return 0.;
  }

  // Get cells global coordinate "system"
  dd4hep::DDSegmentation::CellID cID = aCellId;
  unsigned cellSystem = m_decoder->get(cID, "system");
  // cell noise in system
  auto it = m_systemNoiseOffsetMap.find(cellSystem);
  if (it == m_systemNoiseOffsetMap.end()) {
    warning() << "No noise offset set for subsystem " << cellSystem << "! Noise offset of cell set to 0. " << endmsg;
    return 0;
  }
  return it->second;
}
