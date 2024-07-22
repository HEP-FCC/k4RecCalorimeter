#include "ConstNoiseTool.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"

// k4FWCore
#include "k4Interface/IGeoSvc.h"

// DD4hep
#include "DD4hep/Detector.h"

DECLARE_COMPONENT(ConstNoiseTool)

ConstNoiseTool::ConstNoiseTool(const std::string& type, const std::string& name, const IInterface* parent)
    : GaudiTool(type, name, parent) {
   declareInterface<INoiseConstTool>(this);
}

StatusCode ConstNoiseTool::initialize() {

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

  // Get GeoSvc
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service. "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }
 
  // loop over the detectors
  for (size_t iDet=0; iDet<m_detectors.size(); iDet++) {
    std::string detName = "DetID_" + m_detectors[iDet];
    try {
      int detID = m_geoSvc->getDetector()->constant<int>(detName);
      m_systemNoiseRMSMap.emplace(detID, m_detectorsNoiseRMS[iDet]);
      debug() << "Set noise RMS for detector " << detName << " (ID=" << detID << ") to: "
	      << m_detectorsNoiseRMS[iDet] << endmsg;
      if (m_setNoiseOffset) {
	m_systemNoiseOffsetMap.emplace(detID, m_detectorsNoiseOffset[iDet]);
	debug() << "Set noise offset for detector " << detName << " (ID=" << detID << ") to: "
	        << m_detectorsNoiseOffset[iDet] << endmsg;
      }
    } catch (const std::exception& e) {
      error() << "Error: detector with name " << detName << " not found, exception raised: " << e.what() << endmsg;
      return StatusCode::FAILURE;
    }
  }

  StatusCode sc = GaudiTool::initialize();
  if (sc.isFailure()) return sc;

  return sc;
}

StatusCode ConstNoiseTool::finalize() {
  StatusCode sc = GaudiTool::finalize();
  return sc;
}

double ConstNoiseTool::getNoiseRMSPerCell(uint64_t aCellId) {

  double noiseRMS = 0.;
  // Get cells global coordinate "system"
  dd4hep::DDSegmentation::CellID cID = aCellId;
  unsigned cellSystem = m_decoder->get(cID, "system");
  // cell noise in system
  if (m_systemNoiseRMSMap.count(cellSystem))
    noiseRMS = m_systemNoiseRMSMap[cellSystem];
  else
    warning() << "No noise RMS set for subsystem " << cellSystem << "! Noise RMS of cell set to 0. " << endmsg;
  return noiseRMS;
}

double ConstNoiseTool::getNoiseOffsetPerCell(uint64_t aCellId) {

  if (!m_setNoiseOffset) {
    return 0.;
  }

  double noiseOffset = 0.;
  // Get cells global coordinate "system"
  dd4hep::DDSegmentation::CellID cID = aCellId;
  unsigned cellSystem = m_decoder->get(cID, "system");
  // cell noise in system
  if (m_systemNoiseOffsetMap.count(cellSystem))
    noiseOffset = m_systemNoiseOffsetMap[cellSystem];
  else
    warning() << "No noise offset set for subsystem " << cellSystem << "! Noise offset of cell set to 0. " << endmsg;
  return noiseOffset;
}
