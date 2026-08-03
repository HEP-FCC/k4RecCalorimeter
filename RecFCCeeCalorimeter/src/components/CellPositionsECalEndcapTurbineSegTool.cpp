#include "CellPositionsECalEndcapTurbineSegTool.h"
#include "k4FWCore/GaudiChecks.h"

// EDM
#include "edm4hep/CalorimeterHitCollection.h"

#include <cmath>

DECLARE_COMPONENT(CellPositionsECalEndcapTurbineSegTool)

StatusCode CellPositionsECalEndcapTurbineSegTool::initialize() {
  K4_GAUDI_CHECK(AlgTool::initialize());
  K4_GAUDI_CHECK(m_geoSvc.retrieve());

  // get segmentation
  m_segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWEndcapTurbine_k4geo*>(
      m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());
  if (m_segmentation == nullptr) {
    error() << "There is no endcap turbine segmentation!!!!" << endmsg;
    return StatusCode::FAILURE;
  }
  debug() << "Found endcap turbine segmentation" << endmsg;

  // Take readout bitfield decoder from GeoSvc
  m_decoder = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder();
  m_volman = m_geoSvc->getDetector()->volumeManager();

  // check if decoder contains "layer"
  std::vector<std::string> fields;
  for (uint itField = 0; itField < m_decoder->size(); itField++) {
    fields.push_back((*m_decoder)[itField].name());
    debug() << "In positioning tool, field is " << (*m_decoder)[itField].name() << endmsg;
  }
  auto iter = std::find(fields.begin(), fields.end(), "layer");
  if (iter == fields.end()) {
    error() << "Readout does not contain field: 'layer'" << endmsg;
  }

  return StatusCode::SUCCESS;
}

void CellPositionsECalEndcapTurbineSegTool::getPositions(const edm4hep::CalorimeterHitCollection& aCells,
                                                         edm4hep::CalorimeterHitCollection& outputColl) const {

  debug() << "Input collection size : " << aCells.size() << endmsg;

  // Loop through input cell collection, call xyzPosition method for each cell
  // and assign position to cloned hit to be saved in outputColl
  for (const auto& cell : aCells) {
    auto outSeg = CellPositionsECalEndcapTurbineSegTool::xyzPosition(cell.getCellID());
    auto edmPos = edm4hep::Vector3f();
    edmPos.x = outSeg.x() / dd4hep::mm;
    edmPos.y = outSeg.y() / dd4hep::mm;
    edmPos.z = outSeg.z() / dd4hep::mm;

    auto positionedHit = cell.clone();
    positionedHit.setPosition(edmPos);
    outputColl.push_back(positionedHit);

    // Debug information about cell position
    debug() << "Cell energy (GeV) : " << positionedHit.getEnergy() << "\tcellID " << positionedHit.getCellID()
            << endmsg;
    debug() << "Position of cell (mm) : \t" << outSeg.x() / dd4hep::mm << "\t" << outSeg.y() / dd4hep::mm << "\t"
            << outSeg.z() / dd4hep::mm << "\n"
            << endmsg;
  }
  debug() << "Output positions collection size: " << outputColl.size() << endmsg;
}

dd4hep::Position CellPositionsECalEndcapTurbineSegTool::xyzPosition(const CellID aCellId) const {

  // find position of volume corresponding to first of group of merged cells
  debug() << "cellID: " << aCellId << endmsg;

  auto detelement = m_volman.lookupDetElement(aCellId);

  // get local position of readout cell wrt parent calibration layer
  dd4hep::DDSegmentation::Vector3D inSeg = m_segmentation->position(aCellId);
  debug() << "Local position of cell (mm) : \t" << inSeg.x() / dd4hep::mm << "\t" << inSeg.y() / dd4hep::mm << "\t"
          << inSeg.z() / dd4hep::mm << endmsg;
  // translate to global position
  dd4hep::Position outSeg = detelement.nominal().localToWorld(dd4hep::Position(inSeg));
  CellID iSide = m_decoder->get(aCellId, "side");
  if (iSide != 1) {
    // account for the fact that -z endcap is mirrored from the +z one
    outSeg.SetZ(-outSeg.z());
    outSeg.SetY(-outSeg.y());
  }

  debug() << "Position of cell (mm) : \t" << outSeg.x() / dd4hep::mm << "\t" << outSeg.y() / dd4hep::mm << "\t"
          << outSeg.z() / dd4hep::mm << "\n"
          << endmsg;

  return outSeg;
}

int CellPositionsECalEndcapTurbineSegTool::layerId(const CellID aCellId) const {
  return m_decoder->get(aCellId, "layer");
}
