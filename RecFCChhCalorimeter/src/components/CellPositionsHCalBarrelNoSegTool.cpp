#include "CellPositionsHCalBarrelNoSegTool.h"
#include "k4FWCore/GaudiChecks.h"

#include "edm4hep/CalorimeterHitCollection.h"

DECLARE_COMPONENT(CellPositionsHCalBarrelNoSegTool)

StatusCode CellPositionsHCalBarrelNoSegTool::initialize() {
  K4_GAUDI_CHECK(AlgTool::initialize());
  K4_GAUDI_CHECK(m_geoSvc.retrieve());

  // get PhiEta segmentation
  m_segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo*>(
      m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());
  if (m_segmentation == nullptr) {
    error() << "There is no phi-eta segmentation for readout " << m_readoutName << "!!!!" << endmsg;
    // return StatusCode::FAILURE;
  }
  // Take readout bitfield decoder from GeoSvc
  m_decoder = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder();
  m_volman = m_geoSvc->getDetector()->volumeManager();
  // check if decoder contains "layer"
  std::vector<std::string> fields;
  for (uint itField = 0; itField < m_decoder->size(); itField++) {
    fields.push_back((*m_decoder)[itField].name());
  }
  auto iter = std::find(fields.begin(), fields.end(), "layer");
  if (iter == fields.end()) {
    error() << "Readout does not contain field: 'layer'" << endmsg;
  }
  return StatusCode::SUCCESS;
}

void CellPositionsHCalBarrelNoSegTool::getPositions(const edm4hep::CalorimeterHitCollection& aCells,
                                                    edm4hep::CalorimeterHitCollection& outputColl) const {
  debug() << "Input collection size : " << aCells.size() << endmsg;
  // Loop through cell collection
  for (const auto& cell : aCells) {
    auto outPos = CellPositionsHCalBarrelNoSegTool::xyzPosition(cell.getCellID());

    auto edmPos = edm4hep::Vector3f();
    edmPos.x = outPos.x() / dd4hep::mm;
    edmPos.y = outPos.y() / dd4hep::mm;
    edmPos.z = outPos.z() / dd4hep::mm;

    auto positionedHit = cell.clone();
    positionedHit.setPosition(edmPos);
    outputColl.push_back(positionedHit);

    // Debug information about cell position
    debug() << "Cell energy (GeV) : " << positionedHit.getEnergy() << "\tcellID " << positionedHit.getCellID()
            << endmsg;
    debug() << "Position of cell (mm) : \t" << outPos.x() / dd4hep::mm << "\t" << outPos.y() / dd4hep::mm << "\t"
            << outPos.z() / dd4hep::mm << endmsg;
  }
  debug() << "Output positions collection size: " << outputColl.size() << endmsg;
}

dd4hep::Position CellPositionsHCalBarrelNoSegTool::xyzPosition(const CellID aCellId) const {
  // global cartesian coordinates calculated from r,phi,eta, for r=1
  auto detelement = m_volman.lookupDetElement(aCellId);
  double local[] = {0, 0, 0};
  const auto global = detelement.nominal().localToWorld(local);
  double zPos = global.Z();

  dd4hep::DDSegmentation::CellID volumeId = aCellId;
  m_decoder->set(volumeId, "phi", 0);
  m_decoder->set(volumeId, "eta", 0);
  int layer = m_decoder->get(volumeId, "layer");
  double radius = m_radii[layer];

  // x and y calculated with phi position, and radius
  double phi = m_segmentation->phi(aCellId);
  double xPos = cos(phi) * radius;
  double yPos = sin(phi) * radius;

  return dd4hep::Position(xPos, yPos, zPos);
}

int CellPositionsHCalBarrelNoSegTool::layerId(const CellID aCellId) const { return m_decoder->get(aCellId, "layer"); }
