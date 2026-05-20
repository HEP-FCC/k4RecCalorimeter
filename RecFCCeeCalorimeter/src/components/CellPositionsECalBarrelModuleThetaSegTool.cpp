#include "CellPositionsECalBarrelModuleThetaSegTool.h"

// EDM
#include "edm4hep/CalorimeterHitCollection.h"
#include "k4FWCore/GaudiChecks.h"

#include <cmath>

DECLARE_COMPONENT(CellPositionsECalBarrelModuleThetaSegTool)

StatusCode CellPositionsECalBarrelModuleThetaSegTool::initialize() {
  K4_GAUDI_CHECK(AlgTool::initialize());
  K4_GAUDI_CHECK(m_geoSvc.retrieve());

  // get segmentation
  m_segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo*>(
      m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());
  if (m_segmentation == nullptr) {
    error() << "There is no module-theta segmentation!!!!" << endmsg;
    return StatusCode::FAILURE;
  }
  debug() << "Found merged module-theta segmentation" << endmsg;
  for (int iLayer = 0; iLayer < m_segmentation->nLayers(); iLayer++) {
    info() << "Layer : " << iLayer << " theta merge : " << m_segmentation->mergedThetaCells(iLayer)
           << " module merge : " << m_segmentation->mergedModules(iLayer) << endmsg;
    if (m_segmentation->mergedThetaCells(iLayer) < 1) {
      error() << "Number of cells merged along theta should be >= 1!!!!" << endmsg;
    }
    if (m_segmentation->mergedModules(iLayer) < 1) {
      error() << "Number of modules merged should be >= 1!!!!" << endmsg;
    }
  }

  m_volman = m_geoSvc->getDetector()->volumeManager();

  return StatusCode::SUCCESS;
}

void CellPositionsECalBarrelModuleThetaSegTool::getPositions(const edm4hep::CalorimeterHitCollection& aCells,
                                                             edm4hep::CalorimeterHitCollection& outputColl) const {

  debug() << "Input collection size : " << aCells.size() << endmsg;

  // Loop through input cell collection, call xyzPosition method for each cell
  // and assign position to cloned hit to be saved in outputColl
  for (const auto& cell : aCells) {
    auto outSeg = CellPositionsECalBarrelModuleThetaSegTool::xyzPosition(cell.getCellID());
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

dd4hep::Position CellPositionsECalBarrelModuleThetaSegTool::xyzPosition(const CellID aCellId) const {

  // find position of volume corresponding to first of group of merged cells
  dd4hep::DDSegmentation::CellID volumeId = m_segmentation->volumeID(aCellId);
  dd4hep::VolumeManagerContext* vc = m_volman.lookupContext(volumeId);
  dd4hep::DDSegmentation::Vector3D inSeg = m_segmentation->position(aCellId);
  dd4hep::Position outSeg = vc->localToWorld(dd4hep::Position(inSeg));
  if (this->msgLevel(MSG::DEBUG)) {
    [[unlikely]] debug() << "cellID: " << aCellId << endmsg;
    debug() << "volumeID: " << volumeId << endmsg;
    debug() << "Local position of cell (mm) : \t" << inSeg.x() / dd4hep::mm << "\t" << inSeg.y() / dd4hep::mm << "\t"
            << inSeg.z() / dd4hep::mm << endmsg;
    debug() << "Position of cell (mm) : \t" << outSeg.x() / dd4hep::mm << "\t" << outSeg.y() / dd4hep::mm << "\t"
            << outSeg.z() / dd4hep::mm << "\n"
            << endmsg;
  }
  return outSeg;
}

int CellPositionsECalBarrelModuleThetaSegTool::layerId(const CellID aCellId) const {
  return m_segmentation->layer(aCellId);
}
