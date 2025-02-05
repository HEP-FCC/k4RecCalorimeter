#include "CellPositionsECalBarrelModuleThetaSegTool.h"

// EDM
#include "edm4hep/CalorimeterHitCollection.h"

#include <cmath>

DECLARE_COMPONENT(CellPositionsECalBarrelModuleThetaSegTool)

CellPositionsECalBarrelModuleThetaSegTool::CellPositionsECalBarrelModuleThetaSegTool(const std::string& type, const std::string& name,
										     const IInterface* parent)
    : AlgTool(type, name, parent) {
  declareInterface<ICellPositionsTool>(this);
}

StatusCode CellPositionsECalBarrelModuleThetaSegTool::initialize() {
  StatusCode sc = AlgTool::initialize();
  if (sc.isFailure()) return sc;
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry service." << endmsg;
    return StatusCode::FAILURE;
  }

  // get segmentation
  m_segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo*>(m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());
  if (m_segmentation == nullptr) {
    error() << "There is no module-theta segmentation!!!!" << endmsg;
    return StatusCode::FAILURE;
  }
  debug() << "Found merged module-theta segmentation" << endmsg;
  for (int iLayer=0; iLayer<m_segmentation->nLayers(); iLayer++) {
    info() << "Layer : " << iLayer
	   << " theta merge : " << m_segmentation->mergedThetaCells(iLayer)
	   << " module merge : " << m_segmentation->mergedModules(iLayer) << endmsg;
    if (m_segmentation->mergedThetaCells(iLayer)<1) {
      error() << "Number of cells merged along theta should be >= 1!!!!" << endmsg;
    }
    if (m_segmentation->mergedModules(iLayer)<1) {
      error() << "Number of modules merged should be >= 1!!!!" << endmsg;
    }
  }
    

  m_volman = m_geoSvc->getDetector()->volumeManager();

  return sc;
}

void CellPositionsECalBarrelModuleThetaSegTool::getPositions(const edm4hep::CalorimeterHitCollection& aCells,
                                               edm4hep::CalorimeterHitCollection& outputColl) {

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
    debug() << "Cell energy (GeV) : " << positionedHit.getEnergy() << "\tcellID " << positionedHit.getCellID() << endmsg;
    debug() << "Position of cell (mm) : \t" << outSeg.x() / dd4hep::mm << "\t" << outSeg.y() / dd4hep::mm << "\t"
            << outSeg.z() / dd4hep::mm << "\n"
            << endmsg;
  }
  debug() << "Output positions collection size: " << outputColl.size() << endmsg;
}

dd4hep::Position CellPositionsECalBarrelModuleThetaSegTool::xyzPosition(const uint64_t& aCellId) const {

  // find position of volume corresponding to first of group of merged cells
  debug() << "cellID: " << aCellId << endmsg;
  dd4hep::DDSegmentation::CellID volumeId = m_segmentation->volumeID(aCellId);
  debug() << "volumeID: " << volumeId << endmsg;
  dd4hep::VolumeManagerContext* vc = m_volman.lookupContext(volumeId);
  dd4hep::DDSegmentation::Vector3D inSeg = m_segmentation->position(aCellId);
  debug() << "Local position of cell (mm) : \t" 
	  << inSeg.x() / dd4hep::mm << "\t" 
	  << inSeg.y() / dd4hep::mm << "\t"
	  << inSeg.z() / dd4hep::mm << endmsg;
  dd4hep::Position outSeg = vc->localToWorld(dd4hep::Position (inSeg));
  debug() << "Position of cell (mm) : \t" << outSeg.x() / dd4hep::mm << "\t" << outSeg.y() / dd4hep::mm << "\t"
	  << outSeg.z() / dd4hep::mm << "\n"
	  << endmsg;

  return outSeg;
}

int CellPositionsECalBarrelModuleThetaSegTool::layerId(const uint64_t& aCellId) {
  return m_segmentation->layer(aCellId);
}

StatusCode CellPositionsECalBarrelModuleThetaSegTool::finalize() { return AlgTool::finalize(); }
