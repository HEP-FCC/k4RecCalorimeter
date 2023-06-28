#include "CellPositionsECalBarrelTool.h"

// EDM
#include "edm4hep/CalorimeterHitCollection.h"

#include <cmath>

DECLARE_COMPONENT(CellPositionsECalBarrelTool)

CellPositionsECalBarrelTool::CellPositionsECalBarrelTool(const std::string& type, const std::string& name,
                                                         const IInterface* parent)
    : GaudiTool(type, name, parent) {
  declareInterface<ICellPositionsTool>(this);
}

StatusCode CellPositionsECalBarrelTool::initialize() {
  StatusCode sc = GaudiTool::initialize();
  if (sc.isFailure()) return sc;
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry service." << endmsg;
    return StatusCode::FAILURE;
  }

  // get segmentation
  auto segmentation = m_geoSvc->lcdd()->readout(m_readoutName).segmentation().segmentation();
  m_segmentationType = -1;

  // determine segmentation type
  // could probably do it using method segmentation::type()...
  m_segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridPhiEta*>(segmentation);
  if (m_segmentation != nullptr) {
    info() << "Found phi-eta segmentation" << endmsg;
    m_segmentationType = 0;
  }
  else {
    warning() << "There is no phi-eta segmentation, trying phi-theta" << endmsg;
    m_segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridPhiTheta*>(segmentation);
    if (m_segmentation != nullptr) {
      info() << "Found phi-theta segmentation" << endmsg;
      m_segmentationType = 1;
    }
    else {
      warning() << "There is no phi-theta segmentation, trying merged module-theta merged" << endmsg;
      m_segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged*>(segmentation);
      if (m_segmentation != nullptr) {
	m_segmentationType = 2;
	info() << "Found merged module-theta segmentation" << endmsg;
	dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged* seg = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged*>(m_segmentation);
	for (int iLayer=0; iLayer<12; iLayer++) {
	  info() << "Layer : " << iLayer << " theta merge : " << seg->mergedThetaCells(iLayer) << " module merge : " << seg->mergedModules(iLayer) << endmsg;
	  if (seg->mergedThetaCells(iLayer)<1) {
	    error() << "Number of cells merged along theta should be >= 1!!!!" << endmsg;
	  }
	  if (seg->mergedModules(iLayer)<1) {
	    error() << "Number of modules merged should be >= 1!!!!" << endmsg;
	  }
	}
      }
      else {
	error() << "There is no module-theta segmentation!!!!" << endmsg;
	return StatusCode::FAILURE;
      }
    }
  }
    

  // Take readout bitfield decoder from GeoSvc
  m_decoder = m_geoSvc->lcdd()->readout(m_readoutName).idSpec().decoder();
  m_volman = m_geoSvc->lcdd()->volumeManager();

  // check if decoder contains "layer"
  std::vector<std::string> fields;
  for (uint itField = 0; itField < m_decoder->size(); itField++) {
    fields.push_back((*m_decoder)[itField].name());
  }
  auto iter = std::find(fields.begin(), fields.end(), "layer");
  if (iter == fields.end()) {
    error() << "Readout does not contain field: 'layer'" << endmsg;
  }

  return sc;
}

void CellPositionsECalBarrelTool::getPositions(const edm4hep::CalorimeterHitCollection& aCells,
                                               edm4hep::CalorimeterHitCollection& outputColl) {

  debug() << "Input collection size : " << aCells.size() << endmsg;

  // Loop through input cell collection, call xyzPosition method for each cell
  // and assign position to cloned hit to be saved in outputColl
  for (const auto& cell : aCells) {
    auto outSeg = CellPositionsECalBarrelTool::xyzPosition(cell.getCellID());
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

dd4hep::Position CellPositionsECalBarrelTool::xyzPosition(const uint64_t& aCellId) const {

  dd4hep::Position outSeg;
  double radius;

  if (m_segmentationType==2) {

    // for module-theta merged segmentation, the local position returned
    // by the segmentation class is theta of group of merged cells 
    // and relative phi wrt first module in group of merged modules
    // in a vector3 scaled such that R_xy=1

    // find position of volume corresponding to first of group of merged cells
    debug() << "cellID: " << aCellId << endmsg;
    dd4hep::DDSegmentation::CellID volumeId = aCellId;
    m_decoder->set(volumeId, "theta", 0);
    debug() << "volumeID: " << volumeId << endmsg;
    auto detelement = m_volman.lookupDetElement(volumeId);
    const auto& transformMatrix = detelement.nominal().worldTransformation();
    double outGlobal[3];
    double inLocal[] = {0, 0, 0};
    transformMatrix.LocalToMaster(inLocal, outGlobal);
    debug() << "Position of volume (mm) : \t" 
	    << outGlobal[0] / dd4hep::mm << "\t" 
	    << outGlobal[1] / dd4hep::mm << "\t"
            << outGlobal[2] / dd4hep::mm << endmsg;

    // get R, phi of volume
    radius = std::sqrt(std::pow(outGlobal[0], 2) + std::pow(outGlobal[1], 2));
    double phi = std::atan2(outGlobal[1], outGlobal[0]);
    debug() << "R (mm), phi of volume : \t" << radius / dd4hep::mm << " , " << phi << endmsg;

    // now get offset in theta and in phi from cell ID (due to theta grid + merging in theta/modules)
    // the local position is normalised to r_xy=1 so theta is atan(1/z)
    dd4hep::DDSegmentation::Vector3D inSeg = m_segmentation->position(aCellId);
    debug() << "Local position of cell (mm) : \t" 
	    << inSeg.x() / dd4hep::mm << "\t" 
	    << inSeg.y() / dd4hep::mm << "\t"
            << inSeg.z() / dd4hep::mm << endmsg;
    double dphi = std::atan2(inSeg.y(), inSeg.x());
    debug() << "Local phi of cell : \t"  << dphi << endmsg;
    phi += dphi;
    outSeg = dd4hep::Position(radius * std::cos(phi),
			      radius * std::sin(phi),
			      radius * inSeg.z());
    /*
    double dtheta = (inSeg.z() == 0.) ? M_PI_2 : std::atan(1./inSeg.z());
    if (inSeg.z() == 0.)
      outSeg = dd4hep::Position(radius * std::cos(phi), 
				radius * std::sin(phi), 
				0.0);
    else
      outSeg = dd4hep::Position(radius * std::cos(phi), 
				radius * std::sin(phi), 
				radius * std::cos(dtheta)/std::sin(dtheta));
    */
  }
  else {
    // for grid eta-phi and theta-phi segmentations, the local position returned
    // by the segmentation class is actually the global position but
    // scaled such that R_xy = 1
    // so first need to find the proper value of R and then to rescale the
    // vector

    // find radius of volume at phi bin=0, eta/theta bin=0
    dd4hep::DDSegmentation::CellID volumeId = aCellId;
    m_decoder->set(volumeId, "phi", 0);
    if (m_segmentationType==0)
      m_decoder->set(volumeId, "eta", 0);
    else if (m_segmentationType==1)
      m_decoder->set(volumeId, "theta", 0);
    auto detelement = m_volman.lookupDetElement(volumeId);
    const auto& transformMatrix = detelement.nominal().worldTransformation();
    double outGlobal[3];
    double inLocal[] = {0, 0, 0};
    transformMatrix.LocalToMaster(inLocal, outGlobal);
    debug() << "Position of volume (mm) : \t" 
	    << outGlobal[0] / dd4hep::mm << "\t" 
	    << outGlobal[1] / dd4hep::mm << "\t"
	    << outGlobal[2] / dd4hep::mm << endmsg;
    radius = std::sqrt(std::pow(outGlobal[0], 2) + std::pow(outGlobal[1], 2));

    // get position of cell from segmentation class (for a radius of 1)
    auto inSeg = m_segmentation->position(aCellId);

    // rescale position by actual radius
    outSeg = dd4hep::Position(inSeg.x() * radius, inSeg.y() * radius, inSeg.z() * radius);
  }

  debug() << "Position of cell (mm) : \t" << outSeg.x() / dd4hep::mm << "\t" << outSeg.y() / dd4hep::mm << "\t"
	  << outSeg.z() / dd4hep::mm << "\n"
	  << endmsg;

  return outSeg;
}

int CellPositionsECalBarrelTool::layerId(const uint64_t& aCellId) {
  int layer;
  dd4hep::DDSegmentation::CellID cID = aCellId;
  layer = m_decoder->get(cID, "layer");
  return layer;
}

StatusCode CellPositionsECalBarrelTool::finalize() { return GaudiTool::finalize(); }
