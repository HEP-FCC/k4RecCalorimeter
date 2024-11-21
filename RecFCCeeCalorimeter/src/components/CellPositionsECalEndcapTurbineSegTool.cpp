#include "CellPositionsECalEndcapTurbineSegTool.h"

// EDM
#include "edm4hep/CalorimeterHitCollection.h"

#include <cmath>

DECLARE_COMPONENT(CellPositionsECalEndcapTurbineSegTool)

CellPositionsECalEndcapTurbineSegTool::CellPositionsECalEndcapTurbineSegTool(const std::string& type, const std::string& name,
										     const IInterface* parent)
    : AlgTool(type, name, parent) {
  declareInterface<ICellPositionsTool>(this);
}

StatusCode CellPositionsECalEndcapTurbineSegTool::initialize() {
  StatusCode sc = AlgTool::initialize();
  if (sc.isFailure()) return sc;
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry service." << endmsg;
    return StatusCode::FAILURE;
  }

  // get segmentation
  m_segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWEndcapTurbine_k4geo*>(m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());
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

  return sc;
}

void CellPositionsECalEndcapTurbineSegTool::getPositions(const edm4hep::CalorimeterHitCollection& aCells,
                                               edm4hep::CalorimeterHitCollection& outputColl) {

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
    debug() << "Cell energy (GeV) : " << positionedHit.getEnergy() << "\tcellID " << positionedHit.getCellID() << endmsg;
    debug() << "Position of cell (mm) : \t" << outSeg.x() / dd4hep::mm << "\t" << outSeg.y() / dd4hep::mm << "\t"
            << outSeg.z() / dd4hep::mm << "\n"
            << endmsg;
  }
  debug() << "Output positions collection size: " << outputColl.size() << endmsg;
}

dd4hep::Position CellPositionsECalEndcapTurbineSegTool::xyzPosition(const uint64_t& aCellId) const {

  dd4hep::Position outSeg;
  double radius;

  // find position of volume corresponding to first of group of merged cells
  debug() << "cellID: " << std::hex << aCellId << std::dec << endmsg;
  dd4hep::DDSegmentation::CellID volumeId = aCellId;
  debug() << "volumeId: " << std::hex << volumeId << std::dec << endmsg;
  //m_decoder->set(volumeId, "rho", 0);
  //debug() << "volumeId: " << volumeId << endmsg;
  m_decoder->set(volumeId, "z", 0);
  //  m_decoder->set(volumeId, "module", 0); 
  debug() << "volumeId: " << std::hex << volumeId << std::dec << endmsg;  
  auto detelement = m_volman.lookupDetElement(volumeId);
  dd4hep::DDSegmentation::CellID iModule = m_decoder->get(aCellId, "module");
  dd4hep::DDSegmentation::CellID iWheel = m_decoder->get(aCellId, "wheel");
  dd4hep::DDSegmentation::CellID iSide = m_decoder->get(aCellId, "side");
  debug() << "Detector wheel, module = " << iWheel << " " << iModule << endmsg;
  debug() << "Detector element placement path " << detelement.placementPath() << endmsg;
  const auto& transformMatrix = detelement.nominal().worldTransformation();
  double rhoHalfWidth = detelement.placement().volume().boundingBox().z();
  double outGlobal[3];
  double inLocal[] = {0, 0, 0};
  transformMatrix.LocalToMaster(inLocal, outGlobal);
  debug() << "Position of volume (mm) : \t" 
	  << outGlobal[0] / dd4hep::mm << "\t" 
	  << outGlobal[1] / dd4hep::mm << "\t"
	  << outGlobal[2] / dd4hep::mm << endmsg;
  
  // get R, phi of volume
  radius = std::sqrt(std::pow(outGlobal[0], 2) + std::pow(outGlobal[1], 2)) + rhoHalfWidth;
  double phi = std::atan2(outGlobal[1], outGlobal[0]);
  debug() << "R (mm), phi of volume : \t" << radius / dd4hep::mm << " , " << phi << endmsg;
  
  dd4hep::DDSegmentation::Vector3D inSeg = m_segmentation->position(aCellId);
  debug() << "Local position of cell (mm) : \t" 
	  << inSeg.x() / dd4hep::mm << "\t" 
	  << inSeg.y() / dd4hep::mm << "\t"
	  << inSeg.z() / dd4hep::mm << endmsg;
  outSeg.SetZ(outGlobal[2]);

  double z = inSeg.z(); // inSeg contains phi, y, z
  double x = inSeg.y();
  //  phi = inSeg.x();
  // this will place rho at the outer edge of each layer, not in the middle.
  // But I think it's the best I can do without importing details of the geometry
  double rho = std::sqrt(std::pow(outGlobal[0], 2) + std::pow(outGlobal[1], 2)) + rhoHalfWidth;
  double y = std::sqrt(rho*rho-x*x);

  //now rotate by phi of the blade, and change from local to global x,y 
  double yprime = x*TMath::Cos(phi) +y*TMath::Sin(phi);
  double xprime = y*TMath::Cos(phi) -x*TMath::Sin(phi);
  
  outSeg.SetX(xprime);
  outSeg.SetY(((int)iSide)*yprime);

  
  debug() << "Position of cell (mm) : \t" << outSeg.x() / dd4hep::mm << "\t" << outSeg.y() / dd4hep::mm << "\t"
	  << outSeg.z() / dd4hep::mm << "\n"
	  << endmsg;

  return outSeg;
}

int CellPositionsECalEndcapTurbineSegTool::layerId(const uint64_t& aCellId) {
  int layer;
  dd4hep::DDSegmentation::CellID cID = aCellId;
  layer = m_decoder->get(cID, "layer");
  return layer;
}

StatusCode CellPositionsECalEndcapTurbineSegTool::finalize() { return AlgTool::finalize(); }
