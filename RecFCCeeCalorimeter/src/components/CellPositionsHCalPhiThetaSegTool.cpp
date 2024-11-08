#include "CellPositionsHCalPhiThetaSegTool.h"

#include "edm4hep/CalorimeterHitCollection.h"
#include <DDRec/DetectorData.h>

using dd4hep::DetElement;

DECLARE_COMPONENT(CellPositionsHCalPhiThetaSegTool)

CellPositionsHCalPhiThetaSegTool::CellPositionsHCalPhiThetaSegTool(
    const std::string &type,
    const std::string &name,
    const IInterface *parent)
    : AlgTool(type, name, parent)
{
  declareInterface<ICellPositionsTool>(this);
}

StatusCode CellPositionsHCalPhiThetaSegTool::initialize()
{
  StatusCode sc = AlgTool::initialize();
  if (sc.isFailure())
    return sc;
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc)
  {
    error() << "Unable to locate Geometry service." << endmsg;
    return StatusCode::FAILURE;
  }

  std::string segmentationType = m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation()->type();
  if(segmentationType == "FCCSWHCalPhiTheta_k4geo")
  {
    // get PhiTheta segmentation
    m_segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWHCalPhiTheta_k4geo *>(
        m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());
  }
  if(segmentationType == "FCCSWHCalPhiRow_k4geo")
  {
    // get PhiRow segmentation
    m_segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWHCalPhiRow_k4geo *>(
        m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());
  }

  if (m_segmentation == nullptr)
  {
    error() << "There is no phi-theta segmentation!!!!" << endmsg;
    return StatusCode::FAILURE;
  }
  // Take readout bitfield decoder from GeoSvc
  m_decoder = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder();

  return sc;
}

void CellPositionsHCalPhiThetaSegTool::getPositions(const edm4hep::CalorimeterHitCollection &aCells,
                                                          edm4hep::CalorimeterHitCollection &outputColl)
{
  debug() << "Input collection size : " << aCells.size() << endmsg;
  // Loop through cell collection
  for (const auto &cell : aCells)
  {
    auto outSeg = CellPositionsHCalPhiThetaSegTool::xyzPosition(cell.getCellID());

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
            << outSeg.z() / dd4hep::mm << endmsg;
  }
  debug() << "Output positions collection size: " << outputColl.size() << endmsg;
}


dd4hep::Position CellPositionsHCalPhiThetaSegTool::xyzPosition(const uint64_t &aCellId) const
{
  // get global position
  auto posCell = m_segmentation->position(aCellId);
  dd4hep::Position outSeg(posCell.x(), posCell.y(), posCell.z());

  debug() << "Layer : " << m_decoder->get(aCellId, "layer") << endmsg;
  debug() << "Global position : x = " << outSeg.x() << " y = "  << outSeg.y() << " z = "  << outSeg.z() << endmsg;

  return outSeg;
}


int CellPositionsHCalPhiThetaSegTool::layerId(const uint64_t &aCellId)
{
  int layer;
  layer = m_decoder->get(aCellId, "layer");
  return layer;
}
StatusCode CellPositionsHCalPhiThetaSegTool::finalize() { return AlgTool::finalize(); }
