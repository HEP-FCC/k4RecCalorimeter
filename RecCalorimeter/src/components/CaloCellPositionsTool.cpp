/**
 * @file k4RecCalorimeter/RecCalorimeter/src/components/CaloCellPositionsTool.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date Apr, 2025
 * @brief Generic tool to find positions of calorimeter cells.
 */


#include "CaloCellPositionsTool.h"
#include "k4FWCore/GaudiChecks.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "DD4hep/Detector.h"
#include "DD4hep/Readout.h"


DECLARE_COMPONENT(CaloCellPositionsTool)


/** Standard Gaudi initialize method.
 */
StatusCode CaloCellPositionsTool::initialize()
{
  K4_GAUDI_CHECK( base_class::initialize() );
  K4_GAUDI_CHECK( m_geoSvc.retrieve() );

  dd4hep::Readout readout = m_geoSvc->getDetector()->readout(m_readoutName);
  m_segmentation = readout.segmentation();
  m_decoder = m_segmentation.decoder();
  m_layerFieldIdx = m_decoder->index(m_layerFieldName);

  const dd4hep::DetElementObject& de = m_segmentation.detector();
  dd4hep::VolumeManager vman_glob = m_geoSvc->getDetector()->volumeManager();
  m_volman = vman_glob.subdetector (de.id);

  return StatusCode::SUCCESS;
}


/** Return the cartesian global coordinates of a cell center
 */
dd4hep::Position
CaloCellPositionsTool::xyzPosition(const uint64_t& aCellId) const
{
  // Find the volume corresponding to this cell.
  dd4hep::DDSegmentation::CellID volumeId = m_segmentation->volumeID(aCellId);
  dd4hep::VolumeManagerContext* vc = m_volman.lookupContext(volumeId);

  // Find the position in the local coordinates of that volume.
  dd4hep::DDSegmentation::Vector3D inSeg = m_segmentation->position(aCellId);

  // Transform to global coordinates.
  return vc->localToWorld(dd4hep::Position (inSeg));
}


/** Copy cells from aCells to outputColl filling in cell positions.
 */
void CaloCellPositionsTool::getPositions(const edm4hep::CalorimeterHitCollection& aCells,
                                         edm4hep::CalorimeterHitCollection& outputColl) const
{
  // Loop through input cell collection, call xyzPosition method for each cell
  // and assign position to cloned hit to be saved in outputColl
  for (const edm4hep::CalorimeterHit& cell : aCells) {
    static constexpr double inv_mm = 1 / dd4hep::mm;
    dd4hep::Position pos = this->xyzPosition(cell.getCellID()) * inv_mm;

    edm4hep::MutableCalorimeterHit positionedHit = cell.clone();
    positionedHit.setPosition({static_cast<float>(pos.x()),
                               static_cast<float>(pos.y()),
                               static_cast<float>(pos.z())});
    outputColl.push_back(positionedHit);
  }
}


/** Return the layer number of a cell.
 */
int CaloCellPositionsTool::layerId(const uint64_t& aCellId) const
{
  return m_decoder->get(aCellId, m_layerFieldIdx);
}
