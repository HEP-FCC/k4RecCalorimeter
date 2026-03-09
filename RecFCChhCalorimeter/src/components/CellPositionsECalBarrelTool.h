#ifndef RECCALORIMETER_CELLPOSITIONSECALBARRELTOOL_H
#define RECCALORIMETER_CELLPOSITIONSECALBARRELTOOL_H

// GAUDI
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ServiceHandle.h"

// FCCSW
#include "detectorCommon/DetUtils_k4geo.h"
#include "detectorSegmentations/FCCSWGridPhiEta_k4geo.h"
#include "k4FWCore/DataHandle.h"
#include "k4Interface/ICellPositionsTool.h"
#include "k4Interface/IGeoSvc.h"

// DD4hep
#include "DD4hep/Readout.h"
#include "DD4hep/Volumes.h"
#include "DDSegmentation/Segmentation.h"
#include "TGeoManager.h"

class IGeoSvc;
namespace DD4hep {
namespace DDSegmentation {
  class Segmentation;
}
} // namespace DD4hep

/** @class CellPositionsECalBarrelTool Reconstruction/RecFCChhCalorimeter/src/components/CellPositionsECalBarrelTool.h
 * CellPositionsECalBarrelTool.h
 *
 *  Tool to determine each Calorimeter cell position.
 *
 *   For the FCChh Barrel ECAL, determined from the placed volumes and the FCCSW eta-phi segmentation.
 *
 *  @author Coralie Neubueser
 */

class CellPositionsECalBarrelTool : public extends<AlgTool, ICellPositionsTool> {
public:
  using base_class::base_class;
  ~CellPositionsECalBarrelTool() = default;

  virtual StatusCode initialize() override final;

  virtual void getPositions(const edm4hep::CalorimeterHitCollection& aCells,
                            edm4hep::CalorimeterHitCollection& outputColl) const override final;

  virtual dd4hep::Position xyzPosition(const uint64_t& aCellId) const override final;

  virtual int layerId(const uint64_t& aCellId) const override final;

private:
  /// Handle to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc { this, "GeoSvc", "GeoSvc" };
  /// Name of the electromagnetic calorimeter readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "ECalBarrelPhiEta"};
  /// Eta-phi segmentation
  dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo* m_segmentation;
  /// Cellid decoder
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  /// Volume manager
  dd4hep::VolumeManager m_volman;
};
#endif /* RECCALORIMETER_CELLPOSITIONSECALBARRELTOOL_H */
