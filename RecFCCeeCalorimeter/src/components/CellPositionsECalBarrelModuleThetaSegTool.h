#ifndef RECCALORIMETER_CELLPOSITIONSECALBARRELMODULETHETASEGTOOL_H
#define RECCALORIMETER_CELLPOSITIONSECALBARRELMODULETHETASEGTOOL_H

// GAUDI
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ServiceHandle.h"

// Interfaces
#include "RecCaloCommon/ICellPositionsTool.h"
#include "k4FWCore/DataHandle.h"
#include "k4Interface/IGeoSvc.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"
#include "detectorSegmentations/FCCSWGridModuleThetaMerged_k4geo.h"

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

/** @class CellPositionsECalBarrelModuleThetaSegTool
 * Reconstruction/RecFCCeeCalorimeter/src/components/CellPositionsECalBarrelModuleThetaSegTool.h
 * CellPositionsECalBarrelModuleThetaSegTool.h
 *
 *  Tool to determine each Calorimeter cell position.
 *
 *   For the FCCee Barrel ECAL, determined from the placed volumes and the FCCSW theta-module segmentation.
 *
 *  @author Giovanni Marchiori
 */

class CellPositionsECalBarrelModuleThetaSegTool : public extends<AlgTool, k4::recCalo::ICellPositionsTool> {
public:
  using base_class::base_class;
  virtual ~CellPositionsECalBarrelModuleThetaSegTool() = default;

  virtual StatusCode initialize() override final;

  virtual void getPositions(const edm4hep::CalorimeterHitCollection& aCells,
                            edm4hep::CalorimeterHitCollection& outputColl) const override final;

  virtual dd4hep::Position xyzPosition(const dd4hep::CellID aCellId) const override final;

  virtual int layerId(const dd4hep::CellID aCellId) const override final;

private:
  /// Handle to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc"};
  /// Name of the electromagnetic calorimeter readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "ECalBarrelModuleThetaMerged"};
  /// Merged module-theta segmentation
  dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo* m_segmentation;
  /// Volume manager
  dd4hep::VolumeManager m_volman;
};
#endif /* RECCALORIMETER_CELLPOSITIONSECALBARRELMODULETHETASEGTOOL_H */
