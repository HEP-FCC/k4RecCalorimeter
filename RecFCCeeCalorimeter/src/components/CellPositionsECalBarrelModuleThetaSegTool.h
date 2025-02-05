#ifndef RECCALORIMETER_CELLPOSITIONSECALBARRELMODULETHETASEGTOOL_H
#define RECCALORIMETER_CELLPOSITIONSECALBARRELMODULETHETASEGTOOL_H

// GAUDI
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ServiceHandle.h"

// k4FWCore
#include "k4FWCore/DataHandle.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/ICellPositionsTool.h"

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
}

/** @class CellPositionsECalBarrelModuleThetaSegTool Reconstruction/RecFCCeeCalorimeter/src/components/CellPositionsECalBarrelModuleThetaSegTool.h
 * CellPositionsECalBarrelModuleThetaSegTool.h
 *
 *  Tool to determine each Calorimeter cell position.
 *
 *   For the FCCee Barrel ECAL, determined from the placed volumes and the FCCSW theta-module segmentation.   
 * 
 *  @author Giovanni Marchiori
 */

class CellPositionsECalBarrelModuleThetaSegTool : public AlgTool, virtual public ICellPositionsTool {
public:
  CellPositionsECalBarrelModuleThetaSegTool(const std::string& type, const std::string& name, const IInterface* parent);
  ~CellPositionsECalBarrelModuleThetaSegTool() = default;

  virtual StatusCode initialize() final;

  virtual StatusCode finalize() final;

  virtual void getPositions(const edm4hep::CalorimeterHitCollection& aCells, edm4hep::CalorimeterHitCollection& outputColl) final;

  virtual dd4hep::Position xyzPosition(const uint64_t& aCellId) const final;

  virtual int layerId(const uint64_t& aCellId) final;

private:
  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;
  /// Name of the electromagnetic calorimeter readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "ECalBarrelModuleThetaMerged"};
  /// Merged module-theta segmentation
  dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo* m_segmentation;
  /// Volume manager
  dd4hep::VolumeManager m_volman;
};
#endif /* RECCALORIMETER_CELLPOSITIONSECALBARRELMODULETHETASEGTOOL_H */
