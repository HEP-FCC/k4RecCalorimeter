#ifndef RECCALORIMETER_CellPositionsHCalBarrelPhiThetaSegTool_H
#define RECCALORIMETER_CellPositionsHCalBarrelPhiThetaSegTool_H

// Gaudi
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/ServiceHandle.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"
#include "detectorSegmentations/FCCSWGridPhiTheta_k4geo.h"

// k4FWCore
#include "k4FWCore/DataHandle.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/ICellPositionsTool.h"

// DD4hep
#include "DD4hep/Readout.h"
#include "DD4hep/Volumes.h"
#include "DDSegmentation/Segmentation.h"

// ROOT
#include "TGeoManager.h"

class IGeoSvc;
namespace DD4hep {
namespace DDSegmentation {
class Segmentation;
}
}

/** @class CellPositionsHCalBarrelPhiThetaSegTool
 *
 */

class CellPositionsHCalBarrelPhiThetaSegTool : public GaudiTool, virtual public ICellPositionsTool {
public:
  CellPositionsHCalBarrelPhiThetaSegTool(const std::string& type, const std::string& name, const IInterface* parent);
  ~CellPositionsHCalBarrelPhiThetaSegTool() = default;

  virtual StatusCode initialize() final;

  virtual StatusCode finalize() final;

  virtual void getPositions(const edm4hep::CalorimeterHitCollection& aCells, edm4hep::CalorimeterHitCollection& outputColl) final;

  virtual dd4hep::Position xyzPosition(const uint64_t& aCellId) const final;

  virtual int layerId(const uint64_t& aCellId) final;

private:
  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;
  /// Name of the electromagnetic calorimeter readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "HCalBarrelReadout"};
  /// segmentation for FCCee cell dimensions: 4x50mm, 6x100mm, 3x200mm
  /// GM: risky to put the numbers here... better to get them from the geometry in memory...
  Gaudi::Property<std::vector<double>> m_radii{this, "radii", {283.55, 288.55, 293.55, 298.55, 306.05, 316.05, 326.05, 336.05, 346.05, 356.05, 371.05, 391.05, 411.05}};
  dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo* m_segmentation;
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  dd4hep::VolumeManager m_volman;
};
#endif /* RECCALORIMETER_CellPositionsHCalBarrelPhiThetaSegTool_H */
