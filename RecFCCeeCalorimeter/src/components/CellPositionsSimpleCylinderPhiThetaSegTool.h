#ifndef RECCALORIMETER_CELLPOSITIONSSIMPLECYLINDERPHITHETASEGTOOL_H
#define RECCALORIMETER_CELLPOSITIONSSIMPLECYLINDERPHITHETASEGTOOL_H

// GAUDI
#include "GaudiKernel/AlgTool.h"
// #include "GaudiKernel/ServiceHandle.h"

// FCCSW
// #include "detectorCommon/DetUtils_k4geo.h"
#include "detectorSegmentations/FCCSWGridPhiTheta_k4geo.h"
#include "k4Interface/IGeoSvc.h"
// #include "k4FWCore/DataHandle.h"
#include "k4Interface/ICellPositionsTool.h"

// DD4hep
// #include "DD4hep/Readout.h"
// #include "DD4hep/Volumes.h"
#include "DDSegmentation/Segmentation.h"
// #include "TGeoManager.h"

class IGeoSvc;
namespace DD4hep {
namespace DDSegmentation {
  class Segmentation;
}
} // namespace DD4hep

/** @class CellPositionsSimpleCylinderPhiThetaSegTool
 * k4RecCalorimeter/RecFCCeeCalorimeter/src/components/CellPositionsSimpleCylinderPhiThetaSegTool.h
 * CellPositionsSimpleCylinderPhiThetaSegTool.h
 *
 *  Tool to determine cell positions for a simple cylinder volume with a phi-theta grid readout
 *  Uses layer information stored as CaloData to determine radius at which to position the hits
 *
 *  @author Archil Durglishvili
 *  @author Giovanni Marchiori
 */

class CellPositionsSimpleCylinderPhiThetaSegTool : public extends<AlgTool, ICellPositionsTool> {
public:
  using base_class::base_class;
  ~CellPositionsSimpleCylinderPhiThetaSegTool() = default;

  virtual StatusCode initialize() override final;

  virtual void getPositions(const edm4hep::CalorimeterHitCollection& aCells,
                            edm4hep::CalorimeterHitCollection& outputColl) const override final;

  virtual dd4hep::Position xyzPosition(const uint64_t& aCellId) const override final;

  virtual int layerId(const uint64_t& aCellId) const override final;

private:
  /// Handle to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc { this, "GeoSvc", "GeoSvc" };
  /// Name of the detector
  Gaudi::Property<std::string> m_detectorName{this, "detectorName", "SimpleCylinderDetector"};
  /// Name of the readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "SimpleCylinderPhiTheta"};
  /// Theta-phi segmentation
  dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo* m_segmentation;
  /// Cellid decoder
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  /// Detector radius (for barrel)
  double m_detRadius;
  /// Detector z (for endcap)
  double m_detZ;
  // vector to store layer positions for Barrel and Endcap
  std::vector<double> m_layerPositions;
};
#endif /* RECCALORIMETER_CELLPOSITIONSSIMPLECYLINDERPHITHETASEGTOOL_H */
