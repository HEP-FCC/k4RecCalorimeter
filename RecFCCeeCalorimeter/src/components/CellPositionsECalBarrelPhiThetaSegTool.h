#ifndef RECCALORIMETER_CELLPOSITIONSECALBARRELPHITHETASEGTOOL_H
#define RECCALORIMETER_CELLPOSITIONSECALBARRELPHITHETASEGTOOL_H

// GAUDI
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ServiceHandle.h"

// FCCSW
#include "detectorCommon/DetUtils_k4geo.h"
#include "detectorSegmentations/FCCSWGridPhiTheta_k4geo.h"
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

/** @class CellPositionsECalBarrelPhiThetaSegTool
 * k4RecCalorimeter/RecFCCeeCalorimeter/src/components/CellPositionsECalBarrelPhiThetaSegTool.h
 * CellPositionsECalBarrelPhiThetaSegTool.h
 *
 *  Tool to determine each Calorimeter cell position.
 *
 *   For the FCCee Barrel ECAL with phi-theta segmentation, determined from the placed volumes and the segmentation.
 *
 *  @author Giovanni Marchiori
 */

class CellPositionsECalBarrelPhiThetaSegTool : public AlgTool, virtual public ICellPositionsTool {
public:
  CellPositionsECalBarrelPhiThetaSegTool(const std::string& type, const std::string& name, const IInterface* parent);
  ~CellPositionsECalBarrelPhiThetaSegTool() = default;

  virtual StatusCode initialize() final;

  virtual StatusCode finalize() final;

  virtual void getPositions(const edm4hep::CalorimeterHitCollection& aCells,
                            edm4hep::CalorimeterHitCollection& outputColl) final;

  virtual dd4hep::Position xyzPosition(const uint64_t& aCellId) const final;

  virtual int layerId(const uint64_t& aCellId) final;

private:
  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;
  /// Name of the electromagnetic calorimeter readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "ECalBarrelPhiTheta"};
  /// Eta-phi segmentation
  dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo* m_segmentation;
  /// Cellid decoder
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  /// Volume manager
  dd4hep::VolumeManager m_volman;
};
#endif /* RECCALORIMETER_CELLPOSITIONSECALBARRELPHITHETASEGTOOL_H */
