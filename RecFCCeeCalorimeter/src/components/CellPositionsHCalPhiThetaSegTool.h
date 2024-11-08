#ifndef RECCALORIMETER_CellPositionsHCalPhiThetaSegTool_H
#define RECCALORIMETER_CellPositionsHCalPhiThetaSegTool_H

// Gaudi
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ServiceHandle.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"
#include "detectorSegmentations/FCCSWHCalPhiTheta_k4geo.h"
#include "detectorSegmentations/FCCSWHCalPhiRow_k4geo.h"

// k4FWCore
#include "k4FWCore/DataHandle.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/ICellPositionsTool.h"

// DD4hep
#include "DD4hep/Readout.h"
#include "DD4hep/Volumes.h"
#include "DDSegmentation/Segmentation.h"
#include <DDRec/DetectorData.h>

// ROOT
#include "TGeoManager.h"

class IGeoSvc;
namespace DD4hep {
namespace DDSegmentation {
class Segmentation;
}
}

/** @class CellPositionsHCalPhiThetaSegTool
 * 
 *  Tool to determine each Calorimeter cell position.
 *
 *  For the FCCee HCal Barrel and EndCap with phi-theta / phi-row segmentation,
 *  retrieved from the FCCSWHCalPhiTheta_k4geo/FCCSWHCalPhiRow_k4geo segmentation.
 *
 *  @author Michaela Mlynarikova
 *  Modified by Archil Durglishvili
 */

class CellPositionsHCalPhiThetaSegTool : public AlgTool, virtual public ICellPositionsTool {
public:
  CellPositionsHCalPhiThetaSegTool(const std::string& type, const std::string& name, const IInterface* parent);
  ~CellPositionsHCalPhiThetaSegTool() = default;

  virtual StatusCode initialize() final;

  virtual StatusCode finalize() final;

  virtual void getPositions(const edm4hep::CalorimeterHitCollection& aCells, edm4hep::CalorimeterHitCollection& outputColl) final;

  virtual dd4hep::Position xyzPosition(const uint64_t& aCellId) const final;

  virtual int layerId(const uint64_t& aCellId) final;

private:
  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;
  /// Name of the hadronic calorimeter readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "HCalBarrelReadout"};
  /// Theta-phi segmentation
  dd4hep::DDSegmentation::Segmentation* m_segmentation = nullptr;
  /// Cellid decoder
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder = nullptr;
};
#endif /* RECCALORIMETER_CellPositionsHCalPhiThetaSegTool_H */
