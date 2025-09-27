#ifndef RECCALORIMETER_CellPositionsHCalPhiThetaSegTool_H
#define RECCALORIMETER_CellPositionsHCalPhiThetaSegTool_H

// Gaudi
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ServiceHandle.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"
#include "detectorSegmentations/FCCSWGridPhiTheta_k4geo.h"
#include "detectorSegmentations/FCCSWHCalPhiRow_k4geo.h"
#include "detectorSegmentations/FCCSWHCalPhiTheta_k4geo.h"

// k4FWCore
#include "k4FWCore/DataHandle.h"
#include "k4Interface/ICellPositionsTool.h"
#include "k4Interface/IGeoSvc.h"

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
} // namespace DD4hep

/** @class CellPositionsHCalPhiThetaSegTool
 *
 *  Tool to determine each Calorimeter cell position.
 *
 *  For the FCCee HCal Barrel and EndCap with FCCSWGridPhiTheta_k4geo segmentation,
 *  determined from the segmentation and the LayeredCalorimeterData extension.
 *  The LayeredCalorimeterData extension is part of the geometry description.
 *
 *  Positions for FCCSWHCalPhiTheta_k4geo and FCCSWHCalPhiRow_k4geo segmentations are retrieved directly from the
 * segmentation class.
 *
 *  @author Michaela Mlynarikova
 */

class CellPositionsHCalPhiThetaSegTool : public AlgTool, virtual public ICellPositionsTool {
public:
  CellPositionsHCalPhiThetaSegTool(const std::string& type, const std::string& name, const IInterface* parent);
  ~CellPositionsHCalPhiThetaSegTool() = default;

  virtual StatusCode initialize() final;

  virtual StatusCode finalize() final;

  virtual void getPositions(const edm4hep::CalorimeterHitCollection& aCells,
                            edm4hep::CalorimeterHitCollection& outputColl) const final;
  virtual void getPositions(const edm4hep::CalorimeterHitCollection& aCells,
                            edm4hep::CalorimeterHitCollection& outputColl) final
  { const auto* cthis = this;  cthis->getPositions(aCells, outputColl); }

  virtual dd4hep::Position xyzPosition(const uint64_t& aCellId) const final;

  virtual int layerId(const uint64_t& aCellId) const final;
  virtual int layerId(const uint64_t& aCellId) final
  { const auto* cthis = this;  return cthis->layerId(aCellId); }

  virtual std::vector<double> calculateLayerRadii(unsigned int startIndex, unsigned int endIndex);

  virtual std::vector<double> calculateLayerRadiiBarrel();

  virtual std::vector<double> calculateLayerRadiiEndcap();

private:
  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;
  /// Name of the hadronic calorimeter readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "HCalBarrelReadout"};
  /// Name of the hadronic calorimeter
  Gaudi::Property<std::string> m_detectorName{this, "detectorName", "HCalBarrel"};
  /// Segmentation class type
  std::string m_segmentationType;
  /// Theta-phi segmentation
  dd4hep::DDSegmentation::Segmentation* m_segmentation = nullptr;
  /// Cellid decoder
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder = nullptr;
  /// vector to store calculated layer radii
  std::vector<double> m_radii;
  /// layers retrieved from the geometry
  const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>* m_layersRetrieved = nullptr;
  /// for the HCal Endcap, one needs to provide the number of layers in each cylinder
  Gaudi::Property<std::vector<int>> m_numLayersHCalThreeParts{this, "numLayersHCalThreeParts", {6, 9, 22}};
};
#endif /* RECCALORIMETER_CellPositionsHCalPhiThetaSegTool_H */
