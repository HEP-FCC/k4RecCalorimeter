#ifndef RECCALORIMETER_CELLPOSITIONSHCALBARRELTOOL_H
#define RECCALORIMETER_CELLPOSITIONSHCALBARRELTOOL_H

// GAUDI
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ServiceHandle.h"

// FCCSW
#include "detectorCommon/DetUtils_k4geo.h"
#include "k4Interface/IGeoSvc.h"
#include "detectorSegmentations/FCCSWGridPhiEta_k4geo.h"
#include "k4FWCore/DataHandle.h"
#include "k4Interface/ICellPositionsTool.h"

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

/** @class CellPositionsHCalBarrelTool
 *
 */

class CellPositionsHCalBarrelTool : public AlgTool, virtual public ICellPositionsTool {
public:
  CellPositionsHCalBarrelTool(const std::string& type, const std::string& name, const IInterface* parent);
  ~CellPositionsHCalBarrelTool() = default;

  virtual StatusCode initialize() final;

  virtual StatusCode finalize() final;

  virtual void getPositions(const edm4hep::CalorimeterHitCollection& aCells, edm4hep::CalorimeterHitCollection& outputColl) final;

  virtual dd4hep::Position xyzPosition(const uint64_t& aCellId) const final;

  virtual int layerId(const uint64_t& aCellId) final;

private:
  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;
  /// Name of the electromagnetic calorimeter readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "BarHCal_Readout_phieta"};
  /// old segmentation used for FCChh cell dimensions
  /// Gaudi::Property<std::vector<double>> m_radii{this, "radii", {291.05, 301.05, 313.55, 328.55, 343.55, 358.55, 378.55, 413.55, 428.55, 453.55}};
  /// segmentation for FCCee cell dimensions: 4x50mm, 6x100mm, 3x200mm
  Gaudi::Property<std::vector<double>> m_radii{this, "radii", {283.55, 288.55, 293.55, 298.55, 306.05, 316.05, 326.05, 336.05, 346.05, 356.05, 371.05, 391.05, 411.05}};
  dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo* m_segmentation;
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  dd4hep::VolumeManager m_volman;
};
#endif /* RECCALORIMETER_CELLPOSITIONSHCALBARRELTOOL_H */
