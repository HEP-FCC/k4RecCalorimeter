#ifndef RECCALORIMETER_CELLPOSITIONSHCALBARRELPHISEGTOOL_H
#define RECCALORIMETER_CELLPOSITIONSHCALBARRELPHISEGTOOL_H

// GAUDI
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ServiceHandle.h"

// FCCSW
#include "detectorCommon/DetUtils_k4geo.h"
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

/** @class CellPositionsHCalBarrelPhiSegTool
 * Reconstruction/RecFCChhCalorimeter/src/components/CellPositionsHCalBarrelPhiSegTool.h
 *  CellPositionsHCalBarrelPhiSegTool.h
 *
 *  Tool to determine each Calorimeter cell position.
 *
 *  For the FCChh Barrel and extended Barrel HCAL, determined from the placed volumes.
 *
 *  @author Coralie Neubueser
 */

class CellPositionsHCalBarrelPhiSegTool : public AlgTool, virtual public ICellPositionsTool {
public:
  CellPositionsHCalBarrelPhiSegTool(const std::string& type, const std::string& name, const IInterface* parent);
  ~CellPositionsHCalBarrelPhiSegTool() = default;

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
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "HCalBarrelReadout"};
  /// Cellid decoder
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  /// Volume manager
  dd4hep::VolumeManager m_volman;
  dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo* m_segmentation;
  Gaudi::Property<std::vector<double>> m_radii{
      this, "radii", {291.05, 301.05, 313.55, 328.55, 343.55, 358.55, 378.55, 413.55, 428.55, 453.55}};
};
#endif /* RECCALORIMETER_CELLPOSITIONSHCALBARRELPHISEGTOOL_H */
