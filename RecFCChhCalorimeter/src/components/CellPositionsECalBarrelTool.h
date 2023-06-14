#ifndef RECCALORIMETER_CELLPOSITIONSECALBARRELTOOL_H
#define RECCALORIMETER_CELLPOSITIONSECALBARRELTOOL_H

// GAUDI
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/ServiceHandle.h"

// FCCSW
#include "DetCommon/DetUtils.h"
#include "k4Interface/IGeoSvc.h"
#include "DetSegmentation/FCCSWGridPhiEta.h"
#include "DetSegmentation/FCCSWGridPhiTheta.h"
#include "DetSegmentation/FCCSWGridModuleThetaMerged.h"
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

/** @class CellPositionsECalBarrelTool Reconstruction/RecFCChhCalorimeter/src/components/CellPositionsECalBarrelTool.h
 * CellPositionsECalBarrelTool.h
 *
 *  Tool to determine each Calorimeter cell position.
 *
 *   For the FCChh Barrel ECAL, determined from the placed volumes and the FCCSW eta-phi segmentation.   
 * 
 *  @author Coralie Neubueser
 */

class CellPositionsECalBarrelTool : public GaudiTool, virtual public ICellPositionsTool {
public:
  CellPositionsECalBarrelTool(const std::string& type, const std::string& name, const IInterface* parent);
  ~CellPositionsECalBarrelTool() = default;

  virtual StatusCode initialize() final;

  virtual StatusCode finalize() final;

  virtual void getPositions(const edm4hep::CalorimeterHitCollection& aCells, edm4hep::CalorimeterHitCollection& outputColl) final;

  virtual dd4hep::Position xyzPosition(const uint64_t& aCellId) const final;

  virtual int layerId(const uint64_t& aCellId) final;

private:
  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;
  /// Name of the electromagnetic calorimeter readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "ECalBarrelPhiEta"};
  /// Eta-phi / theta-phi / merged module-phi segmentation
  dd4hep::DDSegmentation::Segmentation* m_segmentation;
  /// segmentation type: -1 unknown, 0 eta-phi, 1 theta-phi, 2 theta-module merged
  int m_segmentationType;
  /// Cellid decoder
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  /// Volume manager
  dd4hep::VolumeManager m_volman;
};
#endif /* RECCALORIMETER_CELLPOSITIONSECALBARRELTOOL_H */
