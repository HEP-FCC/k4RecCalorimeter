#ifndef RECCALORIMETER_CELLPOSITIONSECALENDCAPTURBINESEGTOOL_H
#define RECCALORIMETER_CELLPOSITIONSECALENDCAPTURBINESEGTOOL_H

// GAUDI
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ServiceHandle.h"

// k4FWCore
#include "k4FWCore/DataHandle.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/ICellPositionsTool.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"
#include "detectorSegmentations/FCCSWEndcapTurbine_k4geo.h"

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

/** @class CellPositionsECalEndcapTurbineSegTool Reconstruction/RecFCCeeCalorimeter/src/components/CellPositionsECalEndcapTurbineSegTool.h
 * CellPositionsECalEndcapTurbineSegTool.h
 *
 *  Tool to determine each Calorimeter cell position.
 *
 *   For the FCCee Endcap ECAL, determined from the placed volumes and the FCCSW endcap turbine segmentation.   
 * 
 *  @author Erich Varnes
 */

class CellPositionsECalEndcapTurbineSegTool : public AlgTool, virtual public ICellPositionsTool {

public:
  CellPositionsECalEndcapTurbineSegTool(const std::string& type, const std::string& name, const IInterface* parent);
  ~CellPositionsECalEndcapTurbineSegTool() = default;

  virtual StatusCode initialize() final;

  virtual StatusCode finalize() final;

  virtual void getPositions(const edm4hep::CalorimeterHitCollection& aCells, edm4hep::CalorimeterHitCollection& outputColl) final;

  virtual dd4hep::Position xyzPosition(const uint64_t& aCellId) const final;

  virtual int layerId(const uint64_t& aCellId) final;

private:
  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;
  /// Name of the electromagnetic calorimeter readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "ECalEndcapTurbine"};
  /// Merged module-theta segmentation
  dd4hep::DDSegmentation::FCCSWEndcapTurbine_k4geo* m_segmentation;
  /// Cellid decoder
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  /// Volume manager
  dd4hep::VolumeManager m_volman;
};
#endif /* RECCALORIMETER_CELLPOSITIONSECALENDCAPTURBINESEGTOOL_H */
