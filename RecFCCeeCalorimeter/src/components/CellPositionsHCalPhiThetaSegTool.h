#ifndef RECCALORIMETER_CellPositionsHCalPhiThetaSegTool_H
#define RECCALORIMETER_CellPositionsHCalPhiThetaSegTool_H

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
 */

class CellPositionsHCalPhiThetaSegTool : public GaudiTool, virtual public ICellPositionsTool {
public:
  CellPositionsHCalPhiThetaSegTool(const std::string& type, const std::string& name, const IInterface* parent);
  ~CellPositionsHCalPhiThetaSegTool() = default;

  virtual StatusCode initialize() final;

  virtual StatusCode finalize() final;

  virtual void getPositions(const edm4hep::CalorimeterHitCollection& aCells, edm4hep::CalorimeterHitCollection& outputColl) final;

  virtual dd4hep::Position xyzPosition(const uint64_t& aCellId) const final;

  virtual int layerId(const uint64_t& aCellId) final;

  virtual std::vector<double> calculateLayerRadii(unsigned int startIndex, unsigned int endIndex); 
  virtual std::vector<double> calculateLayerRadiiBarrel(); 
  virtual std::vector<double> calculateLayerRadiiEndcap(); 

private:
  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;
  /// Name of the calorimeter readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "HCalBarrelReadout"};
  /// Name of the calorimeter 
  Gaudi::Property<std::string> m_detectorName{this, "detectorName", "HCalBarrel"};

  std::vector<double> m_radii;
  dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo* m_segmentation;
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  dd4hep::VolumeManager m_volman;
  // layer radii calculated on the flight from the geometry 
  const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>* m_layersRetrieved;
  Gaudi::Property<std::vector<int>> m_numLayersHCalThreeParts{this, "numLayersHCalThreeParts", {6,9,22}};
};
#endif /* RECCALORIMETER_CellPositionsHCalPhiThetaSegTool_H */