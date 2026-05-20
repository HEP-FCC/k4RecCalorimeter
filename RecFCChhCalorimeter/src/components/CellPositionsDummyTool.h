#ifndef RECCALORIMETER_CELLPOSITIONSDUMMYTOOL_H
#define RECCALORIMETER_CELLPOSITIONSDUMMYTOOL_H

// GAUDI
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ServiceHandle.h"

// FCCSW
#include "RecCaloCommon/ICellPositionsTool.h"
#include "detectorCommon/DetUtils_k4geo.h"
#include "k4FWCore/DataHandle.h"
#include "k4Interface/IGeoSvc.h"

class IGeoSvc;

/** @class CellPositionsDummyTool Reconstruction/RecFCChhCalorimeter/src/components/CellPositionsDummyTool.h
 *   CellPositionsDummyTool.h
 *
 *  Tool that returns (0,0,0) for each Calorimeter cell.
 *
 *  @author Coralie Neubueser
 */

class CellPositionsDummyTool : public extends<AlgTool, k4::recCalo::ICellPositionsTool> {

public:
  using base_class::base_class;
  ~CellPositionsDummyTool() = default;

  virtual StatusCode initialize() override final;

  virtual void getPositions(const edm4hep::CalorimeterHitCollection& aCells,
                            edm4hep::CalorimeterHitCollection& outputColl) const override final;

  virtual dd4hep::Position xyzPosition(const CellID aCellId) const override final;

  virtual int layerId(const CellID aCellId) const override final;

private:
  /// Handle to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc"};
};
#endif /* RECCALORIMETER_CELLPOSITIONSDUMMYTOOL_H */
