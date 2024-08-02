#ifndef RECCALORIMETER_CELLPOSITIONSDUMMYTOOL_H
#define RECCALORIMETER_CELLPOSITIONSDUMMYTOOL_H

// GAUDI
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ServiceHandle.h"

// FCCSW
#include "detectorCommon/DetUtils_k4geo.h"
#include "k4Interface/IGeoSvc.h"
#include "k4FWCore/DataHandle.h"
#include "k4Interface/ICellPositionsTool.h"

class IGeoSvc;

/** @class CellPositionsDummyTool Reconstruction/RecFCChhCalorimeter/src/components/CellPositionsDummyTool.h
 *   CellPositionsDummyTool.h
 *
 *  Tool that returns (0,0,0) for each Calorimeter cell.
 * 
 *  @author Coralie Neubueser
*/

class CellPositionsDummyTool : public AlgTool, virtual public ICellPositionsTool {

public:
  CellPositionsDummyTool(const std::string& type, const std::string& name, const IInterface* parent);
  ~CellPositionsDummyTool() = default;

  virtual StatusCode initialize() final;

  virtual StatusCode finalize() final;

  virtual void getPositions(const edm4hep::CalorimeterHitCollection& aCells, edm4hep::CalorimeterHitCollection& outputColl) final;

  virtual dd4hep::Position xyzPosition(const uint64_t& aCellId) const final;

  virtual int layerId(const uint64_t& aCellId) final;

private:
  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;
};
#endif /* RECCALORIMETER_CELLPOSITIONSDUMMYTOOL_H */
