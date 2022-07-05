#ifndef DETCOMPONENTS_CREATECELLPOSITIONSFCCEE_H
#define DETCOMPONENTS_CREATECELLPOSITIONSFCCEE_H

// FCCSW
#include "k4FWCore/DataHandle.h"
#include "k4Interface/ICellPositionsTool.h"

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"

#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"

#include <unordered_map>

class IGeoSvc;

/** @class CreateCaloCellPositions Reconstruction/RecCalorimeter/src/components/CreateCaloCellPositions.h
 * CreateCaloCellPositions.h
 *
 *  Retrieve positions of the cells from cell ID.
 *  This algorithm saves the centre position of the volume. Defined for all Calo-Subsystems within tools.
 *
 *  @author Coralie Neubueser
 *
 */

class CreateCaloCellPositionsFCCee : public GaudiAlgorithm {

public:
  CreateCaloCellPositionsFCCee(const std::string& name, ISvcLocator* svcLoc);
  /**  Initialize.
   *   @return status code
   */
  StatusCode initialize();
  /**  Execute.
   *   @return status code
   */
  StatusCode execute();
  /**  Finalize.
   *   @return status code
   */
  StatusCode finalize();

private:
  /// Handle for tool to get positions in ECal Barrel
  ToolHandle<ICellPositionsTool> m_cellPositionsECalBarrelTool;
  /// Handle for tool to get positions in HCal Barrel and Ext Barrel, no Segmentation
  ToolHandle<ICellPositionsTool> m_cellPositionsHCalBarrelTool;
  /// Handle for tool to get positions in HCal Barrel and Ext Barrel, no Segmentation
  ToolHandle<ICellPositionsTool> m_cellPositionsHCalExtBarrelTool;
  /// Handle for tool to get positions in Calo Discs
  ToolHandle<ICellPositionsTool> m_cellPositionsEMECTool;
  /// Handle for tool to get positions in Calo Discs
  ToolHandle<ICellPositionsTool> m_cellPositionsHECTool;
  /// Handle for tool to get positions in Calo Discs
  ToolHandle<ICellPositionsTool> m_cellPositionsEMFwdTool;
  /// Handle for tool to get positions in Calo Discs
  ToolHandle<ICellPositionsTool> m_cellPositionsHFwdTool;
  /// Decoder for system ID
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder = new dd4hep::DDSegmentation::BitFieldCoder("system:4");
  /// Input collection
  DataHandle<edm4hep::CalorimeterHitCollection> m_hits{"hits/hits", Gaudi::DataHandle::Reader, this};
  /// Output collection
  DataHandle<edm4hep::CalorimeterHitCollection> m_positionedHits{"hits/positionedHits", Gaudi::DataHandle::Writer, this};
  PodioDataSvc* m_podioDataSvc;
  ServiceHandle<IDataProviderSvc> m_eventDataSvc;

  int m_system_id = m_decoder->index("system");
  std::unordered_map<dd4hep::DDSegmentation::CellID, edm4hep::Vector3f> m_positions_cache{};
};

#endif /* DETCOMPONENTS_CREATECELLPOSITIONSFCCEE_H */
