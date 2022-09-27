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
  /// ECal Barrel System ID
  Gaudi::Property<size_t> m_systemIdECalBarrel {
      this, "systemIdECalBarrel", 4, "System ID of ECal Barrel"
  };
  /// HCal Barrel System ID
  Gaudi::Property<size_t> m_systemIdHCalBarrel {
      this, "systemIdHCalBarrel", 10, "System ID of HCal Barrel"
  };
  /// HCal Extended Barrel System ID
  Gaudi::Property<size_t> m_systemIdHCalExtBarrel {
      this, "systemIdHCalExtBarrel", 9, "System ID of HCal Extended Barrel"
  };
  /// EMEC System ID
  Gaudi::Property<size_t> m_systemIdEMEC {
      this, "systemIdEMEC", 6, "System ID of EMEC"
  };
  /// HEC System ID
  Gaudi::Property<size_t> m_systemIdHEC {
      this, "systemIdHEC", 7, "System ID of HEC"
  };
  /// EMFwd System ID
  Gaudi::Property<size_t> m_systemIdEMFwd {
      this, "systemIdEMFwd", 999, "System ID of EMFwd"
  };
  /// HFwd System ID
  Gaudi::Property<size_t> m_systemIdHFwd {
      this, "systemIdHFwd", 999, "System ID of HFwd"
  };

  /// Handle for tool to get positions in ECal Barrel
  ToolHandle<ICellPositionsTool> m_cellPositionsECalBarrelTool{"CellPositionsECalBarrelTool/CellPositionsTool"};
  /// Handle for tool to get positions in HCal Barrel and Ext Barrel, no Segmentation
  ToolHandle<ICellPositionsTool> m_cellPositionsHCalBarrelTool{"CellPositionsHCalBarrelTool/CellPositionsTool"};
  /// Handle for tool to get positions in HCal Barrel and Ext Barrel, no Segmentation
  ToolHandle<ICellPositionsTool> m_cellPositionsHCalExtBarrelTool{"CellPositionsHCalBarrelNoSegTool/CellPositionsTool"};
  /// Handle for tool to get positions in Calo Discs
  ToolHandle<ICellPositionsTool> m_cellPositionsEMECTool{"CellPositionsCaloDiscsTool/CellPositionsTool"};
  /// Handle for tool to get positions in Calo Discs
  ToolHandle<ICellPositionsTool> m_cellPositionsHECTool{"CellPositionsCaloDiscsTool/CellPositionsTool"};
  /// Handle for tool to get positions in Calo Discs
  ToolHandle<ICellPositionsTool> m_cellPositionsEMFwdTool{"CellPositionsCaloDiscsTool/CellPositionsTool"};
  /// Handle for tool to get positions in Calo Discs
  ToolHandle<ICellPositionsTool> m_cellPositionsHFwdTool{"CellPositionsCaloDiscsTool/CellPositionsTool"};
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
