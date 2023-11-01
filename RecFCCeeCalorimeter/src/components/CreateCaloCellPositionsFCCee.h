#ifndef DETCOMPONENTS_CREATECELLPOSITIONSFCCEE_H
#define DETCOMPONENTS_CREATECELLPOSITIONSFCCEE_H

// FCCSW
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/MetaDataHandle.h"
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
  /// Handle for tool to get positions width
  ToolHandle<ICellPositionsTool> m_cellPositionsTool{};
  /// Input collection
  DataHandle<edm4hep::CalorimeterHitCollection> m_hits{"hits/hits", Gaudi::DataHandle::Reader, this};
  /// Input collection metadata handle. The name MUST be "CellIDEncoding"
  MetaDataHandle<std::string> m_hitsCellIDEncoding{m_hits, "CellIDEncoding", Gaudi::DataHandle::Reader};
  /// Output collection
  DataHandle<edm4hep::CalorimeterHitCollection> m_positionedHits{"hits/positionedHits", Gaudi::DataHandle::Writer, this};
  /// Output collection metadata handle. The name MUST be "CellIDEncoding"
  MetaDataHandle<std::string> m_positionedHitsCellIDEncoding{m_positionedHits, "CellIDEncoding", Gaudi::DataHandle::Writer};

  // Cache
  std::unordered_map<dd4hep::DDSegmentation::CellID, edm4hep::Vector3f> m_positions_cache{};
};

#endif /* DETCOMPONENTS_CREATECELLPOSITIONSFCCEE_H */
