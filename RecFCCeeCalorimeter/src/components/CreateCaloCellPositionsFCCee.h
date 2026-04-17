#ifndef DETCOMPONENTS_CREATECELLPOSITIONSFCCEE_H
#define DETCOMPONENTS_CREATECELLPOSITIONSFCCEE_H

// k4FWCore
#include "k4FWCore/DataHandle.h"
#include "RecCaloCommon/ICellPositionsTool.h"

// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"

// edm4hep
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

class CreateCaloCellPositionsFCCee : public Gaudi::Algorithm {

public:
  CreateCaloCellPositionsFCCee(const std::string& name, ISvcLocator* svcLoc);
  /**  Initialize.
   *   @return status code
   */
  StatusCode initialize();
  /**  Execute.
   *   @return status code
   */
  StatusCode execute(const EventContext&) const;
  /**  Finalize.
   *   @return status code
   */
  StatusCode finalize();

private:
  /// Handle for tool to get positions
  ToolHandle<k4::recCalo::ICellPositionsTool> m_cellPositionsTool{};
  /// Input collection
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_hits{"hits/hits", Gaudi::DataHandle::Reader, this};
  /// Output collection
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_positionedHits{"hits/positionedHits",
                                                                                   Gaudi::DataHandle::Writer, this};

  // Cache
  mutable std::unordered_map<dd4hep::DDSegmentation::CellID, edm4hep::Vector3f> m_positions_cache{};
};

#endif /* DETCOMPONENTS_CREATECELLPOSITIONSFCCEE_H */
