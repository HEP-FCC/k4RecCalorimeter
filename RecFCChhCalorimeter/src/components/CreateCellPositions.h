#ifndef DETCOMPONENTS_CREATECELLPOSITIONS_H
#define DETCOMPONENTS_CREATECELLPOSITIONS_H

// FCCSW
#include "k4FWCore/DataHandle.h"
#include "k4Interface/ICellPositionsTool.h"

// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"

#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"

class IGeoSvc;

/** @class CreateCellPositions Reconstruction/RecCalorimeter/src/components/CreateCellPositions.h CreateCellPositions.h
 *
 *  Retrieve positions of the cells from cell ID.
 *  This algorithm saves the centre position of the volume. The segmentation of volume is taken into account.
 *  Transformation matrix from global coordinates to local is taken from DD4hep::Geometry::DetElement.
 *  Full hierarchy of DetElements (for each sensitive volume) is required.
 *
 *  @author Anna Zaborowska, Coralie Neubueser
 *
 */

class CreateCellPositions : public Gaudi::Algorithm {

public:
  CreateCellPositions(const std::string& name, ISvcLocator* svcLoc);
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
  mutable ToolHandle<ICellPositionsTool> m_cellPositionsTool;
  /// Input collection
  mutable DataHandle<edm4hep::CalorimeterHitCollection> m_hits{"hits/hits", Gaudi::DataHandle::Reader, this};
  /// Output collection
  mutable DataHandle<edm4hep::CalorimeterHitCollection> m_positionedHits{"hits/positionedHits", Gaudi::DataHandle::Writer, this};
};

#endif /* DETCOMPONENTS_CREATECELLPOSITIONS_H */
