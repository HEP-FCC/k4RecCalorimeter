#ifndef RECCALORIMETER_TOPOCALONEIGHBOURS_H
#define RECCALORIMETER_TOPOCALONEIGHBOURS_H

// from Gaudi
#include "GaudiKernel/AlgTool.h"

// Interfaces
#include "RecCaloCommon/ICaloReadNeighboursMap.h"

class IGeoSvc;

/** @class TopoCaloNeighbours Reconstruction/RecCalorimeter/src/components/TopoCaloNeighbours.h
 *TopoCaloNeighbours.h
 *
 *  Tool that reads a ROOT file containing the TTree with branch "cellId" and branch "neighbours".
 *  This tools reads the tree, creates a map, and allows a lookup of all neighbours of a cell.
 *
 *  @author Anna Zaborowska
 *  @author Coralie Neubueser
 */

class TopoCaloNeighbours : public extends<AlgTool, k4::recCalo::ICaloReadNeighboursMap> {
public:
  using base_class::base_class;

  /** Read a map of cellIDs to vector of cellIDs (neighbours).
   */
  virtual StatusCode initialize() final override;

  /** Function to be called for the neighbours of a cell.
   *   @param[in] aCellId, cellid of the cell of interest.
   *   @return vector of cellIDs, corresponding to the cells neighbours.
   */
  virtual const std::vector<CellID>& neighbours(CellID aCellId) const final override;

private:
  /// Name of input root file that contains the TTree with cellID->vec<neighboursCellID>
  Gaudi::Property<std::string> m_fileName{this, "fileName", "neighbours_map.root"};
  /// Output map to be used for the fast lookup in the topo-clusering algorithm
  std::unordered_map<CellID, std::vector<CellID>> m_map;
};

#endif /* RECCALORIMETER_TOPOCALONEIGHBOURS_H */
