#ifndef RECINTERFACE_ICALORIMETERTOOL_H
#define RECINTERFACE_ICALORIMETERTOOL_H

#include "RecCaloInterface/ICaloIndexer.h"

// Gaudi
#include "GaudiKernel/IAlgTool.h"

#include <memory>

/** @class ICalorimeterTool RecInterface/RecInterface/ICalorimeterTool.h ICalorimeterTool.h
 *
 *  Abstract interface to calorimeter geometry tool
 *
 *  @author Anna Zaborowska
 */

class ICalorimeterTool : virtual public IAlgTool {
public:
  DeclareInterfaceID(ICalorimeterTool, 1, 0);

  /** Prepare a map of all existing cells in current geometry.
   *   @param[out] aCells map of existing cells (and deposited energy, set to 0)
   *   return Status code.
   */
  virtual StatusCode prepareEmptyCells(std::unordered_map<uint64_t, double>& aCells) = 0;

  /** Return a new indexer object for this subdetector.
   *
   * Returns a null pointer if indexing is not implemented.
   */
  virtual std::unique_ptr<ICaloIndexer> indexer() const { return nullptr; }
};

#endif /* RECINTERFACE_ICALORIMETERTOOL_H */
