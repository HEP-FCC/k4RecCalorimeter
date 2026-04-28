// This file's extension implies that it's C, but it's really -*- C++ -*-.
/**
 * @file RecCaloCommon/include/RecCaloCommon/CalorimeterToolBase.h
 * @author scott snyder <snyder@bnl.gov>
 * @date Aug, 2025
 * @brief Common base for ICalorimeterTool implementations.
 */

#ifndef RECCALOCOMMON_CALORIMETERTOOLBASE_H
#define RECCALOCOMMON_CALORIMETERTOOLBASE_H

#include "DD4hep/Readout.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ServiceHandle.h"
#include "RecCaloCommon/ICaloCellConstantsSvc.h"
#include "RecCaloCommon/ICalorimeterTool.h"
#include "k4Interface/IGeoSvc.h"
#include <vector>

/** @class CalorimeterToolBase RecCaloCommon/include/CalorimeterToolBase.h
 *
 * This factors out the implementations of prepareEmptyCells/cellIDs
 * in terms of a common method collectCells().
 */
class CalorimeterToolBase : public extends<AlgTool, k4::recCalo::ICalorimeterTool> {
public:
  using base_class::base_class;

  /** Standard Gaudi initialize method.
   */
  virtual StatusCode initialize() override;

  /** Return a vector of all existing cells in the current geometry.
   *
   * The result is sorted and unique.
   * Returns an empty vector on error.
   */
  virtual const std::vector<CellID>& cellIDs() const final override;

  /** Prepare a map of all existing cells in current geometry.
   *   @param[out] aCells map of existing cells (and deposited energy, set to 0)
   *   return Status code.
   */
  virtual StatusCode prepareEmptyCells(std::unordered_map<CellID, double>& aCells) const final override;

  /** Return the segmentation associated with this geometry.
   */
  virtual const dd4hep::DDSegmentation::Segmentation* segmentation() const final override;

  /** Return the name specified for the readout.
   */
  virtual const std::string& readoutName() const final override;

  /** Return the subdetector ID.
   */
  virtual int id() const final override;

protected:
  /// Return the resolved readout.
  const dd4hep::Readout readout() const { return m_readout; }

  /** Fill vector with all existing cells for this geometry.
   */
  virtual StatusCode collectCells(std::vector<CellID>& cells) const = 0;

  /// Create the list of cells and store with the constants service.
  StatusCode makeCells();

private:
  /// Name of the detector readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", ""};

  /// Handle to the geometry service.
  ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc"};

  /// Handle to the cell constants service.
  /// Used to hold the set of cell IDs.
  ServiceHandle<k4::recCalo::ICaloCellConstantsSvc> m_constantsSvc{this, "CaloCellConstantsSvc",
                                                                   "k4::recCalo::CaloCellConstantsSvc", ""};

  /// Resolved detector readout.
  dd4hep::Readout m_readout;

  // Pointer to the vector of cells.  The vector itself is stored in the
  // constants service; we create it if it's not already there.
  const std::vector<CellID>* m_cells;
};

#endif // not RECCALOCOMMON_CALORIMETERTOOLBASE_H
