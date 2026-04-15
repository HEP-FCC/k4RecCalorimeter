// This file's extension implies that it's C, but it's really -*- C++ -*-.
/**
 * @file RecCaloCommon/include/RecCaloCommon/CalorimeterToolBase.h
 * @author scott snyder <snyder@bnl.gov>
 * @date Aug, 2025
 * @brief Common base for ICalorimeterTool implementations.
 */


#ifndef RECCALOCOMMON_CALORIMETERTOOLBASE_H
#define RECCALOCOMMON_CALORIMETERTOOLBASE_H


#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ServiceHandle.h"
#include "RecCaloCommon/ICalorimeterTool.h"
#include "k4Interface/IGeoSvc.h"
#include "DD4hep/Readout.h"
#include <vector>
#include <mutex>


/** @class CalorimeterToolBase RecCaloCommon/include/CalorimeterToolBase.h
 *
 * This factors out the implementations of prepareEmptyCells/cellIDs
 * in terms of a common method collectCells().
 */
class CalorimeterToolBase : public extends<AlgTool, k4::recCalo::ICalorimeterTool>
{
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


private:
  /// Name of the detector readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", ""};

  /// Handle to the geometry service.
  ServiceHandle<IGeoSvc> m_geoSvc { this, "GeoSvc", "GeoSvc" };

  /// Resolved detector readout.
  dd4hep::Readout m_readout;

  // The vector of cells is filled once, the first time we need it.
  mutable std::mutex m_mutex;
  mutable bool m_filledCells = false;
  mutable std::vector<CellID> m_cells;
};


#endif // not RECCALOCOMMON_CALORIMETERTOOLBASE_H
