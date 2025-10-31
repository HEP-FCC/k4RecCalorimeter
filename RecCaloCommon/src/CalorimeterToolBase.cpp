/**
 * @file RecCaloCommon/src/CalorimeterToolBase.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date Aug, 2025
 * @brief Common base for ICalorimeterTool implementations.
 */

#include "RecCaloCommon/CalorimeterToolBase.h"
#include "k4Interface/IGeoSvc.h"
#include "RecCaloCommon/k4RecCalorimeter_check.h"
#include "DD4hep/Detector.h"
#include <algorithm>


/** Standard Gaudi initialize method.
 */
StatusCode CalorimeterToolBase::initialize()
{
  K4RECCALORIMETER_CHECK( AlgTool::initialize() );
  K4RECCALORIMETER_CHECK( m_geoSvc.retrieve() );

  // Look up the readout.
  if (!m_readoutName.empty()) {
    // Check if readouts exist
    info() << "Readout: " << m_readoutName << endmsg;
    const dd4hep::Detector* det = m_geoSvc->getDetector();
    auto it = det->readouts().find(m_readoutName);
    if (it == det->readouts().end()) {
      error() << "Readout <<" << m_readoutName << ">> does not exist." << endmsg;
      return StatusCode::FAILURE;
    }
    m_readout = it->second;
  }

  return StatusCode::SUCCESS;
}


/** Return a vector of all existing cells in the current geometry.
 *
 * The result is sorted and unique.
 * Returns an empty vector on error.
 */
const std::vector<uint64_t>& CalorimeterToolBase::cellIDs() const
{
  {
    std::lock_guard lock (m_mutex);
    if (!m_filledCells) {
      if (collectCells(m_cells).isSuccess()) {
        // Sort and make unique.
        std::ranges::sort(m_cells);
        const auto ret = std::ranges::unique(m_cells);
        m_cells.erase(ret.begin(), ret.end());
      }
      else {
        m_cells.clear();
      }
      m_filledCells = true;
    }
  }
  return m_cells;
}


/** Prepare a map of all existing cells in current geometry.
 *   @param[out] aCells map of existing cells (and deposited energy, set to 0)
 *   return Status code.
 */
StatusCode CalorimeterToolBase::prepareEmptyCells(std::unordered_map<uint64_t, double>& aCells) const
{
  aCells.clear();
  for (uint64_t cellId : cellIDs()) {
    aCells.emplace (cellId, 0);
  }
  if (aCells.empty()) return StatusCode::FAILURE;
  return StatusCode::SUCCESS;
}


/** Return the segmentation assoicated with this geometry.
 */
const dd4hep::DDSegmentation::Segmentation*
CalorimeterToolBase::segmentation() const
{
  if (m_readout.isValid()) {
    return m_readout.segmentation().segmentation();
  }
  return nullptr;
}
