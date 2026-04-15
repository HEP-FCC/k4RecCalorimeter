/**
 * @file RecCaloCommon/src/CalorimeterToolBase.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date Aug, 2025
 * @brief Common base for ICalorimeterTool implementations.
 */

#include "RecCaloCommon/CalorimeterToolBase.h"
#include "k4Interface/IGeoSvc.h"
#include "k4FWCore/GaudiChecks.h"
#include "DD4hep/Detector.h"
#include <algorithm>
#include <string>


/** Standard Gaudi initialize method.
 */
StatusCode CalorimeterToolBase::initialize()
{
  K4_GAUDI_CHECK( AlgTool::initialize() );
  K4_GAUDI_CHECK( m_geoSvc.retrieve() );

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

  // Need to do this after m_readout is defined, since it may ask us
  // for our id.
  K4_GAUDI_CHECK( m_constantsSvc.retrieve() );

  K4_GAUDI_CHECK( makeCells() );

  return StatusCode::SUCCESS;
}


/** Return a vector of all existing cells in the current geometry.
 *
 * The result is sorted and unique.
 * Returns an empty vector on error.
 */
auto CalorimeterToolBase::cellIDs() const -> const std::vector<CellID>&
{
  return *m_cells;
}


/** Prepare a map of all existing cells in current geometry.
 *   @param[out] aCells map of existing cells (and deposited energy, set to 0)
 *   return Status code.
 */
StatusCode CalorimeterToolBase::prepareEmptyCells(std::unordered_map<CellID, double>& aCells) const
{
  aCells.clear();
  for (CellID cellId : cellIDs()) {
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


/** Return the name specified for the readout.
 */
const std::string& CalorimeterToolBase::readoutName() const
{
  return m_readoutName;
}


/** Return the subdetector ID.
 */
int CalorimeterToolBase::id() const
{
  if (!m_readout) {
    error() << name() << ": " << "Readout not found; can't find detector ID" << endmsg;
  }
  return m_readout.segmentation().detector()->id;
}


/** Create the list of cells and store with the constants service.
 */
StatusCode CalorimeterToolBase::makeCells()
{
  int detid = id();
  std::string keyName = "cellIDs-";
  if (detid >= 0)
    keyName += std::to_string (detid);
  else
    keyName += "dummy";

  m_cells = m_constantsSvc->getObj<std::vector<uint64_t> > (keyName);
  if (m_cells) {
    return StatusCode::SUCCESS;
  }

  std::vector<uint64_t> cells;
  if (detid >= 0) {
    K4_GAUDI_CHECK( collectCells(cells) );
    // Sort and make unique.
    std::ranges::sort(cells);
    const auto ret = std::ranges::unique(cells);
    cells.erase(ret.begin(), ret.end());
  }

  K4_GAUDI_CHECK( m_constantsSvc->putObj (keyName, std::move(cells)) );
  m_cells = m_constantsSvc->getObj<std::vector<uint64_t> > (keyName);
  K4_GAUDI_CHECK( m_cells != nullptr );

  return StatusCode::SUCCESS;
}
