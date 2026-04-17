/**
 * @file RecCaloCommon/src/MultiIndexer.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date Feb, 2026
 * @brief Implementation of ICaloIndexer for multiple detectors.
 */


#include "RecCaloCommon/MultiIndexer.h"
#include <string>
#include <format>


namespace k4::recCalo {


/**
 * @brief For reporting errors from the constructor.
 */
MultiIndexer::MultiIndexerException::MultiIndexerException (const std::string& what)
  : std::runtime_error ("MultiIndexerException: " + what)
{
}


/**
 * @brief Dummy indexer that always returns INVALID.
 */
class DummyIndexer : public ICaloIndexer
{
  virtual index_t index (uint64_t /*id*/) const override final
  { return ICaloIndexer::INVALID; }
  virtual std::span<const uint64_t> cellIDs() const override final
  { return std::span<const uint64_t>(); }
  virtual std::span<const int> detIDs() const override final
  { return std::span<const int>(); }
  virtual size_t detIDBits() const override final { return 4; }
};


/**
 * @brief Constructor.
 */
MultiIndexer::MultiIndexer (size_t detIDBits,
                            std::span<const ICaloIndexer* const> indexers,
                            ICaloCellConstantsSvc& constsSvc)
  : m_detIDBits (detIDBits)
{
  // Check that the number of detector ID bits is reasonable.
  if (detIDBits < 1 || detIDBits > 8) {
    throw MultiIndexerException (std::format ("detIDBits {} out of range; should be between 1 and 8",
                                              detIDBits));
  }

  // Mask for extracting detector ID.
  m_detIDMask = (static_cast<uint64_t>(1) << detIDBits) - 1;

  // Resize the indexer vector to the number of possible detector IDs,
  // and fill it with pointers to a dummy indexer that always
  // returns INVALID.
  static const DummyIndexer dummyIndexer;
  m_indexers.resize (m_detIDMask+1, std::make_pair (&dummyIndexer, 0));

  std::string keyName = "cellIDs";
  size_t totcells = 0;

  // Loop over indexers.
  for (const ICaloIndexer* indexer : indexers) {
    // Check that the detector ID from this indexer is OK.
    std::span<const int> detIDs = indexer->detIDs();
    if (detIDs.size() != 1) {
      throw MultiIndexerException (std::format ("bad indexer returns detIDs of size {}",
                                                detIDs.size()));
    }
    if (detIDs[0] > static_cast<int>(m_detIDMask)) {
      throw MultiIndexerException (std::format ("bad indexer returns out-of-range detID {}",
                                                detIDs[0]));
    }

    // Remember this ID.
    m_detIDs.push_back (detIDs[0]);

    // Fill in the indexer and offset.
    m_indexers[detIDs[0]].first = indexer;
    m_indexers[detIDs[0]].second = totcells;

    // Append the ID to the key, and count the total number of cells.
    keyName += "-" + std::to_string(detIDs[0]);
    totcells += indexer->cellIDs().size();
  }

  // Now we get the list of cellIDs.
  // Don't need to do anything if it would be empty.
  if (totcells == 0) return;

  // See if we've already made this combination.
  const std::vector<uint64_t>* cellsptr =
    constsSvc.getObj<std::vector<uint64_t> > (keyName);
  if (!cellsptr) {
    // Nope.  Concatenate the cell ID lists from each indexer.
    std::vector<uint64_t> cells;
    cells.reserve (totcells);
    for (const ICaloIndexer* indexer : indexers) {
#if __cpp_lib_containers_ranges // c++23
      cells.append_range (indexer->cellIDs());
#else
      std::span<const uint64_t> ids = indexer->cellIDs();
      cells.insert (cells.end(), ids.begin(), ids.end());
#endif
    }

    // And store in the constants service.
    constsSvc.putObj (keyName, std::move (cells));
    cellsptr = constsSvc.getObj<std::vector<uint64_t> > (keyName);
  }

  // Some validity checking.
  if (!cellsptr) {
    throw MultiIndexerException (std::format ("cannot get cells vector"));
  }
  if (cellsptr->size() != totcells) {
    throw MultiIndexerException (std::format ("cells vector size mismatch.  Got {}; expected {}",
                                              cellsptr->size(), totcells));
  }

  // Save the span over cells.
  m_cellIDs = *cellsptr;
}


} // namespace k4::recCalo
