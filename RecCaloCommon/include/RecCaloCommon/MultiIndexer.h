// This file's extension implies that it's C, but it's really -*- C++ -*-.
/**
 * @file RecCaloCommon/MultiIndexer.h
 * @author scott snyder <snyder@bnl.gov>
 * @date Feb, 2026
 * @brief Implementation of ICaloIndexer for multiple detectors.
 */


#ifndef RECCALOCOMMON_MULTIINDEXER_H
#define RECCALOCOMMON_MULTIINDEXER_H


#include "RecCaloCommon/ICaloCellConstantsSvc.h"
#include "RecCaloCommon/ICaloIndexer.h"
#include <vector>
#include <stdexcept>


namespace k4::recCalo {


/**
 * @brief Implementation of ICaloIndexer for multiple detectors.
 *
 * This is an ICaloIndexer implementation that logically joins together
 * several detectors.  It is initialized with the ICaloIndexer objects
 * for each individual detector, along with the number of bits used
 * in the identifiers for the detector ID.  We reference the ICaloIndexer
 * instances given to us, but we need to make a vector of the concatenated lists
 * of IDs.  We store this in a given ICaloCellConstantsSvc also given
 * at initialization.
 */
class MultiIndexer : public ICaloIndexer
{
public:
  /// Type of an index.
  using index_t = ICaloIndexer::index_t;

  /// Flag indicating an invalid index.
  static constexpr index_t INVALID = static_cast<index_t>(-1);


  /**
   * @brief Constructor.
   * @param detIDBits Number of bits in the cell identifier used for the
   *                  detector ID.  (We assume that this starts at the
   *                  least significant bit.)
   *                  We require this to be no more than 8 to avoid
   *                  unreasonable memory allocation.s.
   * @param indexers  Collection of indexers for each detector.
   * @param constsSvc Used to store the concatenated vector of cell IDs.
   */
  MultiIndexer (size_t detIDBits,
                std::span<const ICaloIndexer* const> indexers,
                ICaloCellConstantsSvc& constsSvc);


  /**
   * @brief Return the index of an identifier.
   * @param id The identifier to look for.
   *
   * Returns the index of @c id in @c cellIDs(), or @c INVALID.
   */
  virtual index_t index (uint64_t id) const override final;


  /**
   * @brief Return the set of all identifiers that we index.
   */
  virtual std::span<const uint64_t> cellIDs() const override final;


  /**
   * @brief Return the IDs of the detector(s) that we index.
   */
  virtual std::span<const int> detIDs() const override final;


  /**
   * @brief Number of bits in the cell ID used for the detector ID.
   */
  virtual size_t detIDBits() const override;


  /**
   * @brief Exceptions thrown by the ctor.
   */
  class MultiIndexerException : public std::runtime_error
  {
  public:
    MultiIndexerException (const std::string& what);
  };


private:
  /// Number of detector ID bits.
  size_t m_detIDBits;

  /// Mask to extract the detector ID from a cell ID.
  uint64_t m_detIDMask;

  /// indexer/offset pairs. indexed by detector ID.
  /// To avoid having to do checks on the accesses, we size this vector
  /// to the maximum possible number of detector IDs, given the number
  /// of bits, and fill in empty entries with pointers to a dummy
  /// indexer which always returns INVALID.
  std::vector<std::pair<const ICaloIndexer*, size_t> > m_indexers;

  /// The IDs of the detectors that we index.
  std::vector<int> m_detIDs;

  /// List of cell IDs, concatenated over all detectors.
  /// The actual vector is held by the ICaloCellConstantsSvc.
  std::span<const uint64_t> m_cellIDs;
};


/**
 * @brief Return the index of an identifier.
 */
inline
auto MultiIndexer::index (uint64_t id) const -> index_t
{
  auto& p = m_indexers[id & m_detIDMask];
  index_t ndx = p.first->index (id);
  if (ndx != INVALID) [[likely]]
    return ndx + p.second;
  return INVALID;
}


/**
 * @brief Return the set of all identifiers that we index.
 */
inline
std::span<const uint64_t> MultiIndexer::cellIDs() const
{
  return m_cellIDs;
}


/**
 * @brief Return the IDs of the detector(s) that we index.
 */
inline
std::span<const int> MultiIndexer::detIDs() const
{
  return m_detIDs;
}


/**
 * @brief Number of bits in the cell ID used for the detector ID.
 */
inline
size_t MultiIndexer::detIDBits() const
{
  return m_detIDBits;
}


} // namespace k4::recCalo


#endif // not RECCALOCOMMON_MULTIINDEXER_H
