// This file's extension implies that it's C, but it's really -*- C++ -*-.
/**
 * @file RecCaloCommon/IDMapIndexer.h
 * @author scott snyder <snyder@bnl.gov>
 * @date Jan, 2026
 * @brief Implementation of ICaloIndexer using IDMap.
 */


#ifndef RECCALOCOMMON_IDMAPINDEXER_H
#define RECCALOCOMMON_IDMAPINDEXER_H


#include "RecCaloCommon/IDMap.h"
#include "RecCaloCommon/ICaloIndexer.h"
#include <stdexcept>


namespace k4::recCalo {


/**
 * @brief Implementation of ICaloIndexer using IDMap.
 *
 * This is an implementation of @c ICaloIndexer for a specific subdetector
 * based on @c IDMap with a fixed number of fields.
 */
template <unsigned NFIELDS>
class IDMapIndexer : public ICaloIndexer
{
public:
  /// Type of an index.
  using index_t = ICaloIndexer::index_t;

  /// Flag indicating an invalid index.
  static constexpr index_t INVALID = static_cast<index_t>(-1);

  /// Type of the mapping.
  using IDMap_t = IDMapN<index_t, NFIELDS>;

  /// Type describing a field; used to initialize the mapping.
  using FieldDesc_t = typename IDMap_t::FieldDesc_t;


  /**
   * @brief Constructor.
   * @param detID ID of the detector that we index.
   * @param detIDBits Number of bits in the cell IDs for the detector ID.
   * @param fields Set of fields to use from the identifiers.
   *               Must be sufficient, along with @c ignoredFields,
   *               to make identifiers unique.
   *               For best results, should be listed in order of increasing
   *               bit width.
   *               The size must be exactly @c NFIELDS.
   * @param ids Set of all identifiers to be indexed.
   *            Should be sorted for best results.
   * param sizeHint If non-zero, this is an estimate of the total size,
   *                in bytes, required by this mapping.  This will be
   *                used to reserve an appropriate size for the data vector.
   * @param ignoredFields Additional fields to ignore in order to make
   *               the identifiers unique.
   *
   * All entries in @c ids should be identical once the fields described
   * in @c fields and @c ignoredFields have been masked off; that is, any
   * additional fields must be identical for all ids.  When we try to find
   * an index, we first check that the extra bits match what we expect.
   */
  IDMapIndexer (int detID,
                size_t detIDBits,
                std::span<const FieldDesc_t> fields,
                std::span<const uint64_t> ids,
                size_t sizeHint = 0,
                std::span<const FieldDesc_t> ignoredFields = {});


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
  virtual size_t detIDBits() const override final;


  /**
   * @brief Return the allocated size of the mapping data, in bytes.
   */
  size_t byteSize() const;


  /**
   * @brief Format some statistics on the mapping.
   */
  void printStats(std::ostream& s) const;


private:
  /// The mapping.
  IDMap_t m_map;

  uint64_t m_otherFieldsMask;
  uint64_t m_otherFieldsVal;

  /// The ID of the detector we index.
  int m_detID;

  /// The number of bits in cell IDs for the detector ID.
  size_t m_detIDBits;

  /// The set of cells that we index.
  std::span<const uint64_t> m_cellIDs;
};


/**
 * @brief Constructor.
 */
template <unsigned NFIELDS>
IDMapIndexer<NFIELDS>::IDMapIndexer (int detID,
                                     size_t detIDBits,
                                     std::span<const FieldDesc_t> fields,
                                     std::span<const uint64_t> ids,
                                     size_t sizeHint /*= 0*/,
                                     std::span<const FieldDesc_t> ignoredFields /* = {}*/)
  : m_map (fields, INVALID, ids,
           [](size_t i) { return i; },
           sizeHint),
    m_detID (detID),
    m_detIDBits (detIDBits),
    m_cellIDs (ids)
{
  m_otherFieldsMask = ~ static_cast<uint64_t>(0);
  auto unmaskFields = [&] (std::span<const FieldDesc_t> ff) {
    for (const FieldDesc_t& f : ff) {
      unsigned offset = f.first;
      unsigned width = f.second;
      uint64_t fmask = (static_cast<uint64_t>(1) << width) - 1;
      m_otherFieldsMask &= ~ (fmask << offset);
    }
  };
  unmaskFields (fields);
  unmaskFields (ignoredFields);

  m_otherFieldsVal = 0;
  if (ids.size() > 0) {
    m_otherFieldsVal = ids[0] & m_otherFieldsMask;
    for (uint64_t id : ids) {
      if ((id & m_otherFieldsMask) != m_otherFieldsVal) {
        throw std::runtime_error ("IDMapIndexer: Inconsistent ID list");
      }
    }
  }
}


/**
 * @brief Return the index of an identifier.
 */
template <unsigned NFIELDS>
inline
auto IDMapIndexer<NFIELDS>::index (uint64_t id) const -> index_t
{
  if ((id & m_otherFieldsMask) != m_otherFieldsVal) return INVALID;
  return m_map.lookup (id);
}


/**
 * @brief Return the set of all identifiers that we index.
 */
template <unsigned NFIELDS>
inline
std::span<const uint64_t> IDMapIndexer<NFIELDS>::cellIDs() const
{
  return m_cellIDs;
}


/**
 * @brief Return the IDs of the detector(s) that we index.
 */
template <unsigned NFIELDS>
inline
std::span<const int> IDMapIndexer<NFIELDS>::detIDs() const
{
  return std::span<const int> (&m_detID, 1);
}


/**
 * @brief Number of bits in the cell ID used for the detector ID.
 */
template <unsigned NFIELDS>
inline
size_t IDMapIndexer<NFIELDS>::detIDBits() const
{
  return m_detIDBits;
}


/**
 * @brief Return the allocated size of the mapping data, in bytes.
 */
template <unsigned NFIELDS>
inline
size_t IDMapIndexer<NFIELDS>:: byteSize() const
{
  return m_map.byteSize();
}


/**
 * @brief Format some statistics on the mapping.
 */
template <unsigned NFIELDS>
inline
void IDMapIndexer<NFIELDS>::printStats(std::ostream& s) const
{
  return m_map.printStats(s);
}


} // namespace k4::recCalo


#endif // not RECCALOCOMMON_IDMAPINDEXER_H
