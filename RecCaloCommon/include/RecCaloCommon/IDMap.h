// This file's extension implies that it's C, but it's really -*- C++ -*-.
/**
 * @file RecCaloCommon/IDMap.h
 * @author scott snyder <snyder@bnl.gov>
 * @date Jan, 2026
 * @brief Fast lookup of identifiers.
 *
 * Each addressable component of the detector (calorimeter cells for example)
 * is assigned an `identifier', a collection of bit fields (representing
 * for example layer, eta, phi, etc) collected together into a single
 * 64-bit value.  For a given subdetector, the valid identifiers are
 * thus a sparse set of integers.  During reconstruction, we spend a lot
 * of time looking up identifiers in maps.  However, the highly-granular
 * FCC detectors have a very large number of cells (~ 2 million for the
 * Allegro Barrel Ecal for example), and using a std::map is very slow.
 * Even using a std::unordered_map tends to be slow --- although lookups
 * are O(1), it has poor locality of reference and thus poor cache behavior.
 * Further, usual implementations of both std::map and std::unordered_map
 * require a separate memory allocation for each contained item, resulting
 * in a large memory-management overhead.  IDMap takes advantage
 * of the field structure of identifiers to achieve very fast lookup,
 * especially in the important case where we access identifiers sequentially,
 * while at the same time trying to minimize memory requirements.
 *
 * IDMap acts as a map from identifiers to some arbitrary POD type.
 * In practice, however, if one has several ID->T mappings, it may be better
 * to have a single IDMap mapping the identifiers to a dense set of integers,
 * and then keep the T's in std::vector's.
 *
 * The payload objects must be `POD' types; that is, no constructor
 * or destructors, and the objects may be copied via a bitwise copy.
 *
 * An IDMap is constructed from the following:
 *  - A list of the fields from the identifier to be used for indexing,
 *    each specified as an offset and a width.  The identifiers for
 *    this subdetector should be unique within this set of fields.
 *    For best performance, fields should be listed in order from smaller
 *    to larger widths.
 *  - One specific value of the payload type that is invalid.  (This will
 *    be used to indicate non-existent entries.)
 *  - The list of all identifiers for this subdetector.  Should be sorted
 *    for best results.
 *  - A function to return the payload value for a given index in the
 *    identifier list.  This will be used to initialize the mapping.
 *
 * One can then call lookup() with an identifier key.  This will return
 * either the mapped value or the invalid value.
 *
 * In principle, it would be possible to allow inserting/changing mappings
 * after the map has been created, but I currently don't have a use-case
 * for that, so let's just keep things simple for now.
 *
 * IDMap<PAYLOAD> will work for a variable number of fields (up to a
 * compile-time maximum).  However, for a small number of fields, performance
 * can be significantly better if the number of fields is fixed
 * at compile time.  That can be implemented with IDMapN<PAYLOAD, N>
 * (deriving from IDMap<PAYLOAD>).  Here, N is the number of fields,
 * and the constructor will abort if an incorrect number of fields
 * is provided.  Otherwise it works the same.
 *
 *
 * Implementation notes:
 *
 * We use a trie data structure, with each level of node corresponding
 * to one field.  For example, with three fields we have:
 *
 * Field 0
 *                         +-------------+
 *                         |  Root node  |
 *                         |0123456789...|
 *                         +--|----|-----+
 *                            |    |
 * Field 1   +----------------+    +---------------+
 *           v                                     v
 *     +-------------+                       +-------------+
 *     |   Node 2    |          ...          |   Node 7    |
 *     |0123456789...|                       |0123456789...|
 *     +-------------+                       +----|--|-----+
 *                            +-------------------+  +------------+
 * Field 2                    v                                   v
 *                      +-----------+                            ...
 *                      |   Node 7,4|
 *                      | 01234     |
 *                      | vvxxv ... |
 *                      | aaxxa     |
 *                      | llxxl     |
 *                      +-----------+
 *
 * The root node has size given by the number of possible values for field 0,
 * and contains pointers to level-1 nodes.  Similarly, level-1 nodes have
 * sizes given by the number of possible values for field 1 and contain
 * pointers to level-2 nodes.  The level-2 (left) nodes contain instances
 * of the payload type, set to the invalid value for nonexistent entries.
 * The sizes of these nodes are found by scanning the provided list of
 * identifiers for the maximum used value of field 2.
 *
 * For compactness and locality, the nodes are not allocated individually,
 * but rather stored in a single std::vector<char>.  The pointers between
 * nodes are stored as 32-bit indices into this vector (with 0 then being
 * usable as a null value, since it would otherwise reference the root node).
 *
 *
 * Performance measurements:
 *
 * Some tests were done to quantify the performance of the lookup, implemented
 * as part of the IDMap_test unit test when run with the --perf switch.
 * The test makes a map of identifiers corresponding to the Allegro
 * ECal barrel as of this writing (about 2 million identifiers) and then
 * does lookups the identifiers, both sequentially and in random order,
 * each repeated the default 100 times.
 *
 * A number of different types of map were compared:
 *
 *  - null: Loops through and selects identifiers, but doesn't actually
 *          do any lookup.  To get an idea of the overhead due to the
 *          test itself.
 *  - IDMap: IDMap<> with fields layer, theta, module.
 *  - IDMapN: The same, but with the number of fields fixed to 3
 *            at compile-time.
 *  - map: std::map
 *  - unordered_map: std::unordered_map
 *  - array: A simple 3-dimensional array.
 *           Let N1, N2, and N3 be the maximum values observed in the set
 *           of identifiers for the three fields.  Then given an identifier
 *           with field values x1, x2, and x3 we construct an index as
 *                     N3 * (x1 * N2 + x2) + x3
 *           into a flat array.
 *
 * Tests were run on an intel i7-1360P cpu, compiled with a pre-release
 * version of gcc16 (20260103), testing combinations of
 * -march=x86_64-v2 / -v3 and -O2 / -O3.
 *
 * Results:
 *
 *                        -v2,-O2   -v3,-O2   -v2,-O3   -v3,-O3
 *    null          seq    0.08      0.08      0.06      0.05
 *                  rand   0.71      0.65      0.69      0.59
 *
 *    IDMap         seq    0.50      0.46      0.39      0.39
 *                  rand   6.65      6.76      5.83      5.64
 *
 *    IDMapN        seq    0.31      0.26      0.31      0.25
 *                  rand   4.73      4.44      4.73      4.49
 *
 *    map           seq   25.77     27.60     27.51     27.55
 *                  rand 205.06    207.93    207.82    197.45
 *
 *    unordered_map seq    3.34     3.20       3.20      3.30
 *                  rand   9.89     9.78       9.97      9.94
 *
 *    array         seq    0.31     0.23       0.31      0.23
 *                  rand   4.46     4.07       4.49      4.04
 *
 * Some observations:
 *
 *  - std::map is really too slow to even be considered.
 *  - std::unordered_map is about twice as slow is IDMap for random access,
 *    but about 10 times slower for sequential access.  This is consistent
 *    with the hash table implementation of std::unordered_map not giving
 *    good locality of reference.
 *  - Fixing the number of fields at compile time gives about a 20-40%
 *    improvement.
 *  - Playing with the optimization/code generation options can give
 *    a ~ 5-20% improvement.
 *  - The simple array case is somewhat faster then IDMap and slightly
 *    faster than IDMapN.  However, it takes more than three times
 *    as much memory (~15M for IDMap vs. ~50M for array).
 *
 * Overall, for most cases, IDMapN seems like a good option.
 */


#ifndef RECCALOCOMMON_IDMAP_H
#define RECCALOCOMMON_IDMAP_H


#include <type_traits>
#include <concepts>
#include <vector>
#include <span>
#include <utility>
#include <ostream>
#include <cstdlib>
#include <cstdint>
#include <iostream>


namespace k4::recCalo {


// Check if a type is usable as a payload.
// Basically, is it a POD type.  But both std::is_pod and is_trivial
// are deprecated, so do it this way.
template <class PAYLOAD>
concept IDMapPayload =
  std::is_trivially_default_constructible_v<PAYLOAD> &&
  std::is_trivially_copyable_v<PAYLOAD> &&
  std::is_trivially_assignable_v<PAYLOAD&, PAYLOAD>;


/**
 * @brief Fast lookup of identifiers.
 *        PAYLOAD is the mapped (value) type.
 */
template <IDMapPayload PAYLOAD>
class IDMap
{
public:
  /// Maximum supported number of fields.
  constexpr static size_t MAXFIELDS = 6;

  /// Type of keys (identifiers).
  using key_t = uint64_t;

  /// Type of payload (mapped value).
  using payload_t = PAYLOAD;

  /// Describe a single field as bit offset and bit width.
  using FieldDesc_t = std::pair<size_t, size_t>;


  /**
   * @brief Helper for making a @c FieldDesc_t from a @c BitFieldElement.
   *
   * Makes a @c FieldDesc_t from something with @c offset() and @c width()
   * methods.
   */
  template <class T>
  static FieldDesc_t makeDesc (const T& f);


  /**
   * @brief Constructor.
   * @param fields Set of fields to use from the identifiers.
   *               Must be sufficient to make identifiers unique.
   *               For best results, should be listed in order of increasing
   *               bit width.
   *               The number must not be more than MAXFIELDS; if nfields
   *               is specified, then the number must match that.
   * @param invalid An instance of the payload type to use to indicate
   *                an invalid value.
   * @param ids Set of all identifiers to be mapped.
   *            Should be sorted for best results.
   * @param valfunc The value associated with ids[i] is found by
   *                calling valfunc(i).
   * param sizeHint If non-zero, this is an estimate of the total size,
   *                in bytes, required by this mapping.  This will be
   *                used to reserve an appropriate size for the data vector.
   * @param nfields If nonzero, give an error if the size of the fields
   *                argument does not match nfields.  Intended for use
   *                by the IDMapN derived classes.
   */
  template <std::invocable<size_t> VALFUNC>
  requires std::convertible_to<std::result_of_t<VALFUNC(size_t)>, PAYLOAD>
  IDMap (std::span<const FieldDesc_t> fields,
         payload_t invalid,
         std::span<const key_t> ids,
         VALFUNC valfunc,
         size_t sizeHint = 0,
         size_t nfields = 0);


  /**
   * @brief Look up a value in the mapping.
   * @param k The value to look up.
   * @returns The mapped value, or the invalid value if not found.
   */
  payload_t lookup (key_t k) const;


  /**
   * @brief Return the number of keys in the mapping.
   */
  size_t size() const;


  /**
   * @brief Return the allocated size of the mapping data, in bytes.
   */
  size_t byteSize() const;


  /**
   * @brief Format some statistics on the mapping.
   */
  void printStats(std::ostream& s) const;


protected:
  /**
   * @brief Look up a value in the mapping.
   * @param k The value to look up.
   * @param nfields Number of fields in the mapping.
   * @returns The mapped value, or the invalid value if not found.
   *
   * This is to allow code to be generated with the number of fields
   * fixed at compile-time.
   */
  payload_t lookup (key_t k, size_t nfields) const;


private:
  /// Type used for a node index.
  using index_t = uint32_t;


  /**
   * @brief Insert a new value in the map.
   * @param k The key (identifier) of the new mapping.
   * @param v The value of the new mapping.
   */
  void insert (key_t k, payload_t v);


  /**
   * @brief Describe one field.
   *        For each field, we store a bit mask, a shift count,
   *        and the size of the corresponding trie node.
   */
  struct Field
  {
    Field (unsigned pos=0, unsigned nbits=0)
      : m_mask ( ((1ull<<nbits)-1) << pos),
        m_shift (pos)
    {
    }
    unsigned size() const { return m_size; }
    unsigned extract(key_t x) const { return ((x & m_mask) >> m_shift); }
    void updateSize(key_t x) { m_size = std::max (m_size, extract(x)+1); }

    key_t m_mask;
    unsigned m_shift;
    unsigned m_size = 0;
  };


  /// Number of fields.
  size_t m_nfields;

  /// All fields.
  Field m_fields[6];

  /// The mapping data.
  std::vector<char> m_data;

  /// The invalid value.
  payload_t m_invalid;

  /// Some basic statistics, giving how many blocks and filled entries
  /// are present at each level of the trie.
  struct Stat
  {
    unsigned nblock = 0;
    unsigned filled = 0;
  };
  Stat m_stats[MAXFIELDS];
};


/**
 * @brief Version of IDMap specialized for a fixed number of fields.
 *        PAYLOAD is the mapped (value) type.
 *        NFIELDS is the number of fields.
 *
 *        Fixing the number of fields at compile time results
 *        in faster code.
 */
template <IDMapPayload PAYLOAD, size_t NFIELDS>
class IDMapN
  : public IDMap<PAYLOAD>
{
public:
  /// Bring in types from the base IDMap.
  using base_class = IDMap<PAYLOAD>;
  using typename base_class::key_t;
  using typename base_class::payload_t;
  using typename base_class::FieldDesc_t;


  /**
   * @brief Constructor.
   * @param fields Set of fields to use from the identifiers.
   *               Must be sufficient to make identifiers unique.
   *               For best results, should be listed in order of increasing
   *               bit width.
   *               The number must not be more than MAXFIELDS; if nfields
   *               is specified, then the number must match that.
   * @param invalid An instance of the payload type to use to indicate
   *                an invalid value.
   * @param ids Set of all identifiers to be mapped.
   *            Should be sorted for best results.
   * @param valfunc The value associated with ids[i] is found by
   *                calling valfunc(i).
   * param sizeHint If non-zero, this is an estimate of the total size,
   *                in bytes, required by this mapping.  This will be
   *                used to reserve an appropriate size for the data vector.
   */
  template <std::invocable<size_t> VALFUNC>
  requires std::convertible_to<std::result_of_t<VALFUNC(size_t)>, PAYLOAD>
  IDMapN (std::span<const FieldDesc_t> fields,
          payload_t invalid,
          std::span<const key_t> ids,
          VALFUNC valfunc,
          size_t sizeHint = 0)
    : base_class (fields, invalid, ids, valfunc, sizeHint, NFIELDS)
  {
    if (fields.size() != NFIELDS) std::abort();
  }


  /**
   * @brief Look up a value in the mapping.
   * @param k The value to look up.
   * @returns The mapped value, or the invalid value if not found.
   */
  payload_t lookup (key_t k) const { return base_class::lookup (k, NFIELDS); }
};


/**
 * @brief Helper for making a @c FieldDesc_t from a @c BitFieldElement.
 */
template <IDMapPayload PAYLOAD>
template <class T>
inline
auto IDMap<PAYLOAD>::makeDesc (const T& f) -> FieldDesc_t
{
  return std::make_pair (f.offset(), f.width());
}


/**
 * @brief Constructor.
 */
template <IDMapPayload PAYLOAD>
template <std::invocable<size_t> VALFUNC>
requires std::convertible_to<std::result_of_t<VALFUNC(size_t)>, PAYLOAD>
IDMap<PAYLOAD>::IDMap (std::span<const FieldDesc_t> fields,
                       payload_t invalid,
                       std::span<const key_t> ids,
                       VALFUNC valfunc,
                       size_t sizeHint /*= 0*/,
                       size_t nfields /*= 0*/)
  : m_invalid (invalid)
{
  if (fields.size() > MAXFIELDS || fields.size() == 0 ||
      (nfields != 0 && fields.size() != nfields)) {
    std::cerr << "IDMap: Bad field list length; got " << fields.size()
              << "; expected ";
    if (nfields && nfields <= MAXFIELDS && nfields > 0)
      std::cerr << nfields;
    else
      std::cerr << "< " << MAXFIELDS << " and > 0";
    std::cerr << "\n";
    std::abort();
  }

  // Set the field offsets and widths.
  m_nfields = fields.size();
  for (size_t i = 0; i < fields.size(); ++i) {
    if (fields[i].second > 20 || fields[i].first + fields[i].second >= 64) {
      std::cerr << "IDMap: Bad field offset/width: " << fields[i].first
                << " " << fields[i].second << "\n";
      std::abort();
    }
    m_fields[i] = Field(fields[i].first, fields[i].second);
  }

  // For all fields but the last, set the size based on the bit width.
  for (size_t i = 0; i < m_nfields-1; ++i)
    m_fields[i].m_size = 1u << fields[i].second;

  // Set the size of the last field based on the maximum observed value
  // in the identifier list.
  Field& lastfield = m_fields[m_nfields-1];
  for (key_t id : ids) {
    lastfield.updateSize (id);
  }

  // Allocate the root block.
  if (sizeHint > 0) {
    m_data.reserve (sizeHint);
  }
  m_data.resize (m_fields[0].size() * sizeof(index_t));
  m_stats[0].nblock = 1;

  // Initialize the mapping.
  for (size_t i = 0; i < ids.size(); ++i) {
    insert (ids[i], valfunc(i));
  }
}


/**
 * @brief Look up a value in the mapping.
 */
template <IDMapPayload PAYLOAD>
inline
auto IDMap<PAYLOAD>::lookup (key_t k) const -> payload_t
{
  return lookup (k, m_nfields);
}


/**
 * @brief Return the number of keys in the mapping.
 */
template <IDMapPayload PAYLOAD>
inline
size_t IDMap<PAYLOAD>::size() const
{
  // Size is just the number of insertions we've made in leaf nodes.
  return m_stats[m_nfields-1].filled;
}


/**
 * @brief Return the allocated size of the mapping data, in bytes.
 */
template <IDMapPayload PAYLOAD>
inline
size_t IDMap<PAYLOAD>::byteSize() const
{
  return m_data.size();
}


/**
 * @brief Format some statistics on the mapping.
 */
template <IDMapPayload PAYLOAD>
void IDMap<PAYLOAD>::printStats (std::ostream& s) const
{
  unsigned tot = 0;
  for (size_t i = 0; i < m_stats.size(); ++i) {
    s << "trie level " << i << " "
      << m_stats[i].nblock << " blocks of size "
      << m_stats[i].size << " for total "
      << m_stats[i].nblock*m_stats[i].size << " with "
      << m_stats[i].filled << " ("
      << (m_stats[i].filled ? static_cast<float>(m_stats[i].filled)/(m_stats[i].nblock*m_stats[i].size)*100 : 0)
      << "%) filled" << "\n";
     size_t this_sz = i == m_stats.size()-1 ? sizeof(payload_t) : sizeof(void*);
     tot += m_stats[i].nblock*m_stats[i].size * this_sz;
  }
  s << "tot " << (float)tot/1024/1024 << "MB\n";
}


/**
 * @brief Look up a value in the mapping.
 */
template <IDMapPayload PAYLOAD>
auto IDMap<PAYLOAD>::lookup (key_t k, size_t nfields) const -> payload_t
{
  // Start of the mapping data.
  const char* data = m_data.data();

  // Index of the node we're currently looking at.
  index_t pos = 0;

  // Start with the first field.
  const Field* f = m_fields;

  // Extract it from the identifier.
  unsigned ndx = f->extract (k);

  // Looping over fields.
  for (size_t ifield = 1; ifield < nfields; ++ifield) {
    ++f;

    // Pointer to the node we're currently looking at.
    const index_t* node = reinterpret_cast<const index_t*> (data + pos);

    // Index of the node in the next trie level; return if it's null.
    index_t next_pos = node[ndx];
    if (!next_pos) [[unlikely]] {
      return m_invalid;
    }

    // Looking at the next node, and extract the value of the next field.
    pos = next_pos;
    ndx = f->extract (k);
  }

  // Check the bounds for the leaf node.
  if (ndx >= f->size()) [[unlikely]] return m_invalid;

  // Return the value from the leaf node.
  const payload_t* leaf = reinterpret_cast<const payload_t*> (data + pos);
  return leaf[ndx];
}


/**
 * @brief Insert a new value in the map.
 */
template <IDMapPayload PAYLOAD>
void IDMap<PAYLOAD>::insert (key_t k, payload_t v)
{
  if (v == m_invalid) {
    std::cerr << "IDMap: Attempt to insert invalid value.\n";
    std::abort();
  }

  // Index of the node we're currently looking at.
  index_t pos = 0;

  // Extract the value of the first field.
  unsigned ndx = m_fields[0].extract (k);

  for (size_t ifield = 1; ifield < m_nfields; ++ifield) {
    const Field& f = m_fields[ifield];

    // Pointer to the (non-leaf) node we're looking at.
    index_t* node = reinterpret_cast<index_t*> (m_data.data() + pos);

    // Get the pointer from that node.
    index_t next_pos = node[ndx];
    if (!next_pos) {
      // It hasn't been set yet.  We need to make a new node.
      // It goes at the end of the data vector.
      next_pos = m_data.size();
      if (ifield < m_nfields-1) {
        // A non-leaf node.
        m_data.resize (m_data.size() + f.size() * sizeof(index_t));
      }
      else {
        // A leaf node.
        m_data.resize (m_data.size() + f.size() * sizeof(payload_t));
        payload_t* data = reinterpret_cast<payload_t*> (m_data.data() + next_pos);
        std::fill (data, data+f.size(), m_invalid);
      }

      // Get the pointer to our current node again; it may have moved.
      node = reinterpret_cast<index_t*> (m_data.data() + pos);

      // Fill in the node pointer.
      node[ndx] = next_pos;

      // Maintain some statistics.
      m_stats[ifield].nblock += 1;
      m_stats[ifield-1].filled += 1;
    }

    // Move to the next field.
    pos = next_pos;
    ndx = f.extract (k);
  }

  // Pointer to the leaf field.
  payload_t* leaf = reinterpret_cast<payload_t*> (m_data.data() + pos);

  // Fill in the value.
  if (leaf[ndx] == m_invalid) {
    m_stats[m_nfields-1].filled += 1;
  }
  leaf[ndx] = v;
}


} // namespace k4::recCalo


#endif // not RECCALOCOMMON_IDMAP_H
