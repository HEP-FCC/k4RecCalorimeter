/**
 * @file RecCaloCommon/tests/IDMap_test.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date Jan, 2026
 * @brief Unit test for IDMap.
 */

#undef NDEBUG
#include "RecCaloCommon/IDMap.h"
#include "DD4hep/IDDescriptor.h"
#include "boost/timer/timer.hpp"
#include <vector>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <cstdint>
#include <cassert>


// Key type.
using mapkey_t = uint64_t; // libc defines key_t...
using mapkey_span = std::span<const mapkey_t>;

using payload_t = uint32_t;


// Helpers to make lists of IDS.
#include "make_ecalb_ids.icc"
#include "make_ecale_ids.icc"


//************************************************************************
// Very simple RNG that should be repeatable across architectures.
//

static const uint32_t rngmax = static_cast<uint32_t> (-1);


inline
uint32_t rng_seed (uint32_t& seed)
{
  seed = (1664525*seed + 1013904223);
  return seed;
}


inline
float randf_seed (uint32_t& seed, float rmax, float rmin = 0)
{
  return static_cast<float>(rng_seed(seed)) / static_cast<float>(rngmax) * (rmax-rmin) + rmin;
}


inline
int randi_seed (uint32_t& seed, int rmax, int rmin = 0)
{
  return static_cast<int> (randf_seed (seed, rmax, rmin));
}


//*******************************************************************
// Basic test of map lookup.
//


void test1 (mapkey_span ids)
{
  using Map_t = k4::recCalo::IDMap<payload_t>;
  using FieldDesc_t = Map_t::FieldDesc_t;
  dd4hep::IDDescriptor desc ("desc", ecalb_descstr);

  // Make the fields used by the map.
  std::vector<FieldDesc_t> fielddescs;
  fielddescs.push_back (Map_t::makeDesc (*desc.field ("layer")));
  fielddescs.push_back (Map_t::makeDesc (*desc.field ("theta")));
  fielddescs.push_back (Map_t::makeDesc (*desc.field ("module")));

  // Make and fill the map.
  const unsigned INVALID = static_cast<unsigned>(-1);
  Map_t map (fielddescs, INVALID, ids, [](size_t i) { return i+3; });
  assert (map.size() == ids.size());

  // Test lookup of each input identifier.
  for (size_t i = 0; i < ids.size(); ++i) {
    assert (map.lookup (ids[i]) == i+3);
  }

  // Test lookup failure for keys not in the input list.
  uint32_t seed = 1234;
  dd4hep::BitFieldCoder* decoder = desc.decoder();
  size_t layer_index = decoder->index ("layer");
  size_t module_index = decoder->index ("module");
  size_t theta_index  = decoder->index ("theta");
  unsigned nmodule = 1 << desc.field ("module")->width();
  unsigned ntheta = 1 << desc.field ("theta")->width();
  unsigned ntry = 0;
  for (size_t i = 0; i < 1000; i++) {
    mapkey_t id = ids[0];
    decoder->set (id, layer_index, randi_seed (seed, ecalb_numLayers + 2));
    decoder->set (id, module_index, randi_seed (seed, nmodule));
    decoder->set (id, theta_index, randi_seed (seed, ntheta));
    if (!std::ranges::binary_search (ids, id)) {
      assert (map.lookup (id) == INVALID);
      ++ntry;
    }
  }
  assert (ntry > 0);
}


//*******************************************************************
// Run performance tests comparing IDMap with other map types.


// -- IDMap test jigs


using FieldDesc_t = k4::recCalo::IDMap<payload_t>::FieldDesc_t;

template<class IDMAP>
class IDMapLookupB
  : public IDMAP
{
public:
  constexpr static payload_t INVALID = static_cast<payload_t> (-1);
  using FieldDesc_t = typename IDMAP::FieldDesc_t;
  static std::vector<FieldDesc_t> fieldDescs();
  IDMapLookupB (mapkey_span ids)
    : IDMAP (fieldDescs(), INVALID, ids, [](size_t i) { return i; })
  { }
};


template<class IDMAP>
auto IDMapLookupB<IDMAP>::fieldDescs() -> std::vector<FieldDesc_t>
{
  std::vector<FieldDesc_t> fielddescs;
  dd4hep::IDDescriptor desc ("desc", ecalb_descstr);
  auto pushdesc = [&] (const std::string s) {
    const dd4hep::BitFieldElement* bfe = desc.field (s);
    fielddescs.emplace_back (bfe->offset(), bfe->width());
  };
  pushdesc ("layer");
  pushdesc ("theta");
  pushdesc ("module");
  return fielddescs;
}


template<class IDMAP>
class IDMapLookupE
  : public IDMAP
{
public:
  constexpr static payload_t INVALID = static_cast<payload_t> (-1);
  using FieldDesc_t = typename IDMAP::FieldDesc_t;
  static std::vector<FieldDesc_t> fieldDescs();
  IDMapLookupE (mapkey_span ids)
    : IDMAP (fieldDescs(), INVALID, ids, [](size_t i) { return i; })
  { }
};


template<class IDMAP>
auto IDMapLookupE<IDMAP>::fieldDescs() -> std::vector<FieldDesc_t>
{
  std::vector<FieldDesc_t> fielddescs;
  dd4hep::IDDescriptor desc ("desc", ecale_descstr);
  auto pushdesc = [&] (const std::string s) {
    const dd4hep::BitFieldElement* bfe = desc.field (s);
    fielddescs.emplace_back (bfe->offset(), bfe->width());
  };


  // Combine side+wheel into a single field.
  const dd4hep::BitFieldElement* bfe_side = desc.field ("side");
  const dd4hep::BitFieldElement* bfe_wheel = desc.field ("wheel");
  fielddescs.emplace_back (std::min(bfe_side->offset(), bfe_wheel->offset()),
                           bfe_side->width() + bfe_wheel->width());


  pushdesc ("z");
  pushdesc ("rho");
  pushdesc ("module");
  return fielddescs;
}


// -- std::map test jig


class MapLookup
{
public:
  static constexpr payload_t INVALID = static_cast<payload_t>(-1);

  MapLookup (mapkey_span ids)
  {
    for (mapkey_t i = 0; i < ids.size(); ++i)
      m_map[ids[i]] = i;
  }

  payload_t lookup (mapkey_t id) const
  {
    auto it = m_map.find (id);
    if (it != m_map.end()) return it->second;
    return INVALID;
  }

  size_t byteSize() const { return 0; }

private:
  std::map<mapkey_t,  payload_t> m_map;
};


// -- std::unordered_map test jig


class UOMapLookup
{
public:
  static constexpr payload_t INVALID = static_cast<payload_t>(-1);

  UOMapLookup (mapkey_span ids)
  {
    for (mapkey_t i = 0; i < ids.size(); ++i)
      m_map[ids[i]] = i;
  }

  payload_t lookup (mapkey_t id) const
  {
    auto it = m_map.find (id);
    if (it != m_map.end()) return it->second;
    return INVALID;
  }

  size_t byteSize() const { return 0; }

private:
  std::unordered_map<mapkey_t,  payload_t> m_map;
};


// -- Array lookup test jig


class ArrLookup
{
public:
  static constexpr payload_t INVALID = static_cast<payload_t>(-1);
  static constexpr unsigned NFIELDS = 3;

  ArrLookup (mapkey_span ids);
  uint32_t index (mapkey_t id) const;
  payload_t lookup (mapkey_t id) const;
  size_t byteSize() const { return m_arr.size() * sizeof(payload_t); }


private:
  struct Field
  {
    Field (unsigned pos=0, unsigned nbits=1)
      : m_mask ( ((1ull<<nbits)-1) << pos),
        m_shift (pos),
        m_size (0)
    {
    }
    size_t size() const { return m_size; }
    size_t extract(size_t x) const { return ((x & m_mask) >> m_shift); }
    void updateSize(size_t x) { m_size = std::max (m_size, (unsigned)extract(x)+1); }
    size_t m_mask;
    unsigned m_shift;
    unsigned m_size;
  };

  std::vector<payload_t> m_arr;
  Field m_fields[NFIELDS];
  uint32_t m_tot_sz;
};


inline
uint32_t ArrLookup::index (mapkey_t k) const
{
  const Field* f = m_fields;
  uint32_t ndx = f->extract(k);
  if (ndx >= f->size()) [[unlikely]] return m_tot_sz;
  for (size_t ifield = 1; ifield < NFIELDS; ++ifield)
  {
    ++f;
    uint32_t this_ndx = f->extract(k);
    if (this_ndx >= f->size()) [[unlikely]] return m_tot_sz;
    ndx = ndx * f->size() + this_ndx;
  }
  return ndx;
}


ArrLookup::ArrLookup (mapkey_span ids)
{
  m_fields[0] = Field (11, 8); // layer
  m_fields[1] = Field (30, 10); // theta
  m_fields[2] = Field (19, 11); // module

  m_fields[NFIELDS-1].m_size = 0;
  for (mapkey_t id : ids) {
    for (Field& f : m_fields)
      f.updateSize (id);
  }

  size_t tot_sz = 1;
  for (Field& f : m_fields)
    tot_sz *= f.size();
  m_tot_sz = tot_sz;

  m_arr.resize (tot_sz+1, INVALID);

  for (size_t i = 0; i < ids.size(); ++i)
    m_arr[index(ids[i])] = i;
}


inline
payload_t ArrLookup::lookup (mapkey_t k) const
{
  return m_arr[index(k)];
}


// -- skeleton


class Timer
{
public:
  Timer();

  class RunTimer
  {
  public:
    RunTimer (boost::timer::cpu_timer& timer) : m_timer (&timer)
    { timer.resume(); }
    RunTimer (RunTimer&& other) : m_timer (other.m_timer) { other.m_timer = nullptr; }
    ~RunTimer() { if (m_timer) m_timer->stop(); }
  private:
    boost::timer::cpu_timer* m_timer;
  };
  RunTimer run() { return RunTimer (m_timer); }

  std::string format() const { return m_timer.format(3); }


private:
  boost::timer::cpu_timer m_timer;
};


Timer::Timer()
{
  m_timer.stop();
}


class TesterBase
{
public:
  TesterBase (const std::string& name, mapkey_span ids);

  // Loop choosing IDs but not doing any lookup, in order to measure
  // the overhead of the tests.
  size_t test (size_t n);
  size_t test_null_seq (size_t n);
  size_t test_null_rand (size_t n);

  void report() const;


protected:
  size_t m_byteSize = 0;
  std::string m_name;
  mapkey_span m_ids;
  Timer m_seq_timer;
  Timer m_rand_timer;
};


TesterBase::TesterBase (const std::string& name, mapkey_span ids)
  : m_name (name),
    m_ids (ids)
{
}


void TesterBase::report() const
{
  std::cout << m_name;
  if (m_byteSize) {
    std::cout << " " << m_byteSize/1024./1024 << "MB";
  }
  std::cout << "\n";
  std::cout << "seq:  " << m_seq_timer.format();
  std::cout << "rand: " << m_rand_timer.format();
}


size_t TesterBase::test (size_t n)
{
  size_t out = 0;
  out += test_null_seq (n);
  out += test_null_rand (n);
  return out;
}


size_t TesterBase::test_null_seq (size_t n)
{
  size_t out = 0;
  size_t sz = m_ids.size();
  auto timer = m_seq_timer.run();
  for (size_t j = 0; j < n; ++j) {
    for (size_t i = 0; i < sz; ++i) {
      out += m_ids[i];
    }
  }
  return out;
}


size_t TesterBase::test_null_rand (size_t n)
{
  size_t out = 0;
  size_t sz = m_ids.size();
  uint32_t seed = 1234;
  auto timer = m_rand_timer.run();
  for (size_t i = 0; i < n*sz; i++) {
    size_t pos = randi_seed (seed, sz-1);
    out += m_ids[pos];
  }
  return out;
}


template <class LOOKUP>
class Tester
  : public TesterBase
{
public:
  Tester (const std::string& name, mapkey_span ids);
  size_t test (size_t n);
  void test_seq (size_t n);
  void test_rand (size_t n);


private:
  LOOKUP m_lookup;
};


template <class LOOKUP>
Tester<LOOKUP>::Tester (const std::string& name, mapkey_span ids)
  : TesterBase (name, ids),
    m_lookup (ids)
{
  m_byteSize = m_lookup.byteSize();
}


template <class LOOKUP>
size_t Tester<LOOKUP>::test (size_t n)
{
  test_seq (n);
  test_rand (n);
  return 0;
}


template <class LOOKUP>
void Tester<LOOKUP>::test_seq (size_t n)
{
  size_t sz = m_ids.size();
  auto timer = m_seq_timer.run();
  for (size_t j = 0; j < n; ++j) {
    for (size_t i = 0; i < sz; ++i) {
      if (m_lookup.lookup (m_ids[i]) != i) [[unlikely]] {
        std::cout << m_name << " seq " << i << " " << m_ids[i] << " " << m_lookup.lookup (m_ids[i]) << "\n";
        std::abort();
      }
    }
  }
}


template <class LOOKUP>
void Tester<LOOKUP>::test_rand (size_t n)
{
  size_t sz = m_ids.size();
  uint32_t seed = 1234;
  auto timer = m_rand_timer.run();
  for (size_t i = 0; i < n*sz; i++) {
    size_t pos = randi_seed (seed, sz-1);
    if (m_lookup.lookup (m_ids[pos]) != pos) [[unlikely]] {
      std::cout << m_name << " rand " << pos << " " << m_ids[pos] << " " << m_lookup.lookup (m_ids[pos]) << "\n";
      std::abort();
    }
  }
}


template <class TESTER>
size_t dotest (const char* name, mapkey_span ids, size_t n)
{
  TESTER t (name, ids);
  size_t ret = t.test (n);
  t.report();
  return ret;
}


size_t perftest (bool barrel, mapkey_span ids, size_t n)
{
  using IDMap_t = k4::recCalo::IDMap<payload_t>;
  using IDMap3_t = k4::recCalo::IDMapN<payload_t, 3>;
  using IDMap4_t = k4::recCalo::IDMapN<payload_t, 4>;

  size_t ret = 0;  // To prevent tests from being optimized away...
  ret += dotest<TesterBase> ("null", ids, n);
  if (barrel) {
    ret += dotest<Tester<IDMapLookupB<IDMap_t> > > ("IDMapB", ids, n);
    ret += dotest<Tester<IDMapLookupB<IDMap3_t> > > ("IDMapB3", ids, n);
  }
  else {
    ret += dotest<Tester<IDMapLookupE<IDMap_t> > > ("IDMapE", ids, n);
    ret += dotest<Tester<IDMapLookupE<IDMap4_t> > > ("IDMapE4", ids, n);
  }
  ret += dotest<Tester<MapLookup> > ("std::map", ids, n);
  ret += dotest<Tester<UOMapLookup> > ("std::unordered_map", ids, n);
  if (barrel)
    ret += dotest<Tester<ArrLookup> > ("simple array", ids, n);
  return ret;
}


//*******************************************************************


int main (int argc, char** argv)
{
  std::vector<mapkey_t> ecale_ids = make_ecale_ids();

  if (argc >= 2 && std::string(argv[1]).starts_with ("--perf")) {
    size_t n = 0;
    if (argc >= 3) n = atoi (argv[2]);
    if (n == 0) n = 100;
    if (std::string(argv[1]).starts_with ("--perf-endcap")) {
      std::vector<mapkey_t> ids = make_ecale_ids();
      perftest (false, ids, n);
    }
    else {
      std::vector<mapkey_t> ids = make_ecalb_ids();
      perftest (true, ids, n);
    }
  }
  else {
    std::vector<mapkey_t> ids = make_ecalb_ids();
    test1 (ids);
  }

  return 0;
}
