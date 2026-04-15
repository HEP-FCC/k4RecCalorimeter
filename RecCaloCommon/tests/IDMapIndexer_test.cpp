/**
 * @file RecCaloCommon/tests/IDMapIndexer_test.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date Jan, 2026
 * @brief Unit test for IDMapIndexer.
 */

#undef NDEBUG
#include "RecCaloCommon/IDMapIndexer.h"
#include "DD4hep/IDDescriptor.h"
#include <vector>
#include <algorithm>
#include <span>
#include <cassert>


using mapkey_t = uint64_t; // libc defines key_t...
using mapkey_span = std::span<const mapkey_t>;


// Helper to make lists of IDS.
#include "make_ecalb_ids.icc"
#include "make_ecale_ids.icc"


#define EXPECT_EXCEPTION(EXC, CODE) do { \
  bool caught = false;                   \
  try {                                  \
    CODE;                                \
  }                                      \
  catch (const EXC&) {                   \
    caught = true;                       \
  }                                      \
  assert (caught);                       \
} while(0)


//************************************************************************


// ecal barrel
void test1()
{
  std::vector<mapkey_t> ids = make_ecalb_ids();

  using Indexer_t = k4::recCalo::IDMapIndexer<3>;
  using IDMap_t = Indexer_t::IDMap_t;
  using FieldDesc_t = Indexer_t::FieldDesc_t;
  dd4hep::IDDescriptor desc ("desc", ecalb_descstr);

  std::vector<FieldDesc_t> fieldDescs
    {
      IDMap_t::makeDesc (*desc.field ("layer")),
      IDMap_t::makeDesc (*desc.field ("theta")),
      IDMap_t::makeDesc (*desc.field ("module")),
    };

  Indexer_t map (4, 6, fieldDescs, ids);
  assert (map.detIDs().size() == 1);
  assert (map.detIDs()[0] == 4);
  assert (map.detIDBits() == 6);
  assert (map.byteSize() > 100000);

  size_t ncell = ids.size();
  assert (map.cellIDs().size() == ncell);
  for (size_t i = 0; i < ncell; i++) {
    assert (map.cellIDs()[i] == ids[i]);
    assert (map.index(ids[i]) == i);
  }
  assert (map.index (0) == Indexer_t::INVALID);
  assert (map.index (ids[0]+1) == Indexer_t::INVALID);

  std::vector<mapkey_t> ids2 (ids.begin(), ids.end());
  ids2[10] += 1;
  EXPECT_EXCEPTION( std::runtime_error, Indexer_t map2 (4, 6, fieldDescs, ids2) );
}


// ecal endcap.  Here we need to ignore the layers field.
void test2()
{
  std::vector<mapkey_t> ids = make_ecale_ids();

  using Indexer_t = k4::recCalo::IDMapIndexer<4>;
  using IDMap_t = Indexer_t::IDMap_t;
  using FieldDesc_t = Indexer_t::FieldDesc_t;
  dd4hep::IDDescriptor desc ("desc", ecale_descstr);

  const dd4hep::BitFieldElement* bfe_side = desc.field ("side");
  const dd4hep::BitFieldElement* bfe_wheel = desc.field ("wheel");

  std::vector<FieldDesc_t> fieldDescs
    { FieldDesc_t (std::min(bfe_side->offset(), bfe_wheel->offset()),
                   bfe_side->width() + bfe_wheel->width()),
      IDMap_t::makeDesc (*desc.field("rho")),
      IDMap_t::makeDesc (*desc.field("z")),
      IDMap_t::makeDesc (*desc.field("module"))
    };
  std::vector<FieldDesc_t> ignoredDescs
    { IDMap_t::makeDesc (*desc.field("layer"))
    };

  Indexer_t map (4, 6, fieldDescs, ids, 0, ignoredDescs);
  assert (map.detIDs().size() == 1);
  assert (map.detIDs()[0] == 4);
  assert (map.detIDBits() == 6);
  assert (map.byteSize() > 100000);

  size_t ncell = ids.size();
  assert (map.cellIDs().size() == ncell);
  for (size_t i = 0; i < ncell; i++) {
    assert (map.cellIDs()[i] == ids[i]);
    assert (map.index(ids[i]) == i);
  }
}


int main()
{
  test1();
  test2();
  return 0;
}
