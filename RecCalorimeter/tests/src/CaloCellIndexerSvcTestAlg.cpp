/**
 * @file RecCalorimeter/tests/src/CaloCellIndexerSvcTestAlg.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date Jan, 2026
 * @brief Test for CaloCellIndexerSvc
 */

// clang-format off: headers best listed in reverse-dependency order.
#include "RecCaloCommon/ICaloCellIndexerSvc.h"
#include "k4FWCore/GaudiChecks.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/ServiceHandle.h"
#include <span>
// clang-format on

namespace k4::recCalo {

class CaloCellIndexerSvcTestAlg : public Algorithm {
public:
  using Algorithm::Algorithm;

  virtual StatusCode initialize() override;
  virtual StatusCode execute() override;

private:
  // clang-format off
  ServiceHandle<ICaloCellIndexerSvc> m_svc
  { this, "CaloCellIndexerSvc", "k4::recCalo::CaloCellIndexerSvc", "" };
  // clang-format on

  Gaudi::Property<std::vector<int>> m_detIDs{this, "DetIDs", {4}, ""};
};

DECLARE_COMPONENT(k4::recCalo::CaloCellIndexerSvcTestAlg);

StatusCode CaloCellIndexerSvcTestAlg::initialize() {
  K4_GAUDI_CHECK(m_svc.retrieve());

  K4_GAUDI_CHECK(m_svc->indexer(999) == nullptr);
  K4_GAUDI_CHECK(m_detIDs.size() == 2);

  const ICaloIndexer* indexer0 = m_svc->indexer(m_detIDs[0]);
  K4_GAUDI_CHECK(indexer0 != nullptr);
  K4_GAUDI_CHECK(indexer0->detIDs().size() == 1);
  K4_GAUDI_CHECK(indexer0->detIDs()[0] == m_detIDs[0]);

  std::span<const uint64_t> ids0 = indexer0->cellIDs();
  for (size_t i = 0; i < ids0.size(); ++i) {
    K4_GAUDI_CHECK(indexer0->index(ids0[i]) == i);
  }

  const ICaloIndexer* indexer1 = m_svc->indexer(m_detIDs[1]);
  K4_GAUDI_CHECK(indexer1 != nullptr);
  K4_GAUDI_CHECK(indexer1->detIDs().size() == 1);
  K4_GAUDI_CHECK(indexer1->detIDs()[0] == m_detIDs[1]);

  std::span<const uint64_t> ids1 = indexer1->cellIDs();
  for (size_t i = 0; i < ids1.size(); ++i) {
    K4_GAUDI_CHECK(indexer1->index(ids1[i]) == i);
  }

  const ICaloIndexer* indexer01 = m_svc->indexer(m_detIDs);
  K4_GAUDI_CHECK(indexer01 != nullptr);
  K4_GAUDI_CHECK(indexer01->detIDs().size() == 2);
  K4_GAUDI_CHECK(indexer01->detIDs()[0] == m_detIDs[0]);
  K4_GAUDI_CHECK(indexer01->detIDs()[1] == m_detIDs[1]);

  std::span<const uint64_t> ids01 = indexer01->cellIDs();
  K4_GAUDI_CHECK(ids01.size() == ids0.size() + ids1.size());
  K4_GAUDI_CHECK(std::equal(ids0.begin(), ids0.end(), ids01.begin()));
  K4_GAUDI_CHECK(std::equal(ids1.begin(), ids1.end(), ids01.begin() + ids0.size()));
  for (size_t i = 0; i < ids01.size(); ++i) {
    K4_GAUDI_CHECK(indexer01->index(ids01[i]) == i);
  }

  return StatusCode::SUCCESS;
}

StatusCode CaloCellIndexerSvcTestAlg::execute() { return StatusCode::SUCCESS; }

} // namespace k4::recCalo
