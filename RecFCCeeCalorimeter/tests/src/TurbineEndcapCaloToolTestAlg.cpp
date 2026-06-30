/**
 * @file RecFCCeeCalorimeter/tests/src/TurbineEndcapCaloToolTestAlg.cpp
 * @author Giovanni Marchiori <giovanni.marchiori@cern.ch>
 * @date June, 2026
 * @brief Test for TurbineEndcapCaloToolTest
 */

#undef NDEBUG
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "RecCaloCommon/ICalorimeterTool.h"
#include "k4FWCore/GaudiChecks.h"
#include <cassert>
#include <span>

namespace k4::recCalo {

class TurbineEndcapCaloToolTestAlg : public Algorithm {
  using Algorithm::Algorithm;

  virtual StatusCode initialize() override;
  virtual StatusCode execute() override;

private:
  ToolHandle<k4::recCalo::ICalorimeterTool> m_tool{this, "CalorimeterTool", "TurbineEndcapCaloTool", ""};
};

DECLARE_COMPONENT(k4::recCalo::TurbineEndcapCaloToolTestAlg);

StatusCode TurbineEndcapCaloToolTestAlg::initialize() {
  K4_GAUDI_CHECK(m_tool.retrieve());
  K4_GAUDI_CHECK(m_tool->readoutName() == "ECalEndcapTurbine");
  K4_GAUDI_CHECK(m_tool->id() == 5);

  std::span<const uint64_t> ids = m_tool->cellIDs();
  size_t ncells = ids.size();
  K4_GAUDI_CHECK(ncells == 1203200);

  std::unique_ptr<ICaloIndexer> indexer = m_tool->indexer();
  K4_GAUDI_CHECK(indexer != nullptr);
  K4_GAUDI_CHECK(indexer->detIDBits() == 4);

  K4_GAUDI_CHECK(indexer->cellIDs().size() == ncells);
  for (size_t i = 0; i < ncells; i++) {
    K4_GAUDI_CHECK(ids[i] == indexer->cellIDs()[i]);
    K4_GAUDI_CHECK(indexer->index(ids[i]) == i);
  }

  return StatusCode::SUCCESS;
}

StatusCode TurbineEndcapCaloToolTestAlg::execute() { return StatusCode::SUCCESS; }

} // namespace k4::recCalo
