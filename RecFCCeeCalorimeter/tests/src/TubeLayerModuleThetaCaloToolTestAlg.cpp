/**
 * @file RecFCCeeCalorimeter/tests/src/TubeLayerModuleThetaCaloToolTestAlg.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date Jan, 2026
 * @brief Test for TubeLayerModuleThetaCaloToolTest
 */

#undef NDEBUG
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "RecCaloCommon/ICalorimeterTool.h"
#include "k4FWCore/GaudiChecks.h"
#include <cassert>
#include <span>

namespace k4::recCalo {

class TubeLayerModuleThetaCaloToolTestAlg : public Algorithm {
  using Algorithm::Algorithm;

  virtual StatusCode initialize() override;
  virtual StatusCode execute() override;

private:
  ToolHandle<k4::recCalo::ICalorimeterTool> m_tool{this, "CalorimeterTool", "TubeLayerModuleThetaCaloTool", ""};
};

DECLARE_COMPONENT(k4::recCalo::TubeLayerModuleThetaCaloToolTestAlg);

StatusCode TubeLayerModuleThetaCaloToolTestAlg::initialize() {
  K4_GAUDI_CHECK(m_tool.retrieve());
  K4_GAUDI_CHECK(m_tool->readoutName() == "ECalBarrelModuleThetaMerged");
  K4_GAUDI_CHECK(m_tool->id() == 4);

  std::span<const uint64_t> ids = m_tool->cellIDs();
  size_t ncells = ids.size();
  K4_GAUDI_CHECK(ncells == 2041344);

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

StatusCode TubeLayerModuleThetaCaloToolTestAlg::execute() { return StatusCode::SUCCESS; }

} // namespace k4::recCalo
