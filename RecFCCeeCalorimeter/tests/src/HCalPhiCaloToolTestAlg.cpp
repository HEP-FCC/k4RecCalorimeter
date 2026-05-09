/**
 * @file RecFCCeeCalorimeter/tests/src/HCalPhiCaloToolTestAlg.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date May, 2026
 * @brief Test for HCalPhiThetaCaloTool/HCalPhiRowCaloTool
 */

#undef NDEBUG
// Headers should really be in reverse-dependency order, but clang-format
// complains if i do that...
#include "DD4hep/Detector.h"
#include "DD4hep/Readout.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"
#include "RecCaloCommon/ICalorimeterTool.h"
#include "detectorSegmentations/FCCSWHCalPhiRow_k4geo.h"
#include "k4FWCore/GaudiChecks.h"
#include "k4Interface/IGeoSvc.h"
#include <cassert>
#include <fstream>
#include <span>

namespace k4::recCalo {

class HCalPhiCaloToolTestAlg : public Algorithm {
  using Algorithm::Algorithm;

  virtual StatusCode initialize() override;
  virtual StatusCode execute() override;

private:
  StatusCode testTool(const k4::recCalo::ICalorimeterTool& tool, size_t exp_ncells) const;

  Gaudi::Property<int> m_expectedBarrelCells{this, "ExpectedBarrelCells", 0, ""};
  Gaudi::Property<int> m_expectedEndcapCells{this, "ExpectedEndcapCells", 0, ""};

  Gaudi::Property<std::string> m_expectedBarrelReadout{this, "ExpectedBarrelReadout", "", ""};
  Gaudi::Property<std::string> m_expectedEndcapReadout{this, "ExpectedEndcapReadout", "", ""};

  Gaudi::Property<bool> m_groupRows{this, "GroupRows", false, ""};
  Gaudi::Property<bool> m_dumpCells{this, "DumpCells", false, ""};

  ToolHandle<k4::recCalo::ICalorimeterTool> m_barrelTool{this, "HCalBarrelTool", "", ""};
  ToolHandle<k4::recCalo::ICalorimeterTool> m_endcapTool{this, "HCalEndcapTool", "", ""};

  /// Handle to the geometry service.
  ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc"};
};

DECLARE_COMPONENT(k4::recCalo::HCalPhiCaloToolTestAlg);

StatusCode HCalPhiCaloToolTestAlg::initialize() {
  if (m_groupRows) {
    K4_GAUDI_CHECK(m_geoSvc.retrieve());
    const dd4hep::Detector* det = m_geoSvc->getDetector();
    auto it = det->readouts().find("HCalEndcapReadoutPhiRow");
    if (it == det->readouts().end()) {
      error() << "Readout HCalEndcapReadoutPhiRow does not exist." << endmsg;
      return StatusCode::FAILURE;
    }
    dd4hep::Readout readout = it->second;
    auto* seg = dynamic_cast<dd4hep::DDSegmentation::FCCSWHCalPhiRow_k4geo*>(readout.segmentation().segmentation());
    if (!seg) {
      error() << "Unable to cast segmentation pointer to FCCSWHCalPhiRow_k4geo!!!!" << endmsg;
      return StatusCode::FAILURE;
    }
    // This is three groups of 6, but clang-format won't let me write it
    // in a way that shows the actual structure.
    seg->setGroupedRows(std::vector<int>{3, 3, 3, 6, 6, 6, 6, 6, 6, 9, 0, 0, 10, 10, 10, 15, 15, 25});
  }

  K4_GAUDI_CHECK(m_barrelTool.retrieve());
  K4_GAUDI_CHECK(m_endcapTool.retrieve());

  K4_GAUDI_CHECK(m_barrelTool->readoutName() == m_expectedBarrelReadout);
  K4_GAUDI_CHECK(m_barrelTool->id() == 8);
  K4_GAUDI_CHECK(m_endcapTool->readoutName() == m_expectedEndcapReadout);
  K4_GAUDI_CHECK(m_endcapTool->id() == 9);

  K4_GAUDI_CHECK(testTool(*m_barrelTool, m_expectedBarrelCells));
  K4_GAUDI_CHECK(testTool(*m_endcapTool, m_expectedEndcapCells));

  return StatusCode::SUCCESS;
}

StatusCode HCalPhiCaloToolTestAlg::testTool(const k4::recCalo::ICalorimeterTool& tool, size_t exp_ncells) const {
  std::span<const uint64_t> ids = tool.cellIDs();
  size_t ncells = ids.size();
  K4_GAUDI_CHECK(ncells == exp_ncells);

  std::unique_ptr<ICaloIndexer> indexer = tool.indexer();
  K4_GAUDI_CHECK(indexer != nullptr);
  K4_GAUDI_CHECK(indexer->detIDBits() == 4);

  K4_GAUDI_CHECK(indexer->cellIDs().size() == ncells);
  for (size_t i = 0; i < ncells; i++) {
    K4_GAUDI_CHECK(ids[i] == indexer->cellIDs()[i]);
    K4_GAUDI_CHECK(indexer->index(ids[i]) == i);
  }

  if (m_dumpCells) {
    std::string fname = tool.readoutName() + ".cells";
    std::ofstream os(fname);
    for (uint64_t id : ids) {
      os << id << "\n";
    }
    os.close();
  }

  return StatusCode::SUCCESS;
}

StatusCode HCalPhiCaloToolTestAlg::execute() { return StatusCode::SUCCESS; }

} // namespace k4::recCalo
