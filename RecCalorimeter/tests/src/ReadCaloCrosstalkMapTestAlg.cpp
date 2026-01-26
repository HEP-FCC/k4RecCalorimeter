/**
 * @file ReadCaloCrosstalkMapTestAlg.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date Jul, 2026
 * @brief Test for ReadCaloCrosstalkMap.
 */

#include "RecCaloCommon/ICaloCellIndexerSvc.h"
#include "RecCaloCommon/ICaloReadCrosstalkMap.h"

#include "k4FWCore/GaudiChecks.h"

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"

#include "TFile.h"
#include "TTree.h"
#include <memory>

namespace k4::recCalo {

class ReadCaloCrosstalkMapTestAlg : public Algorithm {
public:
  using CellID = k4::recCalo::ICaloIndexer::CellID;
  using Algorithm::Algorithm;

  virtual StatusCode initialize() override;
  virtual StatusCode execute() override;

private:
  ToolHandle<ICaloReadCrosstalkMap> m_tool{this, "Tool", "ReadCaloCrosstalkMap", ""};

  /// Cell indexing service.
  ServiceHandle<k4::recCalo::ICaloCellIndexerSvc> m_indexerSvc{
      this, "CaloCellIndexerSvc", "k4::recCalo::CaloCellIndexerSvc", "The cell indexing service."};

  Gaudi::Property<std::string> m_fileName{this, "fileName", "crosstalkTest.root", ""};

  Gaudi::Property<int> m_detID{this, "detID", -1, ""};

  Gaudi::Property<unsigned> m_ncells{this, "ncells", 100, ""};

  StatusCode writeFile(unsigned int ncells, const std::string& filename,
                       const k4::recCalo::ICaloIndexer& indexer) const;
  StatusCode checkCrosstalk(unsigned int ncells, const ICaloReadCrosstalkMap& tool,
                            const k4::recCalo::ICaloIndexer& indexer) const;
};

DECLARE_COMPONENT(k4::recCalo::ReadCaloCrosstalkMapTestAlg);

StatusCode ReadCaloCrosstalkMapTestAlg::initialize() {
  K4_GAUDI_CHECK(m_indexerSvc.retrieve());

  const k4::recCalo::ICaloIndexer* indexer = m_indexerSvc->indexer(m_detID);
  K4_GAUDI_CHECK(indexer != nullptr);
  K4_GAUDI_CHECK(writeFile(m_ncells, m_fileName, *indexer));
  K4_GAUDI_CHECK(m_tool.retrieve());
  K4_GAUDI_CHECK(checkCrosstalk(m_ncells, *m_tool, *indexer));

  return StatusCode::SUCCESS;
}

StatusCode ReadCaloCrosstalkMapTestAlg::execute() { return StatusCode::SUCCESS; }

StatusCode ReadCaloCrosstalkMapTestAlg::writeFile(unsigned int ncells, const std::string& fileName,
                                                  const k4::recCalo::ICaloIndexer& indexer) const {
  std::unique_ptr<TFile> xtalkFile(TFile::Open(fileName.c_str(), "RECREATE"));
  auto tree = std::make_unique<TTree>("crosstalk_neighbours", "crosstalk_neighbours");

  ULong64_t write_cellId = 0;
  std::vector<uint64_t> write_neighbours;
  std::vector<double> write_crosstalks;
  tree->Branch("cellId", &write_cellId);
  tree->Branch("list_crosstalk_neighbours", &write_neighbours);
  tree->Branch("list_crosstalks", &write_crosstalks);

  std::span<const CellID> cellIDs = indexer.cellIDs();
  for (size_t icell = 0; icell < ncells; icell++) {
    write_cellId = cellIDs[icell];
    unsigned nneigh = 6 + (icell % 4);
    unsigned ineigh = icell * 10 + 100;
    if (ineigh + nneigh >= cellIDs.size()) {
      error() << "Too many cells requested." << endmsg;
      return StatusCode::FAILURE;
    }

    write_neighbours.assign(cellIDs.begin() + ineigh, cellIDs.begin() + ineigh + nneigh);
    write_crosstalks.clear();
    for (size_t i = 0; i < nneigh; i++) {
      write_crosstalks.push_back(ineigh + i + 0.5);
    }
    tree->Fill();
  }

  tree->Write();
  xtalkFile->Close();
  tree.release();
  return StatusCode::SUCCESS;
}

StatusCode ReadCaloCrosstalkMapTestAlg::checkCrosstalk(unsigned int ncells, const ICaloReadCrosstalkMap& tool,
                                                       const k4::recCalo::ICaloIndexer& indexer) const {
  std::span<const CellID> cellIDs = indexer.cellIDs();
  for (size_t icell = 0; icell < ncells; icell++) {
    std::span<const CellID> neighs = tool.getNeighbours(cellIDs[icell]);
    std::span<const double> xtalk = tool.getCrosstalks(cellIDs[icell]);

    unsigned nneigh = 6 + (icell % 4);
    if (neighs.size() != nneigh || xtalk.size() != nneigh) {
      error() << "Size mismatch for cell index " << icell << " ID " << cellIDs[icell] << " expected " << nneigh
              << " got " << neighs.size() << " " << xtalk.size() << endmsg;
      return StatusCode::FAILURE;
    }
    unsigned ineigh = icell * 10 + 100;
    for (unsigned i = 0; i < nneigh; i++) {
      if (neighs[i] != cellIDs[ineigh + i]) {
        error() << "Neighbour mismatch for cell index " << icell << " ID " << cellIDs[icell] << " offset " << i
                << " expected " << cellIDs[ineigh + i] << " got " << neighs[i] << endmsg;
        return StatusCode::FAILURE;
      }
      if (xtalk[i] != ineigh + i + 0.5) {
        error() << "Crosstalk mismatch for cell index " << icell << " ID " << cellIDs[icell] << " offset " << i
                << " expected " << ineigh + i + 0.5 << " got " << xtalk[i] << endmsg;
        return StatusCode::FAILURE;
      }
    }
  }
  return StatusCode::SUCCESS;
}

} // namespace k4::recCalo
