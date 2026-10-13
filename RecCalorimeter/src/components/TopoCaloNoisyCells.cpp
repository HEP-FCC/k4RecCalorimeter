#include "TopoCaloNoisyCells.h"
#include "DDSegmentation/BitFieldCoder.h"
#include "k4FWCore/GaudiChecks.h"

#include "TBranch.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

#include <algorithm>

DECLARE_COMPONENT(TopoCaloNoisyCells)

StatusCode TopoCaloNoisyCells::initialize() {
  K4_GAUDI_CHECK(AlgTool::initialize());
  K4_GAUDI_CHECK(m_constantsSvc.retrieve());
  K4_GAUDI_CHECK(m_indexerSvc.retrieve());

  // setup system decoder
  m_decoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(m_systemEncoding);
  m_indexSystem = m_decoder->index("system");

  m_data = m_constantsSvc->getObj<NoiseData>(m_fileName);
  if (!m_data) {
    // Check if file exists
    if (m_fileName.empty()) {
      error() << "Name of the file with the noisy cells not provided!" << endmsg;
      return StatusCode::FAILURE;
    }
    if (gSystem->AccessPathName(m_fileName.value().c_str())) {
      error() << "Provided file with the noisy cells not found!" << endmsg;
      error() << "File path: " << m_fileName.value() << endmsg;
      return StatusCode::FAILURE;
    }
    std::unique_ptr<TFile> inFile(TFile::Open(m_fileName.value().c_str(), "READ"));
    if (inFile->IsZombie()) {
      error() << "Unable to open the file with the noisy cells!" << endmsg;
      error() << "File path: " << m_fileName.value() << endmsg;
      return StatusCode::FAILURE;
    } else {
      info() << "Using the following file with the noisy cells: " << m_fileName.value() << endmsg;
    }

    NoiseData data = readData(*inFile);
    K4_GAUDI_CHECK(m_constantsSvc->putObj(m_fileName, std::move(data)));
    m_data = m_constantsSvc->getObj<NoiseData>(m_fileName);
    K4_GAUDI_CHECK(m_data != nullptr);
  }

  m_indexer = m_data->m_indexer;
  K4_GAUDI_CHECK(m_indexer != nullptr);

  return StatusCode::SUCCESS;
}

auto TopoCaloNoisyCells::readData(TFile& inFile) const -> NoiseData {
  NoiseData data;

  TTree* tree = nullptr;
  inFile.GetObject("noisyCells", tree);
  ULong64_t readCellId;
  double readNoisyCells;
  double readNoisyCellsOffset;
  tree->SetBranchAddress("cellId", &readCellId);
  tree->SetBranchAddress("noiseLevel",
                         &readNoisyCells); // would be better to call branch noiseRMS rather than noiseLevel
  tree->SetBranchAddress("noiseOffset", &readNoisyCellsOffset);

  // First find the set of detIDs in order to get the proper indexer.
  // Just use a vector, since we only expect a handful.
  std::vector<int> detIDs;
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("cellId", 1);
  for (uint i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    int detID = m_decoder->get(readCellId, m_indexSystem);
    if (std::ranges::find(detIDs, detID) == detIDs.end()) {
      detIDs.push_back(detID);
    }
  }
  std::ranges::sort(detIDs);
  data.m_indexer = m_indexerSvc->indexer(detIDs);
  if (!data.m_indexer)
    return data;

  data.m_noise.resize(data.m_indexer->cellIDs().size());

  tree->SetBranchStatus("*", 1);
  for (uint i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    unsigned ndx = data.m_indexer->index(readCellId);
    data.m_noise.at(ndx) = std::make_pair(readNoisyCells, readNoisyCellsOffset);
  }
  delete tree;
  inFile.Close();

  return data;
}

double TopoCaloNoisyCells::getNoiseRMSPerCell(CellID aCellId) const {
  unsigned ndx = m_indexer->index(aCellId);
  return m_data->m_noise.at(ndx).first;
}

double TopoCaloNoisyCells::getNoiseOffsetPerCell(CellID aCellId) const {
  unsigned ndx = m_indexer->index(aCellId);
  return m_data->m_noise.at(ndx).second;
}

std::pair<double, double> TopoCaloNoisyCells::getNoisePerCell(CellID aCellId) const {
  unsigned ndx = m_indexer->index(aCellId);
  return m_data->m_noise.at(ndx);
}
