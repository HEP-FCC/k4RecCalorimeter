#include "TopoCaloNeighbours.h"
#include "k4FWCore/GaudiChecks.h"

#include "TBranch.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

DECLARE_COMPONENT(TopoCaloNeighbours)

StatusCode TopoCaloNeighbours::initialize() {
  K4_GAUDI_CHECK( base_class::initialize() );

  // Check if neighbours map file exists
  if (m_fileName.empty()) {
    error() << "Proper filepath for the neighbours map not provided!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (gSystem->AccessPathName(m_fileName.value().c_str())) {
    error() << "Provided neighbours map file not found!" << endmsg;
    error() << "File path: " << m_fileName.value() << endmsg;
    return StatusCode::FAILURE;
  }
  std::unique_ptr<TFile> inFile(TFile::Open(m_fileName.value().c_str(), "READ"));
  if (inFile->IsZombie()) {
    error() << "Unable to open the provide file with neighbours map!" << endmsg;
    error() << "File path: " << m_fileName.value() << endmsg;
    return StatusCode::FAILURE;
  } else {
    info() << "Using the following file with neighbours map: " << m_fileName.value() << endmsg;
  }

  TTree* tree = nullptr;
  inFile->GetObject("neighbours", tree);
  ULong64_t readCellId;
  std::vector<uint64_t>* readNeighbours = nullptr;
  tree->SetBranchAddress("cellId", &readCellId);
  tree->SetBranchAddress("neighbours", &readNeighbours);
  for (uint i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    m_map.insert(std::pair<uint64_t, std::vector<uint64_t>>(readCellId, *readNeighbours));
  }
  std::vector<int> counterL;
  counterL.assign(100, 0);
  for (const auto& item : m_map) {
    counterL[item.second.size()]++;
  }
  for (uint iCount = 0; iCount < counterL.size(); iCount++) {
    if (counterL[iCount] != 0) {
      info() << counterL[iCount] << " cells have " << iCount << " neighbours" << endmsg;
    }
  }
  delete tree;
  delete readNeighbours;
  inFile->Close();

  return StatusCode::SUCCESS;
}

auto TopoCaloNeighbours::neighbours(CellID aCellId) const -> const std::vector<CellID>&
{
  auto it = m_map.find(aCellId);
  if (it != m_map.end())
    return it->second;
  static const std::vector<CellID> empty;
  return empty;
}
