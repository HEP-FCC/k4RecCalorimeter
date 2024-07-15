#include "TopoCaloNeighbours.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

DECLARE_COMPONENT(TopoCaloNeighbours)

TopoCaloNeighbours::TopoCaloNeighbours(const std::string& type, const std::string& name,
                                       const IInterface* parent)
    : GaudiTool(type, name, parent) {
  declareInterface<ICaloReadNeighboursMap>(this);
}

StatusCode TopoCaloNeighbours::initialize() {
  {
    StatusCode sc = GaudiTool::initialize();
    if (sc.isFailure()) return sc;
  }

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
    info() << "Using the following file with neighbours map: "
           << m_fileName.value() << endmsg;
  }

  TTree* tree = nullptr;
  inFile->GetObject("neighbours", tree);
  ULong64_t readCellId;
  std::vector<uint64_t>* readNeighbours = nullptr;
  tree->SetBranchAddress("cellId",&readCellId);
  tree->SetBranchAddress("neighbours",&readNeighbours);
  for (uint i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);
      m_map.insert(std::pair<uint64_t, std::vector<uint64_t>>(readCellId, *readNeighbours));
  }
  std::vector<int> counterL;
  counterL.assign(100,0);
  for(const auto& item: m_map) {
    counterL[item.second.size()] ++;
  }
  for(uint iCount = 0; iCount < counterL.size(); iCount++) {
    if (counterL[iCount] != 0) {
      info() << counterL[iCount] << " cells have " << iCount << " neighbours" << endmsg;
    }
  }
  delete tree;
  delete readNeighbours;
  inFile->Close();

  return StatusCode::SUCCESS;
}

StatusCode TopoCaloNeighbours::finalize() { return GaudiTool::finalize(); }

std::vector<uint64_t>& TopoCaloNeighbours::neighbours(uint64_t aCellId) {
  return m_map[aCellId];
}
