#include "TopoCaloNoisyCells.h"
#include "RecCaloCommon/k4RecCalorimeter_check.h"

#include "TBranch.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

DECLARE_COMPONENT(TopoCaloNoisyCells)

StatusCode TopoCaloNoisyCells::initialize() {
  K4RECCALORIMETER_CHECK( AlgTool::initialize() );

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

  TTree* tree = nullptr;
  inFile->GetObject("noisyCells", tree);
  ULong64_t readCellId;
  double readNoisyCells;
  double readNoisyCellsOffset;
  tree->SetBranchAddress("cellId", &readCellId);
  tree->SetBranchAddress("noiseLevel",
                         &readNoisyCells); // would be better to call branch noiseRMS rather than noiseLevel
  tree->SetBranchAddress("noiseOffset", &readNoisyCellsOffset);
  for (uint i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    m_map.insert(std::pair<uint64_t, std::pair<double, double>>(readCellId,
                                                                std::make_pair(readNoisyCells, readNoisyCellsOffset)));
  }
  delete tree;
  inFile->Close();

  return StatusCode::SUCCESS;
}

double TopoCaloNoisyCells::getNoiseRMSPerCell(uint64_t aCellId) const
{
  auto it = m_map.find(aCellId);
  if (it != m_map.end()) {
    return it->second.first;
  }
  return 0;
}


double TopoCaloNoisyCells::getNoiseOffsetPerCell(uint64_t aCellId) const
{
  auto it = m_map.find(aCellId);
  if (it != m_map.end()) {
    return it->second.second;
  }
  return 0;
}
