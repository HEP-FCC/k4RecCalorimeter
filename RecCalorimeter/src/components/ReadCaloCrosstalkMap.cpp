#include "ReadCaloCrosstalkMap.h"

#include "TBranch.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

DECLARE_COMPONENT(ReadCaloCrosstalkMap)

ReadCaloCrosstalkMap::ReadCaloCrosstalkMap(const std::string& type, const std::string& name, const IInterface* parent)
    : AlgTool(type, name, parent) {
  declareInterface<ICaloReadCrosstalkMap>(this);
}

StatusCode ReadCaloCrosstalkMap::initialize() {
  // prevent to initialize the tool if not intended (input file path empty)
  // otherwise things will crash if m_fileName is not available
  // not a perfect solution but tools seems to not be meant to be optional
  if (m_fileName == "") {
    debug() << "Empty 'fileName' provided, it means cross-talk map is not needed, exitting ReadCaloCrosstalkMap "
               "initilization"
            << endmsg;
    return StatusCode::SUCCESS;
  }

  {
    StatusCode sc = AlgTool::initialize();
    info() << "Loading crosstalk map..." << endmsg;
    if (sc.isFailure())
      return sc;
  }

  // Check if crosstalk file exists
  if (gSystem->AccessPathName(m_fileName.value().c_str())) {
    error() << "Provided file with the crosstalk map not found!" << endmsg;
    error() << "File path: " << m_fileName.value() << endmsg;
    return StatusCode::FAILURE;
  }
  std::unique_ptr<TFile> xtalkFile(TFile::Open(m_fileName.value().c_str(), "READ"));
  if (xtalkFile->IsZombie()) {
    error() << "Unable to read the file with the crosstalk map!" << endmsg;
    error() << "File path: " << m_fileName.value() << endmsg;
    return StatusCode::FAILURE;
  } else {
    info() << "Using the following file with the crosstalk map: " << m_fileName.value() << endmsg;
  }

  TTree* tree = nullptr;
  xtalkFile->GetObject("crosstalk_neighbours", tree);
  ULong64_t read_cellId;
  std::vector<uint64_t>* read_neighbours = nullptr;
  std::vector<double>* read_crosstalks = nullptr;

  tree->SetBranchAddress("cellId", &read_cellId);
  tree->SetBranchAddress("list_crosstalk_neighbours", &read_neighbours);
  tree->SetBranchAddress("list_crosstalks", &read_crosstalks);
  for (uint i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    m_mapNeighbours.insert(std::pair<uint64_t, std::vector<uint64_t>>(read_cellId, *read_neighbours));
    m_mapCrosstalks.insert(std::pair<uint64_t, std::vector<double>>(read_cellId, *read_crosstalks));
  }

  info() << "Crosstalk input: " << m_fileName.value().c_str() << endmsg;
  info() << "Total number of cells = " << tree->GetEntries()
         << ", Size of crosstalk neighbours = " << m_mapNeighbours.size()
         << ", Size of coefficients = " << m_mapCrosstalks.size() << endmsg;
  delete tree;
  delete read_neighbours;
  delete read_crosstalks;
  xtalkFile->Close();

  return StatusCode::SUCCESS;
}

StatusCode ReadCaloCrosstalkMap::finalize() { return AlgTool::finalize(); }

std::vector<uint64_t>& ReadCaloCrosstalkMap::getNeighbours(uint64_t aCellId) { return m_mapNeighbours[aCellId]; }

const std::vector<uint64_t>&
ReadCaloCrosstalkMap::getNeighbours(uint64_t aCellId) const {
  auto it = m_mapNeighbours.find(aCellId);
  if (it != m_mapNeighbours.end()) {
    return it->second;
  }
  static const std::vector<uint64_t> empty;
  return empty;
}

std::vector<double>& ReadCaloCrosstalkMap::getCrosstalks(uint64_t aCellId) { return m_mapCrosstalks[aCellId]; }

const std::vector<double>&
ReadCaloCrosstalkMap::getCrosstalks(uint64_t aCellId) const {
  auto it = m_mapCrosstalks.find(aCellId);
  if (it != m_mapCrosstalks.end()) {
    return it->second;
  }
  static const std::vector<double> empty;
  return empty;
}
