#include "ReadCaloCrosstalkMap.h"

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

DECLARE_COMPONENT(ReadCaloCrosstalkMap)

ReadCaloCrosstalkMap::ReadCaloCrosstalkMap(const std::string& type, const std::string& name,
                                                 const IInterface* parent)
    : GaudiTool(type, name, parent) {
  declareInterface<ICaloReadCrosstalkMap>(this);
}

StatusCode ReadCaloCrosstalkMap::initialize() {
  // prevent to initialize the tool if not intended (input file path empty)
  // otherwise things will crash if m_fileName is not available
  // not a perfect solution but tools seems to not be meant to be optional
  if (m_fileName == "") return StatusCode::SUCCESS;

  StatusCode sc = GaudiTool::initialize();
  info() <<"Loading crosstalk map..." << endmsg;
  if (sc.isFailure()) return sc;
  std::unique_ptr<TFile> file(TFile::Open(m_fileName.value().c_str(),"READ"));
  TTree* tree = nullptr;
  file->GetObject("crosstalk_neighbours",tree);
  ULong64_t read_cellId;
  std::vector<uint64_t>  *read_neighbours=0;
  std::vector<double> *read_crosstalks=0;
  
  tree->SetBranchAddress("cellId",&read_cellId);
  tree->SetBranchAddress("list_crosstalk_neighbours", &read_neighbours);
  tree->SetBranchAddress("list_crosstalks", &read_crosstalks);
  for (uint i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);
      m_mapNeighbours.insert(std::pair<uint64_t, std::vector<uint64_t>>(read_cellId, *read_neighbours));
      m_mapCrosstalks.insert(std::pair<uint64_t, std::vector<double>>(read_cellId, *read_crosstalks));
  }
  
  info() <<"Crosstalk input: " << m_fileName.value().c_str() << endmsg;
  info() << "Total number of cells = " << tree->GetEntries() << ", Size of crosstalk neighbours = " << m_mapNeighbours.size() << ", Size of coefficients = " << m_mapCrosstalks.size() << endmsg;
  delete tree;
  delete read_neighbours;
  delete read_crosstalks;
  file->Close();
  return sc;
}

StatusCode ReadCaloCrosstalkMap::finalize() { return GaudiTool::finalize(); }

std::vector<uint64_t>& ReadCaloCrosstalkMap::getNeighbours(uint64_t aCellId) {
  return m_mapNeighbours[aCellId];
}

std::vector<double>& ReadCaloCrosstalkMap::getCrosstalks(uint64_t aCellId) {
  return m_mapCrosstalks[aCellId];
}
