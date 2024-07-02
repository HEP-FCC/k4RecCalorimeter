#include "ClusterJet.h"

// std
#include <vector>
#include <math.h>

namespace k4::recCalo{

ClusterJet::ClusterJet(std::string jetAlg, double jetRadius, int isExclusiveClustering, double minPt, int njets): m_jetAlg(jetAlg), m_jetRadius(jetRadius), m_isExclusiveClustering(isExclusiveClustering), m_minPt(minPt), m_njets(njets){
  m_clustSeq = nullptr;
}


StatusCode ClusterJet::initialize(){
  if (m_jetAlgMap.find(m_jetAlg) == m_jetAlgMap.end()) {
    std::cout << "ERROR: " << " is not in the list of supported jet algorithms" << std::endl;;
    return StatusCode::FAILURE;
  }
  if(m_isExclusiveClustering > 1){
    std::cout << "ERROR: " << "Clustering of " << m_isExclusiveClustering << " is currently not supported" << std::endl;;
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}

std::vector<fastjet::PseudoJet> ClusterJet::cluster(const   std::vector<fastjet::PseudoJet> clustersPJ) {
  std::vector <fastjet::PseudoJet> jets;

  fastjet::JetDefinition* jetDef = new fastjet::JetDefinition( m_jetAlgMap.at(m_jetAlg), m_jetRadius);
  if(m_clustSeq) delete m_clustSeq;
  m_clustSeq = new fastjet::ClusterSequence(clustersPJ, *jetDef);

  // Note: initialize has already checked if m_isExclusiveClustering has the right range
  if(m_isExclusiveClustering == 0){
    jets = fastjet::sorted_by_pt(m_clustSeq->inclusive_jets(m_minPt));
  }
  else if(m_isExclusiveClustering == 1){
    jets = fastjet::sorted_by_pt(m_clustSeq->exclusive_jets(m_njets));
  }

  return jets;

}


}

