//#include "k4RecCalorimeterPlugins/ClusterJet.h"
//#include "include/RecCalorimeter/ClusterJet.h"
#include "ClusterJet.h"

// std
#include <vector>
#include <math.h>

namespace k4::recCalo{

ClusterJet::ClusterJet(std::string jetAlg, double jetRadius, int clustering, double minPt, int njets): m_jetAlg(jetAlg), m_jetRadius(jetRadius), m_clustering(clustering), m_minPt(minPt), m_njets(njets){
  m_clustSeq = nullptr;
}


StatusCode ClusterJet::initialize(){
  if (m_jetAlgMap.find(m_jetAlg) == m_jetAlgMap.end()) {
    std::cout << "ERROR: " << " is not in the list of supported jet algorithms" << std::endl;;
    return StatusCode::FAILURE;
  }
  if(m_clustering > 1){
    std::cout << "ERROR: " << "Clustering of " << m_clustering << " is currently not supported" << std::endl;;
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}

//std::pair<std::vector<fastjet::PseudoJet>, ClusterSequence*> ClusterJet::cluster(const   std::vector<fastjet::PseudoJet> clustersPJ) const{
std::vector<fastjet::PseudoJet> ClusterJet::cluster(const   std::vector<fastjet::PseudoJet> clustersPJ) {
  std::vector <fastjet::PseudoJet> jets;

  fastjet::JetDefinition* jetDef = new fastjet::JetDefinition( m_jetAlgMap.at(m_jetAlg), m_jetRadius);
  //fastjet::ClusterSequence* clustSeq = new fastjet::ClusterSequence(clustersPJ, *jetDef);
  if(m_clustSeq) delete m_clustSeq;
  m_clustSeq = new fastjet::ClusterSequence(clustersPJ, *jetDef);

  // Note: initialize has already checked if m_clustering has the right range
  if(m_clustering == 0){
    jets = fastjet::sorted_by_pt(m_clustSeq->inclusive_jets(m_minPt));
  }
  else if(m_clustering == 1){
    jets = fastjet::sorted_by_pt(m_clustSeq->exclusive_jets(m_njets));
  }

  //return std::make_pair(jets, clustSeq);
  return jets;

}


}

