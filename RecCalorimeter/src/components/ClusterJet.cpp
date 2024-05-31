#include "ClusterJet.h"

// std
#include <vector>
#include <math.h>

ClusterJet::ClusterJet(std::string jetAlg, double jetRadius, double minPt): m_jetAlg(jetAlg), m_jetRadius(jetRadius), m_minPt(minPt){
}


std::vector<fastjet::PseudoJet> ClusterJet::cluster(const   std::vector<fastjet::PseudoJet> clustersPJ) const{
  if (m_jetAlgMap.find(m_jetAlg) == m_jetAlgMap.end()) {
    //error() << m_jetAlg << " is not in the list of supported jet algorithms" << endmsg;
    //return StatusCode::FAILURE;
    return std::vector <fastjet::PseudoJet>();
  }
  fastjet::JetDefinition* jetDef = new fastjet::JetDefinition( m_jetAlgMap.at(m_jetAlg), m_jetRadius);
  fastjet::ClusterSequence clustSeq(clustersPJ, *jetDef);
  std::vector <fastjet::PseudoJet> inclusiveJets = fastjet::sorted_by_pt(clustSeq.inclusive_jets(m_minPt));

  return inclusiveJets;

}



