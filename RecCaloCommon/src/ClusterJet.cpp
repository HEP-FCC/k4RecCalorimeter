#include "RecCaloCommon/ClusterJet.h"

// std
#include <iostream>

namespace k4::recCalo {

ClusterJet::ClusterJet(const std::string& jetAlg, double jetRadius, int isExclusiveClustering, double minPt, int njets)
    : m_jetAlg(jetAlg), m_jetRadius(jetRadius), m_isExclusiveClustering(isExclusiveClustering), m_minPt(minPt),
      m_njets(njets) {}

bool ClusterJet::initialize() {
  if (m_jetAlgMap.find(m_jetAlg) == m_jetAlgMap.end()) {
    std::cout << "ERROR: "
              << " is not in the list of supported jet algorithms" << std::endl;
    return false;
  }

  m_jetDef = new fastjet::JetDefinition(m_jetAlgMap.at(m_jetAlg), m_jetRadius);

  if (m_isExclusiveClustering > 1) {
    std::cout << "ERROR: "
              << "Clustering of " << m_isExclusiveClustering << " is currently not supported" << std::endl;
    ;
    return false;
  }

  return true;
}

std::vector<fastjet::PseudoJet> ClusterJet::cluster(const std::vector<fastjet::PseudoJet>& clustersPJ) {
  // Deleting cluster sequence from previous event
  delete m_clustSeq;

  std::vector<fastjet::PseudoJet> jets;

  m_clustSeq = new fastjet::ClusterSequence(clustersPJ, *m_jetDef);

  // Note: initialize has already checked if m_isExclusiveClustering has the
  // right range
  if (m_isExclusiveClustering == 0) {
    jets = fastjet::sorted_by_pt(m_clustSeq->inclusive_jets(m_minPt));
  } else if (m_isExclusiveClustering == 1) {
    jets = fastjet::sorted_by_pt(m_clustSeq->exclusive_jets(m_njets));
  }

  return jets;
}

} /* namespace k4::recCalo */
