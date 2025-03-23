#ifndef RECCALOCOMMON_CLUSTERJET_H
#define RECCALOCOMMON_CLUSTERJET_H

// std
#include <map>
#include <string>
#include <vector>

// FastJet
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"

namespace k4::recCalo {

class ClusterInfo : public fastjet::PseudoJet::UserInfoBase {
public:
  explicit ClusterInfo(const int& index) : _index(index) {}
  int index() const { return _index; }

protected:
  int _index = 0;
};

/** @class ClusterJet
 * k4RecCalorimeter/RecCalorimeter/src/components/ClusterJet.h
 *
 *  Helper class for reconstructing jets, used to avoid code duplication in the
 * classes that require specific input types. It takes a vector of
 * fastjet::PseudoJets as inputs, and returns the a vector of pseudojets after
 * clustering.
 *
 *  JetAlg: A string corresponding to the jet algorithm for clustering
 *  JetRadius: The radius of the jets being clustered
 *  MinPt: The pT threshold below which jets are ignored
 *  isExclusiveClustering: 1 if jets should use an exclusive clustering, 0
 * otherwise clusterArgs: Any other clustering arguments, such as the number of
 * jets for exclusive clustering. This is not used everywhere, so check the code
 * to make sure it is implemented for your use-case.
 *
 *  @author Jennifer Roloff
 *  @date 2024-7
 */

class ClusterJet {
public:
  ClusterJet(const std::string& jetAlg, double jetRadius, int isExclusiveClustering = 0, double minPt = 0,
             int clusterArgs = 0);
  bool initialize();
  ~ClusterJet() {
    delete m_jetDef;
    delete m_clustSeq;
  }

  std::vector<fastjet::PseudoJet> cluster(const std::vector<fastjet::PseudoJet>& clustersPJ);

private:
  std::map<std::string, fastjet::JetAlgorithm> m_jetAlgMap = {
      {"kt", fastjet::JetAlgorithm::kt_algorithm},         {"cambridge", fastjet::JetAlgorithm::cambridge_algorithm},
      {"antikt", fastjet::JetAlgorithm::antikt_algorithm}, {"genkt", fastjet::JetAlgorithm::genkt_algorithm},
      {"ee_kt", fastjet::JetAlgorithm::ee_kt_algorithm},   {"ee_genkt", fastjet::JetAlgorithm::ee_genkt_algorithm},
  };

  std::string m_jetAlg = "antikt";
  double m_jetRadius = 0.4;
  int m_isExclusiveClustering = 0; // Inclusive clustering by default
  double m_minPt = 10;             // Only relevant for inclusive clustering
  int m_njets = 0;                 // Only relevant for exclusive clustering

  fastjet::JetDefinition* m_jetDef = nullptr;
  fastjet::ClusterSequence* m_clustSeq = nullptr;
};

} /* namespace k4::recCalo */
#endif /* RECCALOCOMMON_CLUSTERJET_H */
