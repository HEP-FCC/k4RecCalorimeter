#ifndef RECCALORIMETER_CLUSTERJET_H
#define RECCALORIMETER_CLUSTERJET_H
//#ifndef k4RecCalorimeterPlugins_CLUSTERJET_H
//#define k4RecCalorimeterPlugins_CLUSTERJET_H

// std
#include <map>

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"

#include "fastjet/JetDefinition.hh"

// EDM4hep
#include "edm4hep/ClusterCollection.h"


namespace k4::recCalo{

class ClusterInfo : public fastjet::PseudoJet::UserInfoBase{
 public:
 ClusterInfo(const int & index) : _index(index){}
  int index() const{ return _index; }
 protected:
  int _index;
};



/** @class ClusterJet k4RecCalorimeter/RecCalorimeter/src/components/ClusterJet.h
 *
 *  @author Jennifer Roloff
 */

class ClusterJet {
public:
  ClusterJet(std::string jetAlg, double jetRadius, int clustering = 0, double minPt = 0, int clusterArgs = 0);
  StatusCode initialize();
  ~ClusterJet(){if(m_clustSeq) delete m_clustSeq;}

  std::vector<fastjet::PseudoJet> cluster(const   std::vector<fastjet::PseudoJet> clustersPJ);

private:

 

  std::map<std::string, fastjet::JetAlgorithm> m_jetAlgMap = {
                        {"kt",                    fastjet::JetAlgorithm::kt_algorithm},
                        {"cambridge",             fastjet::JetAlgorithm::cambridge_algorithm},
                        {"antikt",                fastjet::JetAlgorithm::antikt_algorithm},
                        {"genkt",                 fastjet::JetAlgorithm::genkt_algorithm},
                        {"ee_kt",                 fastjet::JetAlgorithm::ee_kt_algorithm},
                        {"ee_genkt",              fastjet::JetAlgorithm::ee_genkt_algorithm},
  };

  std::string m_jetAlg = "antikt";
  double m_jetRadius = 0.4;
  int m_clustering = 0; // Inclusive clustering by default
  double m_minPt = 10; // Only relevant for inclusive clustering
  int m_njets = 0; // Only relevant for exclusive clustering

  fastjet::ClusterSequence* m_clustSeq;

};

}
#endif /* RECCALORIMETER_CLUSTERJET_H */
