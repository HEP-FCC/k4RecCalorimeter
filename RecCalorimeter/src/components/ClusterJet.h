#ifndef RECCALORIMETER_CLUSTERJET_H
#define RECCALORIMETER_CLUSTERJET_H

// std
#include <map>

/*
// Gaudi
#include "GaudiKernel/ToolHandle.h"
#include "GaudiAlg/Transformer.h"
#include "k4FWCore/BaseClass.h"

// k4FWCore
#include "k4FWCore/DataHandle.h"

#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/MCParticleCollection.h"
*/

#include "fastjet/JetDefinition.hh"


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
  ClusterJet(std::string jetAlg, double jetRadius, double minPt);

  //StatusCode initialize();

  std::vector<fastjet::PseudoJet> cluster(const   std::vector<fastjet::PseudoJet> clustersPJ) const;

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
  double m_minPt = 10;
  double m_jetRadius = 0.4;


};

#endif /* RECCALORIMETER_CLUSTERJET_H */
