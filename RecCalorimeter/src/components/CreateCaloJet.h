#ifndef RECCALORIMETER_CREATECALOJET_H
#define RECCALORIMETER_CREATECALOJET_H

// std
//#include <cstdint>
#include <map>

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"

// k4FWCore
#include "k4FWCore/DataHandle.h"

//EDM4Hep
#include "edm4hep/ReconstructedParticleCollection.h"


//Fastjet
#include "fastjet/JetDefinition.hh"


// Datamodel
namespace edm4hep {
  class ClusterCollection;
  class ReconstructedParticleCollection;
}



// Attach information to pseudojets to be able to map the cluster
// index to the pseudojet
class ClusterInfo : public fastjet::PseudoJet::UserInfoBase{
 public:
 ClusterInfo(const int & index) : _index(index){}
  int index() const{ return _index; }
 protected:
  int _index;
};



/** @class CreateCaloJet k4RecCalorimeter/RecFCCeeCalorimeter/src/components/CreateCaloJet.h
 *
 *  @author Jennifer Roloff
 */

class CreateCaloJet : public GaudiAlgorithm {
public:
  CreateCaloJet(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:

  DataHandle<edm4hep::ClusterCollection> m_inputClusters {"CorrectedCaloClusters", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::ReconstructedParticleCollection> m_jetCollection{"jets", Gaudi::DataHandle::Writer, this};

  double m_minPt = 10;
  double m_jetRadius = 0.4;
  std::string m_jetAlg = "antikt";

  // Map between jet algorithm names and the actual jet clustering algorithm
  std::map<std::string, fastjet::JetAlgorithm> m_jetAlgMap = {
                        {"kt",                    fastjet::JetAlgorithm::kt_algorithm},
                        {"cambridge",             fastjet::JetAlgorithm::cambridge_algorithm},
                        {"antikt",                fastjet::JetAlgorithm::antikt_algorithm},
                        {"genkt",                 fastjet::JetAlgorithm::genkt_algorithm},
                        {"ee_kt",                 fastjet::JetAlgorithm::ee_kt_algorithm},
                        {"ee_genkt",              fastjet::JetAlgorithm::ee_genkt_algorithm},
  };


};
#endif /* RECCALORIMETER_CREATECALOJET_H */
