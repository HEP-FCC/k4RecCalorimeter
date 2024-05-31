#ifndef RECCALORIMETER_CREATECALOJET_H
#define RECCALORIMETER_CREATECALOJET_H

// std
#include <map>

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiAlg/Transformer.h"
#include "k4FWCore/BaseClass.h"

// k4FWCore
#include "k4FWCore/DataHandle.h"

//EDM4Hep
#include "edm4hep/ReconstructedParticleCollection.h"

// EDM4hep
#include "edm4hep/ClusterCollection.h"


#include "ClusterJet.h"



// Datamodel
namespace edm4hep {
  class ClusterCollection;
  class ReconstructedParticleCollection;
}


/** @class CreateCaloJet k4RecCalorimeter/RecFCCeeCalorimeter/src/components/CreateCaloJet.h
 *
 *  @author Jennifer Roloff
 */


using colltype_in  = edm4hep::ClusterCollection;
using colltype_out = edm4hep::ReconstructedParticleCollection;

class CreateCaloJet : public Gaudi::Functional::Transformer <colltype_out(const colltype_in&), BaseClass_t>  {
public:
  CreateCaloJet(const std::string& name, ISvcLocator* svcLoc);
  colltype_out operator()(const colltype_in& input) const override;

  StatusCode initialize() override;

private:


  double m_minPt = 10;
  double m_jetRadius = 0.4;
  std::string m_jetAlg = "antikt";

  // Map between jet algorithm names and the actual jet clustering algorithm
  const std::map<const std::string, const fastjet::JetAlgorithm> m_jetAlgMap = {
                        {"kt",                    fastjet::JetAlgorithm::kt_algorithm},
                        {"cambridge",             fastjet::JetAlgorithm::cambridge_algorithm},
                        {"antikt",                fastjet::JetAlgorithm::antikt_algorithm},
                        {"genkt",                 fastjet::JetAlgorithm::genkt_algorithm},
                        {"ee_kt",                 fastjet::JetAlgorithm::ee_kt_algorithm},
                        {"ee_genkt",              fastjet::JetAlgorithm::ee_genkt_algorithm},
  };

  ClusterJet* clusterer;


};
#endif /* RECCALORIMETER_CREATECALOJET_H */
