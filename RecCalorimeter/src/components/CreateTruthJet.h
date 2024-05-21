#ifndef RECCALORIMETER_CREATETRUTHJET_H
#define RECCALORIMETER_CREATETRUTHJET_H

// std
#include <cstdint>
#include <map>

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiAlg/Transformer.h"
#include "k4FWCore/BaseClass.h"

// k4FWCore
#include "k4FWCore/DataHandle.h"

#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/MCParticleCollection.h"

#include "fastjet/JetDefinition.hh"


// Datamodel
namespace edm4hep {
class MCParticleCollection;
class ReconstructedParticleCollection;
}




class ClusterInfo : public fastjet::PseudoJet::UserInfoBase{
 public:
 ClusterInfo(const int & index) : _index(index){}
  int index() const{ return _index; }
 protected:
  int _index;
};



/** @class CreateTruthJet k4RecCalorimeter/RecFCCeeCalorimeter/src/components/CreateTruthJet.h
 *
 *  @author Jennifer Roloff
 */

using colltype_in  = edm4hep::MCParticleCollection;
using colltype_out = edm4hep::ReconstructedParticleCollection;

class CreateTruthJet : public Gaudi::Functional::Transformer <colltype_out(const colltype_in&), BaseClass_t>  {
public:
  CreateTruthJet(const std::string& name, ISvcLocator* svcLoc);
  colltype_out operator()(const colltype_in& input) const override;

  StatusCode initialize() override;

private:
  //mutable DataHandle<edm4hep::MCParticleCollection> m_inputParticles {
  //  "MCParticles", Gaudi::DataHandle::Reader, this
  //};

  std::map<std::string, fastjet::JetAlgorithm> m_jetAlgMap = {
                        {"kt",                    fastjet::JetAlgorithm::kt_algorithm},
                        {"cambridge",             fastjet::JetAlgorithm::cambridge_algorithm},
                        {"antikt",                fastjet::JetAlgorithm::antikt_algorithm},
                        {"genkt",                 fastjet::JetAlgorithm::genkt_algorithm},
                        {"ee_kt",                 fastjet::JetAlgorithm::ee_kt_algorithm},
                        {"ee_genkt",              fastjet::JetAlgorithm::ee_genkt_algorithm},
  };

  //DataHandle<edm4hep::ReconstructedParticleCollection> m_jetCollection{"truthJets", Gaudi::DataHandle::Writer, this};

  double m_minPt = 10;
  double m_jetRadius = 0.4;
  std::string m_jetAlg = "antikt";


};
#endif /* RECCALORIMETER_CREATETRUTHJET_H */
