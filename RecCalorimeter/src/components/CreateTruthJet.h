#ifndef RECCALORIMETER_CREATETRUTHJET_H
#define RECCALORIMETER_CREATETRUTHJET_H

// Gaudi
#include "GaudiKernel/ToolHandle.h"
#include "GaudiAlg/Transformer.h"
#include "k4FWCore/BaseClass.h"

// k4FWCore
#include "k4FWCore/DataHandle.h"

#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/MCParticleCollection.h"

#include "ClusterJet.h"



// Datamodel
//namespace edm4hep {
//class MCParticleCollection;
//class ReconstructedParticleCollection;
//}



/** @class CreateTruthJet k4RecCalorimeter/RecCalorimeter/src/components/CreateTruthJet.h
 *
 *  @author Jennifer Roloff
 */


class CreateTruthJet : public Gaudi::Functional::Transformer <edm4hep::ReconstructedParticleCollection(const edm4hep::MCParticleCollection&), BaseClass_t>  {
public:
  CreateTruthJet(const std::string& name, ISvcLocator* svcLoc);
  edm4hep::ReconstructedParticleCollection operator()(const edm4hep::MCParticleCollection& input) const override final;

  StatusCode initialize() override final;

private:


  double m_minPt = 10;
  double m_jetRadius = 0.4;
  std::string m_jetAlg = "antikt";
  int m_isExclusive = 0;

  k4::recCalo::ClusterJet* clusterer;


};
#endif /* RECCALORIMETER_CREATETRUTHJET_H */
