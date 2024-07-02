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
 *  Tool for reconstructing jets from calorimeter clusters.
 *  It takes as input a ClusterCollection, and it outputs a ReconstructedParticleCollection with the jets.
 *
 *  JetAlg: A string corresponding to the jet algorithm for clustering
 *  JetRadius: The radius of the jets being clustered
 *  MinPt: The pT threshold below which jets are ignored
 *  isExclusiveClustering: 1 if jets should use an exclusive clustering, 0 otherwise
 *

 *  @author Jennifer Roloff
 *  @date 2024-7
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
