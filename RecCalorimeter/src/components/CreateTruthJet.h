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
#include "edm4hep/MCRecoParticleAssociation.h"
#include "edm4hep/MCRecoParticleAssociationCollection.h"
#include <tuple>

#include "ClusterJet.h"


/** @class CreateTruthJet k4RecCalorimeter/RecCalorimeter/src/components/CreateTruthJet.h
 *
 *  Tool for reconstructing jets from truth particles.
 *  It takes as input an MCParticleCollection, and it outputs a ReconstructedParticleCollection with the jets and 
 *  a MCRecoParticleAssociationCollection with the links between the jets and the truth particles.
 *
 *  JetAlg: A string corresponding to the jet algorithm for clustering
 *  JetRadius: The radius of the jets being clustered
 *  MinPt: The pT threshold below which jets are ignored
 *  isExclusiveClustering: 1 if jets should use an exclusive clustering, 0 otherwise
 *  outputAssociation: Name of the output association collection to link jets to their constituents
 *

 *  @author Jennifer Roloff
 *  @date 2024-7
 */


struct CreateTruthJet final : Gaudi::Functional::MultiTransformer<std::tuple< edm4hep::ReconstructedParticleCollection, edm4hep::MCRecoParticleAssociationCollection> (const edm4hep::MCParticleCollection&), BaseClass_t>  {

public:
  CreateTruthJet(const std::string& name, ISvcLocator* svcLoc);
  std::tuple< edm4hep::ReconstructedParticleCollection, edm4hep::MCRecoParticleAssociationCollection> operator()(const edm4hep::MCParticleCollection& input) const override;

  StatusCode initialize() override final;

private:

  double m_minPt = 1;
  double m_jetRadius = 0.4;
  std::string m_jetAlg = "antikt";
  int m_isExclusive = 0;

  k4::recCalo::ClusterJet* clusterer;


};
#endif /* RECCALORIMETER_CREATETRUTHJET_H */
