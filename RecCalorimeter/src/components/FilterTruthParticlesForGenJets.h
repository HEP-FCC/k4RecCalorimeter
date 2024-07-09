#ifndef RECCALORIMETER_FILTERTRUTHPARTICLESFORGENJETS_H
#define RECCALORIMETER_FILTERTRUTHPARTICLESFORGENJETS_H


// Gaudi
#include "GaudiAlg/Transformer.h"
#include "k4FWCore/BaseClass.h"


#include "edm4hep/MCParticleCollection.h"

#include "fastjet/JetDefinition.hh"


/** @class FilterTruthParticlesForGenJets k4RecCalorimeter/RecCalorimeter/src/components/FilterTruthParticlesForGenJets.h
 *
 *  Tool to filter truth particles based on their status code. 
 *  Only detector stable particles are stored in a new MCParticleCollection.
 *  This can provide inputs to truth jet clustering.
 *
 *  MinPt: the minimum pT of particles that are stored
 *  SaveMuons: 1 if muons should be included in the collection, 0 otherwise.
 *
 *  @author Jennifer Roloff
 *  @date 2024-7
 */


class FilterTruthParticlesForGenJets : public Gaudi::Functional::Transformer <edm4hep::MCParticleCollection(const edm4hep::MCParticleCollection&), BaseClass_t>  {
public:
  FilterTruthParticlesForGenJets(const std::string& name, ISvcLocator* svcLoc);
  edm4hep::MCParticleCollection operator()(const edm4hep::MCParticleCollection& input) const override;


private:

  double m_minPt = 0.0;
  int m_useMuons = 0;


};
#endif /* RECCALORIMETER_FILTERTRUTHPARTICLES_H */
