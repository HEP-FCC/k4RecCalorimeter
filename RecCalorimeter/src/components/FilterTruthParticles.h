#ifndef RECCALORIMETER_FILTERTRUTHPARTICLES_H
#define RECCALORIMETER_FILTERTRUTHPARTICLES_H


// Gaudi
#include "GaudiAlg/Transformer.h"
#include "k4FWCore/BaseClass.h"


#include "edm4hep/MCParticleCollection.h"

#include "fastjet/JetDefinition.hh"


/** @class FilterTruthParticles k4RecCalorimeter/RecCalorimeter/src/components/FilterTruthParticles.h
 *
 *  @author Jennifer Roloff
 */


class FilterTruthParticles : public Gaudi::Functional::Transformer <edm4hep::MCParticleCollection(const edm4hep::MCParticleCollection&), BaseClass_t>  {
public:
  FilterTruthParticles(const std::string& name, ISvcLocator* svcLoc);
  edm4hep::MCParticleCollection operator()(const edm4hep::MCParticleCollection& input) const override;


private:

  double m_minPt = 0.0;
  int m_useMuons = 0;


};
#endif /* RECCALORIMETER_FILTERTRUTHPARTICLES_H */
