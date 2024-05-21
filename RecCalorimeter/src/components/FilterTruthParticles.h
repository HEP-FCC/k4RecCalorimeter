#ifndef RECCALORIMETER_FILTERTRUTHPARTICLES_H
#define RECCALORIMETER_FILTERTRUTHPARTICLES_H


// Gaudi
#include "GaudiAlg/Transformer.h"
#include "k4FWCore/BaseClass.h"


#include "edm4hep/MCParticleCollection.h"

#include "fastjet/JetDefinition.hh"



// Datamodel
namespace edm4hep {
class MCParticleCollection;
}


/** @class FilterTruthParticles k4RecCalorimeter/RecCalorimeter/src/components/FilterTruthParticles.h
 *
 *  @author Jennifer Roloff
 */

using colltype_in  = edm4hep::MCParticleCollection;
using colltype_out = edm4hep::MCParticleCollection;

class FilterTruthParticles : public Gaudi::Functional::Transformer <colltype_out(const colltype_in&), BaseClass_t>  {
public:
  FilterTruthParticles(const std::string& name, ISvcLocator* svcLoc);
  colltype_out operator()(const colltype_in& input) const override;


private:

  double m_minPt = 0.0;
  int m_useMuons = 0;


};
#endif /* RECCALORIMETER_FILTERTRUTHPARTICLES_H */
