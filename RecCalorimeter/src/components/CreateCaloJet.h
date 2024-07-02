#ifndef RECCALORIMETER_CREATECALOJET_H
#define RECCALORIMETER_CREATECALOJET_H


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



//// Datamodel
//namespace edm4hep {
//  class ClusterCollection;
//  class ReconstructedParticleCollection;
//}


/** @class CreateCaloJet k4RecCalorimeter/RecFCCeeCalorimeter/src/components/CreateCaloJet.h
 *
 *  @author Jennifer Roloff
 */



class CreateCaloJet : public Gaudi::Functional::Transformer <edm4hep::ReconstructedParticleCollection(const edm4hep::ClusterCollection&), BaseClass_t>  {
public:
  CreateCaloJet(const std::string& name, ISvcLocator* svcLoc);
  edm4hep::ReconstructedParticleCollection operator()(const edm4hep::ClusterCollection& input) const override;

  StatusCode initialize() override;

private:


  double m_minPt = 10;
  double m_jetRadius = 0.4;
  std::string m_jetAlg = "antikt";
  int m_isExclusive = 0;

  k4::recCalo::ClusterJet* clusterer;


};
#endif /* RECCALORIMETER_CREATECALOJET_H */
