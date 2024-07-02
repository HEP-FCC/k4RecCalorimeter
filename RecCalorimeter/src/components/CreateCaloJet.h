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



/** @class CreateCaloJet k4RecCalorimeter/RecFCCeeCalorimeter/src/components/CreateCaloJet.h
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
