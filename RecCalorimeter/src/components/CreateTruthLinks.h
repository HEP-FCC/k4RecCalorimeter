#ifndef RECCALORIMETER_CreateTruthLinks_H
#define RECCALORIMETER_CreateTruthLinks_H

// k4FWCore
#include "k4FWCore/DataHandle.h"

// Gaudi
#include "Gaudi/Algorithm.h"

// edm4hep
#include "edm4hep/CaloHitSimCaloHitLinkCollection.h"
#include "edm4hep/CaloHitMCParticleLinkCollection.h"
#include "edm4hep/ClusterMCParticleLinkCollection.h"
#include "edm4hep/MCParticleCollection.h"


/** @class CreateTruthLinks
 *
 *  Algorithm for creating links between reco and truth:
 *  - calorimeter cells (CalorimeterHit) -> MC particles (MCParticle)
 *  - calorimeter clusters -> MC particles (weight = fraction of energy of given cluster due to a particular MC particle)
 *  - MC particles -> calorimeter clusters (DOES NOT SEEM POSSIBLE IN EDM4HEP) => same relations as previous map but different weight: w = fraction of total calo hit energy from given particle that goes in particular cluster
 *  Based on https://github.com/iLCSoft/MarlinReco/blob/master/Analysis/RecoMCTruthLink/src/RecoMCTruthLinker.cc
 *
 *  @author Giovanni Marchiori
 *  @date   2025-07
 *
 */

class CreateTruthLinks : public Gaudi::Algorithm {

public:
  CreateTruthLinks(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize();

  StatusCode execute(const EventContext&) const;

  StatusCode finalize();

  virtual ~CreateTruthLinks();

private:
  /// List of input cell<->hits collections
  Gaudi::Property<std::vector<std::string>> m_cell_hit_linkCollections{
    this, "cell_hit_links", {}, "Names of CaloHitSimCaloHitLink collections to read"};
  /// the vector of input k4FWCore::DataHandles for the input hit<->cell link collections
  std::vector<k4FWCore::DataHandle<edm4hep::CaloHitSimCaloHitLinkCollection>*> m_cell_hit_linkCollectionHandles;
  /// List of input cluster collections
  Gaudi::Property<std::vector<std::string>> m_clusterCollections{
    this, "clusters", {}, "Names of Cluster collections to read"};
  /// the vector of input k4FWCore::DataHandles for the input cluster collections
  std::vector<k4FWCore::DataHandle<edm4hep::ClusterCollection>*> m_clusterCollectionHandles;

  /// Handle for mc particles (input collection)
  mutable k4FWCore::DataHandle<edm4hep::MCParticleCollection> m_mcparticles{"mcparticles", Gaudi::DataHandle::Reader, this};

  /// Handle for hit<->particle link (output collection)
  mutable k4FWCore::DataHandle<edm4hep::CaloHitMCParticleLinkCollection> m_cell_mcparticle_links{"cell_mcparticle_links", Gaudi::DataHandle::Writer,
    this};

  /// Handle for cluster<->particle link (output collection)
  mutable k4FWCore::DataHandle<edm4hep::ClusterMCParticleLinkCollection> m_cluster_mcparticle_links{"cluster_mcparticle_links", Gaudi::DataHandle::Writer,
    this};

  /// Handle for particle<->cluster link (output collection) [NOT POSSIBLE]
  //  mutable k4FWCore::DataHandle<edm4hep::MCParticleClusterLinkCollection> m_mcparticle_cluster_links{"mcparticle_cluster_links", Gaudi::DataHandle::Writer,
  //    this};
    
};

#endif /* RECCALORIMETER_CreateTruthLinks_H */
