#ifndef RECCALORIMETER_CreateHitTruthLinks_H
#define RECCALORIMETER_CreateHitTruthLinks_H

// k4FWCore
#include "k4FWCore/DataHandle.h"

// Gaudi
#include "Gaudi/Algorithm.h"

// edm4hep
#include "edm4hep/CaloHitSimCaloHitLinkCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CaloHitMCParticleLinkCollection.h"
#include "edm4hep/MCParticleCollection.h"


/** @class CreateHitTruthLinks
 *
 *  Algorithm for creating links between calorimeter cells (CalorimeterHit) and MC particles (MCParticle)
 *  Based on https://github.com/iLCSoft/MarlinReco/blob/master/Analysis/RecoMCTruthLink/src/RecoMCTruthLinker.cc
 *
 *  @author Giovanni Marchiori
 *  @date   2025-06
 *
 */

class CreateHitTruthLinks : public Gaudi::Algorithm {

public:
  CreateHitTruthLinks(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize();

  StatusCode execute(const EventContext&) const;

  StatusCode finalize();

  virtual ~CreateHitTruthLinks();

private:
  /// List of input sim calo hit collections
  Gaudi::Property<std::vector<std::string>> m_hitCollections{
    this, "hits", {}, "Names of SimCalorimeterHit collections to read"};
  /// the vector of input k4FWCore::DataHandles for the input hit collections
  std::vector<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>*> m_hitCollectionHandles;

  /// List of input cell<->hits collections
  Gaudi::Property<std::vector<std::string>> m_cell_hit_linkCollections{
    this, "cell_hit_links", {}, "Names of CaloHitSimCaloHitLink collections to read"};
  /// the vector of input k4FWCore::DataHandles for the input hit<->cell link collections
  std::vector<k4FWCore::DataHandle<edm4hep::CaloHitSimCaloHitLinkCollection>*> m_cell_hit_linkCollectionHandles;

  /// Handle for mc particles (input collection)
  mutable k4FWCore::DataHandle<edm4hep::MCParticleCollection> m_mcparticles{"mcparticles", Gaudi::DataHandle::Reader, this};

  /// Handle for hit<->particle link (output collection)
  mutable k4FWCore::DataHandle<edm4hep::CaloHitMCParticleLinkCollection> m_cell_mcparticle_links{"cell_mcparticle_links", Gaudi::DataHandle::Writer,
    this};
};

#endif /* RECCALORIMETER_CreateHitTruthLinks_H */
