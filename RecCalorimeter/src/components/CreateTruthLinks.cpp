#include "CreateTruthLinks.h"

DECLARE_COMPONENT(CreateTruthLinks)

CreateTruthLinks::CreateTruthLinks(const std::string& name, ISvcLocator* svcLoc)
    : Gaudi::Algorithm(name, svcLoc) {
  declareProperty("mcparticles", m_mcparticles, "MC particles collection (input)");
  declareProperty("cell_mcparticle_links", m_cell_mcparticle_links, "The links between cells and MC particles (output)");
  declareProperty("cluster_mcparticle_links", m_cluster_mcparticle_links, "The links between clusters and MC particles (output)");
}

CreateTruthLinks::~CreateTruthLinks() { }

StatusCode CreateTruthLinks::initialize() {
  StatusCode sc = Gaudi::Algorithm::initialize();
  if (sc.isFailure())
    return sc;

  // create handles for input collections
  for (const auto& col : m_cell_hit_linkCollections) {
    debug() << "Creating handle for input cell-hit link collection : " << col << endmsg;
    try {
      m_cell_hit_linkCollectionHandles.push_back(
          new k4FWCore::DataHandle<edm4hep::CaloHitSimCaloHitLinkCollection>(col, Gaudi::DataHandle::Reader, this));
    } catch (...) {
      error() << "Error creating handle for input collection: " << col << endmsg;
      return StatusCode::FAILURE;
    }
  }

  for (const auto& col : m_clusterCollections) {
    debug() << "Creating handle for input cluster collection : " << col << endmsg;
    try {
      m_clusterCollectionHandles.push_back(
          new k4FWCore::DataHandle<edm4hep::ClusterCollection>(col, Gaudi::DataHandle::Reader, this));
    } catch (...) {
      error() << "Error creating handle for input collection: " << col << endmsg;
      return StatusCode::FAILURE;
    }
  }


  info() << "CreateTruthLinks initialized" << endmsg;

  return StatusCode::SUCCESS;
}

StatusCode CreateTruthLinks::execute(const EventContext&) const {

  // map from true tracks linked to hit contributions to
  // those that really should have been linked
  std::map< int , int > remap_as_you_go ;
  bool doRemapping = true;  // can turn off for debug

  // This step is needed for clusters, because the first few branchings in the
  // shower are kept in the MCParticle collection.
  // We therefore need to back-track to the particle actually entering the calorimeter.
  // The originator of a hit sometimes IS a generator particle, in which case
  // there is no problem.
  // Otherwise, the criteria to keep a particle as a hit originator is:
  //   The particle is an ancestor of the one initial linked to the hit.
  //      AND
  //   [
  //     The particle is a generator particle
  //         OR
  //     The particle starts in the tracker, ends in calorimeter. That it starts in the tracker
  //     is determined by the fact that it's mother ends there. However, note the possibility of
  //     "non-destructive interactions", where a particles mother does *not* end at the point where
  //     the particle is created !
  //         OR
  //     The particle, or one of it's ancestors is a back-scatterer. This case is further treated in the
  //     loop over clusters, to catch the rather common case that a back-scattered particle re-enters the
  //     calorimeter close to it's start-start point (think 15 MeV charged particle!) and the hits it
  //     produces is in the same cluster as that of the initiator of the shower the backscatter come from.
  //     (This can't be detected in this first loop, since we don't know about clusters here)
  //   ]
  // A map is set up, so obviously the first check is wether the MCP of the hit has already been
  // treated when looking at some previous hit, in which case one just thake that association.

  // create cell<->particle links
  edm4hep::CaloHitMCParticleLinkCollection* caloHitMCParticleLinkCollection = m_cell_mcparticle_links.createAndPut();

  // create cluster<->particle links
  edm4hep::ClusterMCParticleLinkCollection* clusterMCParticleLinkCollection = nullptr;
  if (m_clusterCollections.size()>0) {
    clusterMCParticleLinkCollection = m_cluster_mcparticle_links.createAndPut();
  }

  // retrieve MC particles
  const edm4hep::MCParticleCollection* mcparticles = m_mcparticles.get();

  debug() << "Creating CaloHit <-> MCParticle links : " << endmsg;

  // loop over calo<->sim hit collection
  std::map<int, int> nhits; // map of nhits with given ID in link collections
  for (size_t ih = 0; ih < m_cell_hit_linkCollectionHandles.size(); ih++) {
    debug() << "Processing input sim <-> calo hit link collection " << ih << " : "
            << m_cell_hit_linkCollectionHandles[ih]->objKey() << endmsg;
    const edm4hep::CaloHitSimCaloHitLinkCollection* caloHitSimCaloHitLinks = m_cell_hit_linkCollectionHandles[ih]->get();
    debug() << "Collection size: " << caloHitSimCaloHitLinks->size() << endmsg;

    // Loop over the G4 hits, find the associated calo hits, then loop over
    // the G4 hit contribution to calculate contribution of each MCParticle to the calo hit
    for (const auto &assoc : *caloHitSimCaloHitLinks){
      debug() << "Processing new sim <-> calo hit link: " << endmsg;
      const auto &simHit = assoc.getTo();
      const auto &caloHit = assoc.getFrom();
      debug() << "Sim hit id and index: " << simHit.id() << " " << simHit.id().index << endmsg;
      debug() << "Calo hit id and index: " << caloHit.id() << " " << caloHit.id().index << endmsg;
      if (nhits[simHit.id().index + ih*(1<<24)]>0) {
        warning() << "Sim hit with more than one calo hit!!!" << endmsg;  // GM: maybe this will be possible if we split hits based on timing?
        continue;
      }
      else {
        nhits[simHit.id().index + ih*(1<<24)]=1;
      }

      // extract digitised over deposited energy
      // assuming here there is at most one simhit per calo hit (0 in the case of noise hit)
      double calib_factor = caloHit.getEnergy()/simHit.getEnergy();
      debug() << "Calib factor = " << calib_factor << endmsg;

      // now loop over truth contributions
      int k=-1;
      std::map< size_t , double > simHitMapEnergy ;  //  counts total energy for every MCParticle for given sim hit
      debug() << "Looping over contributions" << endmsg;
      for (auto &simHitContrib: simHit.getContributions())
      {
        k++;
        // mcp: original MCParticle associated to contribution
        const edm4hep::MCParticle& origmcp = simHitContrib.getParticle();
        int index_mcp = origmcp.id().index;
        double e  = simHitContrib.getEnergy() * calib_factor;
        debug() << "initial true contributor "<< k << " id " << origmcp.id()
                <<" (with E = " << origmcp.getEnergy() << " and pdg " << origmcp.getPDG() << " ) e hit: " << e << endmsg;

        if (doRemapping) {
          if ( remap_as_you_go.find(index_mcp) != remap_as_you_go.end() ) {
            // first condition: we have already found from some earlier hit how to remap this particle
            index_mcp = remap_as_you_go.find(index_mcp)->second;
            {
              const edm4hep::MCParticle& amcp = mcparticles->at(index_mcp);
              debug() << "Hit " << simHit.id() << " / " << caloHit.id()
                      << " attributed to " << amcp.id()
                      << " because its origin mcp has already been treated : originator case 0 " << endmsg;
            }
          } else {
            // otherwise: we have not yet decided if this particle has to be remapped or not
            int index_mother = -1;
            int index_this_Kid = index_mcp;

            if ( origmcp.getGeneratorStatus() == 0 && ( origmcp.getParents().empty() == false ) ) {

              // second condition: particles not from the generator that have parents
              // find which true particle this hit should really be attributed to, by
              // tracking back the history.

              // Two cases to treat:
              //   For some reason, the hit is not attributed to the
              //   incoming particle, but to some particle in the shower.
              //   Just track back to the incoming particle, which might (usually)
              //   be a gen-stat 1 particle from the main vertex. It can also
              //   be from a decay or interaction in the tracking made by
              //   Geant : genstat=0, mother isDecayedInTracker, (or
              //   the production vertex is in the track, but the mother
              //   continues on - ticky case, see comment futher down), or a
              //   decyed particle from the generator (genstat 2) that
              //   hits the calo before decaying (lambdas, K^0_S)

              //   Or: for back-scatterers the calorimeter hit is attributed to the
              //   the last particle *even if this particle both
              //   started and ended in the tracker* !!!! Then we back-track
              //   until we find a particle which at least started inside
              //   the calorimeter, and then go on backtracking as above.
              //   This case is triggered by the particle linked to the
              //   hit being DecayedInTracker, hence the case where a
              //   back-scatter actually ends in the calorimeter is
              //   treated as the first case.

              debug() << "simHit " << simHit.id() << " not created by generator particle. backtracking ..." << endmsg;
              debug() << "starting from particle " << origmcp.id() << " gs "<<origmcp.getGeneratorStatus()
                      << " dint "<< origmcp.isDecayedInTracker()
                      << " bs "<< origmcp.isBackscatter()
                      << " ndi "<< origmcp.vertexIsNotEndpointOfParent()
                      << " npar "<< origmcp.getParents().size()
                      << " pdg "<< origmcp.getPDG()
                      << " vtx"
                      << " " << origmcp.getVertex()[0]
                      << " " << origmcp.getVertex()[1]
                      << " " << origmcp.getVertex()[2]
                      << " end"
                      << " " << origmcp.getEndpoint()[0]
                      << " " << origmcp.getEndpoint()[1]
                      << " " << origmcp.getEndpoint()[2] << endmsg;

              index_mother = origmcp.getParents()[0].id().index;
              if (
                index_mother >= 0 &&
                !mcparticles->at(index_this_Kid).isBackscatter() &&
                mcparticles->at(index_mother).getParents().size()>0 &&
                mcparticles->at(index_mother).getGeneratorStatus()==0 &&
                !mcparticles->at(index_mother).isDecayedInTracker()
              ) {
                debug() << "going into into originator loop " << endmsg;
              }
              while (
                index_mother>=0 &&
                !mcparticles->at(index_this_Kid).isBackscatter() &&
                mcparticles->at(index_mother).getParents().size()>0 &&
                mcparticles->at(index_mother).getGeneratorStatus()==0 &&
                !mcparticles->at(index_mother).isDecayedInTracker()
              ) {
                // back-track as long as there is a non-generator
                // mother, or the mother decayed in the tracker
                // (=> the kid is the particle entering the calorimeter.)
                debug() << "in originator loop " << endmsg;
                const edm4hep::MCParticle& mother = mcparticles->at(index_mother);
                debug() << "shower-part mother "<<mother.id()
                        << " gs "<<mother.getGeneratorStatus()
                        << " dint "<<mother.isDecayedInTracker()
                        << " bs "<<mother.isBackscatter()
                        << " ndi "<<mother.vertexIsNotEndpointOfParent()
                        << " npar "<<mother.getParents().size()
                        << " pdg "<<mother.getPDG() << endmsg;
                // case shower-particle (???)
                // if ( this_Kid.vertexIsNotEndpointOfParent() != _invertedNonDestructiveInteractionLogic) {
                if ( mcparticles->at(index_this_Kid).vertexIsNotEndpointOfParent() ) {
                  const edm4hep::MCParticle& gmother = mcparticles->at(index_mother).getParents()[0];
                  if ( gmother.isDecayedInTracker() ) {
                    debug() << "break out : grandmother "<< gmother.id()
                            << " gs "<< gmother.getGeneratorStatus()
                            << " dint "<< gmother.isDecayedInTracker()
                            << " bs "<< gmother.isBackscatter()
                            << " ndi "<< gmother.vertexIsNotEndpointOfParent()
                            << " npar "<< gmother.getParents().size()
                            << " pdg "<< gmother.getPDG() << endmsg;
                    break ;
                  }
                }
                index_this_Kid = index_mother ;
                // mother= dynamic_cast<edm4hep::MCParticle*>(mother.getParents()[0]); // (assume only one...)
                if (mcparticles->at(index_mother).getParents().size()>0)
                  index_mother = mcparticles->at(index_mother).getParents()[0].id().index; // (assume only one...)
                else
                  index_mother = -1;
              } // end while backtracking loop

              // Further treatment (basically determining if it this_Kid or mother that entered the
              // calorimeter) based on why we left the while-loop

              // here one of the while conditions is false, ie. at least one of
              // "kid is back-scatter", "no mother" , "mother has no parents", "mother is from
              // generator", or "mother did decay in tracker" must be true. We know that this_Kid
              // fulfills all the while-conditions except the first: obvious if at least one iteration of
              // the loop was done, since it was the mother in the previous iteration,
              // but also true even if the loop wasn't transversed, due to the conditions
              // to at all enter this block of code (explicitly must be simulator particle with
              // mother, implicitly must have ended in the calo, since it did make calo hits.)
              if (mcparticles->at(index_this_Kid).isBackscatter() ) {
                // case 2: Kid is back-scatterer. It has thus started in a calo, and entered
                // from there into the tracking volume, and did cause hits after leaving the tracker
                // volume again ->  this_Kid started before the calo, and is the one
                remap_as_you_go[index_mcp] = index_this_Kid;
                {
                  const edm4hep::MCParticle& this_Kid = mcparticles->at(index_this_Kid);
                  debug() << "Hit " << simHit.id() << " / " << caloHit.id()
                          << " attributed to kid " << this_Kid.id()
                          << " because it's origin is a back-scatter : originator case 2 " << endmsg;
                  debug() << "Kid: " << this_Kid.id()
                          << " gs " << this_Kid.getGeneratorStatus()
                          << " dint " << this_Kid.isDecayedInTracker()
                          << " bs " << this_Kid.isBackscatter()
                          << " npar " << this_Kid.getParents().size()
                          << " pdg " << this_Kid.getPDG() << endmsg;
                }
              } else if (
                index_mother >= 0 &&
                mcparticles->at(index_mother).isDecayedInTracker()
              ) {
                // the clear-cut case:
                // this_Kid started before the calo, and is the one
                // the hit should be attributed to
                remap_as_you_go[index_mcp] = index_this_Kid;
                {
                  const edm4hep::MCParticle& this_Kid = mcparticles->at(index_this_Kid);
                  debug() << "Hit " << simHit.id() << " / " << caloHit.id()
                          << " attributed to kid " << this_Kid.id()
                          << " because it's origin is in tracker : originator case 3 " << endmsg;
                  debug() << "Kid: " << this_Kid.id()
                          << " gs " << this_Kid.getGeneratorStatus()
                          << " dint " << this_Kid.isDecayedInTracker()
                          << " bs " << this_Kid.isBackscatter()
                          << " npar " << this_Kid.getParents().size()
                          << " pdg " << this_Kid.getPDG() << endmsg;
                }
              } else {
                // the other three cases, ie. one or several of "no mother", "no grand-mother",
                // "generator particle" + that we know that "mother decayed in calo"
                if ( index_mother < 0 ) {
                  // ... which of course implies no grand-mother, and no gen stat
                  // of the mother as well -> should not be possible !
                  warning() << "MCparticle " << mcparticles->at(index_this_Kid).id()
                            << " is a simulation particle, created in the calorimeter by nothing . " << endmsg;
                  remap_as_you_go[index_mcp] = index_this_Kid;
                  index_mcp = index_this_Kid; // can't do better than that.
                } else {
                  // here we know: mother exists, but decayed in calo. In addition, two possibilities:
                  // mother is generator particle, or there was no grand-parents. One or both
                  // must be true here. Here it gets complicated, because what we want to know is
                  // whether this_Kid started in the tracker or not. Unluckily, we don't know that
                  // directly, we only know where the mother ended. IF ithe mother ended in the tracker,
                  // there is no problem, and has already been treated, but if it ended in the calo, it is
                  // still possible that this_Kid came from a "non-destructive interaction" with the tracke-detector
                  // material. This we now try to figure out.

                  if ( mcparticles->at(index_this_Kid).vertexIsNotEndpointOfParent() == false ) {
                    // This bizare condition is due to a bug in LCIO (at least for the DBD samples. this_Kid.vertexIsNotEndpointOfParent()
                    // should be true in the "non-destructive interaction", but it isn't: actually it is "false", but is "true" for the
                    // particles that *do* originate at the end-vertex of their parent. This is a bug in Mokka.
                    // Hence the above ensures that there was NO "non-destructive interaction", and it is clear this_Kid was created at the end-point
                    // of the mother. The mother is either a generator particle (to be saved), or the "Eve" of the decay-chain (or both).
                    // So mother is the one to save and assign the hits to:
                    remap_as_you_go[index_mcp] = index_mother;
                    index_mcp = index_mother;   // mother started at ip, and reached the calo, and is the one
                                  // the hit should be attributed to.

                    debug() << "Hit " << simHit.id() << " / " << caloHit.id()
                            << " attributed to mother " << mcparticles->at(index_mcp).id()
                            << " because it is a generator particle or started in tracker : originator case 4 " << endmsg;

                  } else {
                    // here we DO have a "non-destructive interaction". Unluckily, we cant directly know if this took place in the
                    // tracker (in which case we should keep this_Kid as the mcp to save and assign hits to), or not.
                    // We will play a few clean tricks to find the cases where it either certain that the
                    // "non-destructive interation" was in the tracker, or that it was in the calo. This reduces
                    // the number of uncertain cases to play dirty tricks with.

                    // Clean tricks to play: look at the sisters of this_Kid: with some luck one of them is a promptly decaying
                    // particle, eg. a pi0.
                    // This sister will be flagged as decayed in calo/tracker, and from that we know for certain that the
                    // "non-destructive interaction" was in the calo/tracker.
                    // It can also be that one of the sisters is flagged as a back-scatter, which only happens in the calo.
                    // If the particles grand-mother is decayed in calo, and the mother isn't from a "non-destructive interaction",
                    // the "non-destructive interaction" was in the calo.

                    // Finally, we can check the distance of the end-point of the mother (sure to be in the calo) to
                    // the vertex of this_Kid. If this is small, this_Kid *probably* started in the calo.
                    int starts_in_tracker = 0 ;
                    int has_pi0 = 0 ;
                    int gmother_in_calo = 0;
                    double rdist =0.;
                    int has_bs = 0;

                    const edm4hep::MCParticle& mother = mcparticles->at(index_mother);
                    debug() << "Non destructive interaction, looping over sisters" << endmsg;
                    for ( unsigned kkk=0 ; kkk < mother.getDaughters().size() ; kkk++ ) {
                      // edm4hep::MCParticle* sister = dynamic_cast<edm4hep::MCParticle*>(mother.getDaughters()[kkk]);

                      const edm4hep::MCParticle& sister = mother.getDaughters()[kkk];
                      if ( sister.id() == mcparticles->at(index_this_Kid).id() ) continue;
                      if ( abs(sister.getVertex()[0]-mcparticles->at(index_this_Kid).getVertex()[0]) > 0.1 ||
                            abs(sister.getVertex()[1]-mcparticles->at(index_this_Kid).getVertex()[1]) > 0.1 ||
                            abs(sister.getVertex()[2]-mcparticles->at(index_this_Kid).getVertex()[2]) > 0.1 ) continue;
                          // must check that it is the same vertex:
                          // several "non-destructive interactions" can
                          // take place (think delta-rays !)
                      if ( sister.isBackscatter() ) {
                        has_bs = 1 ;
                        break ;
                      } else if ( sister.isDecayedInTracker() ) {
                        starts_in_tracker = 1 ;
                        break ;
                      }
                      // any pi0:s at all ? (it doesn't matter that we break at the two cases above,
                      // because if we do, it doesn't matter if there are
                      // pi0 sisters or not !)
                      if ( sister.getPDG() == 111 ) {
                        has_pi0 = 1 ;
                      }
                    }
                    // if not already clear-cut, calculate distance vertex to mother end-point
                    if ( starts_in_tracker != 1 && has_bs != 1 && has_pi0 != 1 && gmother_in_calo != 1 ) {
                      rdist = sqrt(pow(mother.getEndpoint()[0]-mcparticles->at(index_this_Kid).getVertex()[0],2)+
                                  pow(mother.getEndpoint()[1]-mcparticles->at(index_this_Kid).getVertex()[1],2)+
                                  pow(mother.getEndpoint()[2]-mcparticles->at(index_this_Kid).getVertex()[2],2));
                      if ( mother.getParents().size() != 0 ) {
                        const edm4hep::MCParticle& gmother = mother.getParents()[0];
                        if ( gmother.isDecayedInCalorimeter() ) {
                          gmother_in_calo = 1 ;
                        }
                      }
                    }
                    debug() << "starts_in_tracker, has_pi0, has_bs, gmother_in_calo : " <<  starts_in_tracker << " " << has_pi0 << " " <<  has_bs << " " << gmother_in_calo << endmsg;
                    if ( starts_in_tracker == 1 ) {
                      // this_Kid is a clear-cut hit-originator
                      // this_Kid started before the calo, and is the one
                      // the hit should be attributed to
                      remap_as_you_go[index_mcp] = index_this_Kid;
                      index_mcp = index_this_Kid;
                      const edm4hep::MCParticle& this_Kid = mcparticles->at(index_this_Kid);
                      debug() << "Hit " << simHit.id() << " / " << caloHit.id()
                              << " attributed to kid " << this_Kid.id()
                              << " because it's origin could be deduced to be in tracker : originator case 5 " << endmsg;

                      debug() << "Details of case 5: "
                              << this_Kid.getVertex()[0] << " "
                              << this_Kid.getVertex()[1] << " "
                              << this_Kid.getVertex()[2] << " "
                              << mother.getEndpoint()[0] << " "
                              << mother.getEndpoint()[1] << " "
                              << mother.getEndpoint()[2] << " "
                              << mother.getGeneratorStatus() <<  " "
                              << this_Kid.vertexIsNotEndpointOfParent() << endmsg;
                    } else if ( has_pi0 != 0 || has_bs != 0 || gmother_in_calo != 0 ) {
                      // clear-cut case of this_Kid starting in the calo.
                      // We do know that the mother
                      // is a generator particle and/or the "Eve" of the cascade,
                      // so we should attribute hits to the
                      // mother and save it.
                      // mother started at ip, and reached the calo, and is the one
                      // the hit should be attributed to.
                      remap_as_you_go[index_mcp] = index_mother;
                      index_mcp = index_mother;
                      const edm4hep::MCParticle& this_Kid = mcparticles->at(index_this_Kid);
                      const edm4hep::MCParticle& mcp = mcparticles->at(index_mcp);
                      debug() << "Hit " << simHit.id() << " / " << caloHit.id()
                              << " attributed to mother " << mcp.id()
                              << " because it's origin could be deduced to be in tracker : originator case 6 " << endmsg;
                      debug() << "Case 6 details: kid starts in calo " << " "
                              << this_Kid.getVertex()[1] << " " << this_Kid.getVertex()[2] << " " << rdist <<  " "
                              << has_pi0  << " " << has_bs << " " <<  gmother_in_calo << endmsg;
                    } else {
                      // un-clear case: no pi0 nor back-scatteres among the sisters to help to decide.
                      // Use distance this_Kid-startpoint to mother-endpoint. We know that the latter is
                      // in the calo, so if this is small, guess that the start point of this_Kid is
                      // also in the calo. Calos are dense, so typically in the case the "non-destructive interaction"
                      // is in the calo, one would guess  that the distance is small, ie. large distance ->
                      // unlikely that it was in the calo.

                      if ( rdist > 200. ) {
                        // guess "non-destructive interaction" not in calo -> this_Kid is originator
                        // this_Kid started before the calo, and is the one
                        // the hit should be attributed to
                        remap_as_you_go[index_mcp] = index_this_Kid;
                        index_mcp = index_this_Kid;
                        const edm4hep::MCParticle& this_Kid = mcparticles->at(index_this_Kid);
                        const edm4hep::MCParticle& mcp = mcparticles->at(index_mcp);

                        debug() << "Hit " << simHit.id() << " / " << caloHit.id()
                                << " attributed to kid " << mcp.id()
                                << " because it's origin is guessed to be in tracker : originator case 7 " << endmsg;
                        debug() << "Case 7 details: guess kid starts in tracker " << " "
                                << this_Kid.getVertex()[1] << " "
                                << this_Kid.getVertex()[2] << " "
                                << this_Kid.getVertex()[3] << " "
                                << rdist <<  " "
                                << has_pi0  << " " << endmsg;
                      } else {
                        // guess in calo -> mother is originator
                        // mother started at ip, and reached the calo, and is the one
                        // the hit should be attributed to.
                        remap_as_you_go[index_mcp] = index_mother;
                        index_mcp = index_mother;
                        const edm4hep::MCParticle& this_Kid = mcparticles->at(index_this_Kid);
                        const edm4hep::MCParticle& mcp = mcparticles->at(index_mcp);

                        debug() << "Hit " << simHit.id() << " / " << caloHit.id()
                                << " attributed to mother " << mcp.id()
                                << " because it's origin is guessed in tracker : originator case 8 " <<endmsg;
                        debug() << "Case 8 details: guess kid starts in calo " << " "
                                << this_Kid.getVertex()[1] << " "
                                << this_Kid.getVertex()[2] << " "
                                << this_Kid.getVertex()[3] << " "
                                << rdist <<  " "
                                << has_pi0  << " " << endmsg;
                      }
                    }
                  }
                }
              }
            } else {
              // first case: the hit generating mcp itself already fulfills the first
              // criterium, ie. it is a generator particle (Gen status == 1 or no parents)
              debug() << "Hit " << simHit.id() << " / " << caloHit.id()
                      << " attributed to " << origmcp.id()
                      << " because its origin is a generator particle : originator case 1 " << endmsg;

              remap_as_you_go[index_mcp] = index_mcp;
            } // genstat if - then - else
          }
        } // end if do Remapping

        const edm4hep::MCParticle& mcp = mcparticles->at(index_mcp);
        debug() << "Final assignment for contribution " << k << " to " << simHit.id() << " / "
                << caloHit.id() << " : " << mcp.id() << endmsg;
        debug() << "gs "<< mcp.getGeneratorStatus()
                << " dint " << mcp.isDecayedInTracker()
                << " bs "<< mcp.isBackscatter()
                << " npar "<< mcp.getParents().size()
                << " pdg "<< mcp.getPDG()
                << " end"
                << " " << mcp.getEndpoint()[0]
                << " " << mcp.getEndpoint()[1]
                << " " << mcp.getEndpoint()[2] << endmsg;

        simHitMapEnergy[ mcp.id().index ] += e;
      } // end loop over contributions
      double sumw=0.0; // for debug
      for (const auto &mcp: *mcparticles) {
        // create calo hit<->mc particle associations
        if (simHitMapEnergy[mcp.id().index]>0) {
          auto link = caloHitMCParticleLinkCollection->create();
          link.setFrom(caloHit);
          link.setTo(mcp);
          double w = simHitMapEnergy[mcp.id().index] / caloHit.getEnergy();
          sumw += w;
          link.setWeight(w);
        }
      }
      debug() << "Sum of weights for this calo hit = " << sumw << endmsg;

      if (nhits[simHit.id().index + ih*(1<<24)]==0) {
          warning() << "Sim hit with no associated  calo hit!!!" << endmsg;
          // might happen if digitiser creates cells only if energy is above a given threshold..
      }
    } // end loop over calo hit <-> sim calo hit links
  } // end loop over calo hit <-> sim calo hit collections

  debug() << "Finished linking calo hits to MC particles" << endmsg;

  debug() << "Creating Cluster<->MCParticle links using CaloHit<->MCParticle links, re-assigning the latter in some rare cases" << endmsg;
  // loop over cluster collections
  for (size_t ih = 0; ih < m_clusterCollectionHandles.size(); ih++) {
    debug() << "Processing input cluster collection " << ih << " : "
            << m_clusterCollectionHandles[ih]->objKey() << endmsg;
    const edm4hep::ClusterCollection* clusters = m_clusterCollectionHandles[ih]->get();
    debug() << "Collection size: " << clusters->size() << endmsg;

    // Loop over the clusters, find the associated calo hits, use the linked MCParticles
    // to calculate the energy contributed by a given particle and set links and weights
    for (const edm4hep::Cluster &cluster : *clusters){

      std::map< int , double > mcpEnergy;  // map of (index_MCParticle, total energy contributed to the cluster)

      // We need to find all seen hits this clutser is made of, which sim hits each
      // of the seen hits came from, and finally which true particles actually created
      // each sim hit. Contrary to the sim tracker hits above, a sim-calo hit can be
      // made by several true particles. They also have a signal size (energy) value.
      // In addition, the true particle creating sometimes needs to be back-tracked
      // to the particle actually entering tha calorimeter.

      debug() << endmsg;
      debug() << "=================" << endmsg;
      debug() << "Processing new cluster: " << endmsg;
      debug() << "Cluster id = "<< cluster.id() << " with E = " << cluster.getEnergy() << " , "  << cluster.hits_size() << " hits " << endmsg;

      double ecalohitsum=0.;
      double ecalohitsum_known=0.;
      double ecalohitsum_unknown=0.;

      // loop over the cluster hits
      for (auto it = cluster.hits_begin(); it != cluster.hits_end(); it++) {
        auto cell = *it;
        float eCell = cell.getEnergy();
        ecalohitsum += eCell;

        // loop over hits -> MCParticle links to calculate contribution from given particle
        // GM, note: original code in https://github.com/iLCSoft/MarlinReco/blob/02a01cfe6154fa42b31081250bf84c8f8718f0b1/Analysis/RecoMCTruthLink/src/RecoMCTruthLinker.cc#L1167
        // is quite more complex, and tries to reassign calo->particle links for some rare cases
        for (const auto &assoc : *caloHitMCParticleLinkCollection){
          const auto &caloHit = assoc.getFrom();
          if (caloHit.id() != cell.id()) continue;
          const auto &mcp = assoc.getTo();
          double w = assoc.getWeight(); // fraction of energy of this hit due to mcp
          mcpEnergy[mcp.id().index] += eCell*w;
          ecalohitsum_known += eCell*w;
        } // end loop over hits -> MCParticle links
      } // end loop over cluster hits
      ecalohitsum_unknown = ecalohitsum - ecalohitsum_known;
      debug() << "Calo energy in hit: all/known/unknown = " << ecalohitsum << " " << ecalohitsum_known << " " << ecalohitsum_unknown << endmsg;

      double sumw=0.0; // for debug
      for (const auto &mcp: *mcparticles) {
        // create cluster<->mc particle associations
        if (mcpEnergy[mcp.id().index]>0) {
          auto link = clusterMCParticleLinkCollection->create();
          link.setFrom(cluster);
          link.setTo(mcp);
          double w = mcpEnergy[mcp.id().index] / ecalohitsum;
          debug() << "Link with weight " << w << " set to particle " << mcp.id()  << " with pdg = " << mcp.getPDG()  << " , energy = " << mcp.getEnergy()  << endmsg;
          sumw += w;
          link.setWeight(w);
        }
      }
      debug() << "Sum of weights for this cluster = " << sumw << endmsg;
    } // end loop over clusters
  } // end loop over cluster collections

  debug() << "Finished linking clusters to MC particles" << endmsg;

  return StatusCode::SUCCESS;
}

StatusCode CreateTruthLinks::finalize() {
  for (size_t ih = 0; ih < m_cell_hit_linkCollectionHandles.size(); ih++)
    delete m_cell_hit_linkCollectionHandles[ih];

  return Gaudi::Algorithm::finalize();
}
