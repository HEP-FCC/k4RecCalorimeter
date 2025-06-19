#include "CreateHitTruthLinks.h"

// edm4hep
// #include "edm4hep/CalorimeterHit.h"

DECLARE_COMPONENT(CreateHitTruthLinks)

CreateHitTruthLinks::CreateHitTruthLinks(const std::string& name, ISvcLocator* svcLoc)
    : Gaudi::Algorithm(name, svcLoc) {
  //declareProperty("hits", m_hits, "Hits collection (input)");
  //declareProperty("cells", m_cells, "Cell collection (input)");
  declareProperty("mcparticles", m_mcparticles, "MC particles collection (input)");
  //declareProperty("cell_hit_links", m_cell_hit_links, "The links between cells and hits (input)");
  declareProperty("cell_mcparticle_links", m_cell_mcparticle_links, "The links between cells and MC particles (output)");
}

CreateHitTruthLinks::~CreateHitTruthLinks() { }

StatusCode CreateHitTruthLinks::initialize() {
  StatusCode sc = Gaudi::Algorithm::initialize();
  if (sc.isFailure())
    return sc;

  // create handles for input collections
  for (const auto& col : m_hitCollections) {
    debug() << "Creating handle for input hit (SimCalorimeterHit) collection : " << col << endmsg;
    try {
      m_hitCollectionHandles.push_back(
          new k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>(col, Gaudi::DataHandle::Reader, this));
    } catch (...) {
      error() << "Error creating handle for input collection: " << col << endmsg;
      return StatusCode::FAILURE;
    }
  }

  for (const auto& col : m_cellCollections) {
    debug() << "Creating handle for input cell (CalorimeterHit) collection : " << col << endmsg;
    try {
      m_cellCollectionHandles.push_back(
          new k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection>(col, Gaudi::DataHandle::Reader, this));
    } catch (...) {
      error() << "Error creating handle for input collection: " << col << endmsg;
      return StatusCode::FAILURE;
    }
  }

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

  info() << "CreateHitTruthLinks initialized" << endmsg;

  return StatusCode::SUCCESS;
}

StatusCode CreateHitTruthLinks::execute(const EventContext&) const {


  // create cell<->particle links
  edm4hep::CaloHitMCParticleLinkCollection* caloHitMCParticleLinkCollection = new edm4hep::CaloHitMCParticleLinkCollection();

  // retrieve MC particles
  const edm4hep::MCParticleCollection* mcparticles = m_mcparticles.get();

  // loop over input sim hit collections
  for (size_t ih = 0; ih < m_hitCollectionHandles.size(); ih++) {
    const edm4hep::SimCalorimeterHitCollection* simCalorimeterHits = m_hitCollectionHandles[ih]->get();
    debug() << "Processing input hit collection " << ih << " : " << m_hitCollectionHandles[ih]->objKey() << endmsg;

    const int nsimhits = (int) simCalorimeterHits->size();
    debug() << "Input Hit collection size: " << nsimhits << endmsg;

    // retrieve corresponding calo<->sim hit collection
    const edm4hep::CaloHitSimCaloHitLinkCollection* caloHitSimCaloHitLinks = m_cell_hit_linkCollectionHandles[ih]->get();

    // Loop over the G4 hits, find the associated calo hits, then loop over the G4 hit contribution to calculate contribution of each MCParticle to the calo hit
    for (const auto &simHit : *simCalorimeterHits){

      std::map< size_t , double > simHitMapEnergy ;  //  counts total energy for every MCParticle
      int nhits = 0;
      for (const auto &assoc : *caloHitSimCaloHitLinks)
      {
        if (assoc.getTo() == simHit) {
          nhits++;
          if (nhits>1) {
            error() << "Sim hit with more than one calo hit!!!" << endmsg;
            return StatusCode::FAILURE;
          }
          auto caloHit = assoc.getFrom();
          debug() << "Found associated calo hit" << endmsg;

          // extract digitised over deposited energy
          // assuming here there is at most one simhit per calo hit (0 in the case of noise hit)
          double calib_factor = caloHit.getEnergy()/simHit.getEnergy();

          // now loop over truth contributions
          int k=0;
          for (const auto &simHitContrib: simHit.getContributions())
          {

            const edm4hep::MCParticle& mcp = simHitContrib.getParticle();
            double e  = simHitContrib.getEnergy() * calib_factor;

            debug() << "initial true contributor "<< k << " id " << mcp.id()
                    <<" (with E = " << mcp.getEnergy() << " and pdg " << mcp.getPDG() << " ) e hit: " << e << std::endl;
            {
  /*
            if ( remap_as_you_go.find(mcp) != remap_as_you_go.end() ) {
              // very first condition: I already know what to do from some earlier hit created by mcp.
              mcp=remap_as_you_go.find(mcp)->second;
              streamlog_out( DEBUG3 ) << "     Hit " << simHit->id() << " / "
                                  << caloHit->id() << " attributed to " << mcp.id() <<
                          " because it's origin mcp has already been treated : originator case 0 " <<std::endl;
            } else {
              MCParticle* mother = 0;
              MCParticle* this_Kid = mcp ;     // ... and a true particle !

              //Particle gun particles dont have parents, but are "created in the simulation" and have genStat 0
              if ( mcp. getGeneratorStatus() == 0 && ( mcp.getParents().empty() == false ) ) {

                // not from generator, find which true particle this
                // hit should really be attributed to, by tracking back the history.
                // (the other case, ie. if the if above is false is "case 1")

                // Two cases to treat:
                //   For some reason, the hit is not attributed to the
                //   incomming particle, but to some particle in the shower.
                //   Just track back to the incomming particle, which might (usually)
                //   be a gen-stat 1 particle from the main vertex. It can also
                //   be from a decay or interaction in the tracking made by
                //   Geant : genstat=0, mother isDecayedInTracker, (or
                //   the production vertex is in the track, but the mother
                //   continues on - ticky case, see comment futher down), or a
                //   decyed particle from the generator (genstat 2) that
                //   hits the calo before decaying (lambdas, K^0_S)

                //   Or: for back-scatterers the calorimiter hit is attributed to the
                //   the last particle *even if this particle both
                //   started and ended in the tracker* !!!! Then we back-track
                //   untill we find a particle which at least started inside
                //   the calorimeter, and then go on backtracking as above.
                //   This case is triggered by the particle linked to the
                //   hit being DecayedInTracker, hence the case where a
                //   back-scatter actually ends in the calorimeter is
                //   treated as the first case.



                streamlog_out( DEBUG2 ) << "        simHit " <<simHit->id()<<","<< j <<
                " not created by generator particle. backtracking ..."
                << std::endl;
                streamlog_out( DEBUG2 ) <<"          "<<mcp.id()<<" gs "<<mcp.getGeneratorStatus()<<
                " dint "<<mcp.isDecayedInTracker()<<
                " bs "<<mcp.isBackscatter()<<
                " ndi "<<mcp.vertexIsNotEndpointOfParent()<<
                " npar "<<mcp.getParents().size()<<
                " pdg "<<mcp.getPDG()<<
                " "<< mcp.getVertex()[0]<<
                " "<<mcp.getVertex()[1]<<
                " "<<mcp.getVertex()[2]<<
                " "<< mcp.getEndpoint()[0]<<
                " "<<mcp.getEndpoint()[1]<<
                " "<<mcp.getEndpoint()[2]<<std::endl;



                mother= dynamic_cast<MCParticle*>(mcp.getParents()[0]);


                if ( !this_Kid->isBackscatter() &&  mother!= 0 &&
                      mother->getParents().size()>0 &&
                      mother->getGeneratorStatus() ==0 &&
                      !mother->isDecayedInTracker() ) {

            streamlog_out( DEBUG2 ) << "        goes into originator loop " << std::endl;
                }
          while ( !this_Kid->isBackscatter() && mother!= 0 &&   mother->getParents().size()>0 &&
                        mother->getGeneratorStatus() ==0 &&
                        !mother->isDecayedInTracker() ) {
                  // back-track as long as there is a non-generator
                  // mother, or the mother decayed in the tracker
                  // (=> the kid is the particle entering the calorimeter.)

                  // case shower-particle
            streamlog_out( DEBUG1 ) <<"          in originator loop " << std::endl;

      if ( this_Kid->vertexIsNotEndpointOfParent() != _invertedNonDestructiveInteractionLogic) {
        MCParticle* oma=dynamic_cast<MCParticle*>(mother->getParents()[0]);
                    if ( oma->isDecayedInTracker() ) {
                      streamlog_out( DEBUG1 ) <<"          break out : gandmother "<<oma->id()<<
                                                " gs "<<oma->getGeneratorStatus()<<
                                                " dint "<<oma->isDecayedInTracker()<<
                                                " bs "<<oma->isBackscatter()<<
                                                " ndi "<<oma->vertexIsNotEndpointOfParent()<<
                                                " npar "<<oma->getParents().size()<<
                                                " pdg "<<oma->getPDG()<<std::endl;
                      break ;
                    }
      }
                  this_Kid=mother ;
                  mother= dynamic_cast<MCParticle*>(mother->getParents()[0]); // (assume only one...)

                  streamlog_out( DEBUG1 ) <<"          shower-part mother "<<mother->id()<<
                  " gs "<<mother->getGeneratorStatus()<<
                  " dint "<<mother->isDecayedInTracker()<<
                  " bs "<<mother->isBackscatter()<<
                  " ndi "<<mother->vertexIsNotEndpointOfParent()<<
                  " npar "<<mother->getParents().size()<<
                  " pdg "<<mother->getPDG()<<std::endl;



                }

                // Further treatment (basically determining if it this_Kid or mother that enetered the
                // calorimeter) based on why we left the while-loop

                // here one of the while conditions is false, ie. at least one of
                // " kid is back-scatter", "no mother" , "mother has no parents", "mother is from
                // generator", or "mother did decay in tracker" must be true. We know that this_Kid
                // fulfills all the while-conditions except the first: obvious if at least one iteration of
                // the loop was done, since it was the mother in the previous iteration,
                // but also true even if the loop wasn't transversed, due to the conditions
                // to at all enter this block of code (explicitly must be simulator particle with
                // mother, implicitly must have ended in the calo, since it did make calo hits.)
          if (this_Kid->isBackscatter() ) {
                  // case 2: Kid is back-scatterer. It has thus started in a calo, and entered
                  //  from there into the tracking volume, and did cause hits after leaving the tracker
                  //  volume again ->  this_Kid started before the calo, and is the one
                  remap_as_you_go[mcp]=this_Kid;
                  mcp=this_Kid;
                  streamlog_out( DEBUG3 ) << "     Hit " << simHit->id() << " / "
                                  << caloHit->id() << " attributed to kid " << mcp.id() <<
                          " because it's origin is a back-scatter : originator case 2 " <<std::endl;
                  streamlog_out( DEBUG2 ) <<"          "<<this_Kid->id()<<
                  " gs "<<this_Kid->getGeneratorStatus()<<
                  " dint "<<this_Kid->isDecayedInTracker()<<
                  " bs "<<this_Kid->isBackscatter()<<
                  " npar "<<this_Kid->getParents().size()<<
                  " pdg "<<this_Kid->getPDG()<<std::endl;


                } else if ( mother->isDecayedInTracker() ) { // the clear-cut case:
                  remap_as_you_go[mcp]=this_Kid;
                  mcp=this_Kid; // this_Kid started before the calo, and is the one
                                // the hit should be attributed to


                  streamlog_out( DEBUG3 ) << "     Hit " << simHit->id() << " / "
                                  << caloHit->id() << " attributed to kid " << mcp.id() <<
                          " because it's origin is in tracker : originator case 3 " <<std::endl;


                } else { // the other three cases, ie. one or sveral of "no mother", "no grand-mother",
                        //  "generator particle" + that we know that "mother decayed in calo"
                  if ( mother == 0 ) { // ... which of course implies no grand-mother, and no gen stat
                                      // of the mother as well -> should not be possible !

                    streamlog_out( WARNING ) << "  MCparticle " << this_Kid->id() <<
                    " is a simulation particle, created in the calorimeter by nothing . "<< std::endl;


                    remap_as_you_go[mcp]=this_Kid;
                    mcp=this_Kid; // can't do better than that.

                  } else { // here we know: "mother exists, but decayed in calo". In addition, two posibilities:
                          // mother is generator particle, or there was no grand-parents. One or both
                          // must be true here. Here it gets complicated, because what we want to know is
                          // whether this_Kid started in the tracker or not. Unluckily, we don't know that
                          // directly, we only know where the mother ended. IF ithe mother ended in the tracker,
                          // there is no problem, and has already been treated, but if it ended in the calo, it is
                          // still possible that this_Kid came from a "non-destructive interaction" with the tracke-detector
                          // material. This we now try to figure out.

        if (   this_Kid->vertexIsNotEndpointOfParent() == _invertedNonDestructiveInteractionLogic ) {
            // This bizare condition is due to a bug in LCIO (at least for the DBD samples. this_Kid->vertexIsNotEndpointOfParent()
                        // should be true in the "non-destructive interaction", but it isn't: actually it is "false", but is "true" for the
                        // for particles that *do* originate at the end-vertex of their parent. This is a bug in Mokka.
                        // Hence the above ensures that there was NO "non-destructive interaction", and it is clear this_Kid was created at the end-point
                        // of the mother. The mother is either a generator particle (to be saved), or the "Eve" of the decay-chain (or both).
                        // So mother is the one to save and assign the hits to:
                      remap_as_you_go[mcp]=mother;
                      mcp=mother;   // mother started at ip, and reached the calo, and is the one
                                    // the hit should be attributed to.

                      streamlog_out( DEBUG3 ) << "     Hit " << simHit->id() << " / "
                                  << caloHit->id() << " attributed to mother " << mcp.id() <<
                          " because it is a generator particle or started in tracker : originator case 4 " <<std::endl;

                    } else { // here we DO have a "non-destructive interaction". Unluckily, we cant directly know if this took place in the
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
                      int oma_in_calo = 0;
                      double rdist =0.;
          int has_bs = 0;

          for ( unsigned kkk=0 ; kkk < mother->getDaughters().size() ; kkk++ ) {
                        MCParticle* sister = dynamic_cast<MCParticle*>(mother->getDaughters()[kkk]);
                        if ( sister == this_Kid ) continue;
            if (  abs(sister->getVertex()[0]-this_Kid->getVertex()[0]) > 0.1 ||
                              abs(sister->getVertex()[1]-this_Kid->getVertex()[1]) > 0.1 ||
                              abs(sister->getVertex()[2]-this_Kid->getVertex()[2]) > 0.1 ) continue;  // must check that it is the same vertex:
                                                                                                      // several "non-destructive interactions" can
                                                                                                      // take place (think delta-rays !)
            if ( sister->isBackscatter()) {
                          has_bs = 1 ;
                          break ;
            } else if ( sister->isDecayedInTracker() ) {
                          starts_in_tracker = 1 ;
                          break ;
                        }
                        // any pi0:s at all ? (it doesn't matter that we break at the two cases above,
                        // because if we do, it doesn't matter if there are
                        // pi0 sisters or not !)
                        if ( sister->getPDG() == 111 ) {
                          has_pi0 = 1 ;
                        }
                      }
                      // if not already clear-cut, calculate distance vertext to mother end-point
                      if ( starts_in_tracker != 1 && has_bs != 1 && has_pi0 != 1 && oma_in_calo != 1 ) {
                        rdist=sqrt(pow(mother->getEndpoint()[0]-this_Kid->getVertex()[0],2)+
                                    pow(mother->getEndpoint()[1]-this_Kid->getVertex()[1],2)+
                                    pow(mother->getEndpoint()[2]-this_Kid->getVertex()[2],2));
                        if ( mother->getParents().size() != 0 ) {
              MCParticle* oma=dynamic_cast<MCParticle*>(mother->getParents()[0]);
                          if ( oma->isDecayedInCalorimeter() ) {
                            oma_in_calo = 1 ;
          streamlog_out( DEBUG1 ) << "          grandmother in calo " << std::endl;
                          }
                        }
                      }
                      streamlog_out( DEBUG1 ) << "          " <<  starts_in_tracker << " " << has_pi0 << " " <<  has_bs << " " << oma_in_calo << std::endl;
          //                    }
                      if ( starts_in_tracker == 1 ) { // this_Kid is a clear-cut hit-originator

                        remap_as_you_go[mcp]=this_Kid;
                        mcp=this_Kid; // this_Kid started before the calo, and is the one
                                      // the hit should be attributed to

                        streamlog_out( DEBUG3 ) << "     Hit " << simHit->id() << " / "
                                  << caloHit->id() << " attributed to kid " << mcp.id() <<
                          " because it's origin could be deduced to be in tracker : originator case 5 " <<std::endl;

            streamlog_out( DEBUG2 ) << "        Details of case 5: "
                                              << this_Kid->getVertex()[0] << " "
                << this_Kid->getVertex()[1] << " " << this_Kid->getVertex()[2] << " "
                                              << mother->getEndpoint()[0] << " "
                << mother->getEndpoint()[1] << " " << mother->getEndpoint()[2] << " "
                                              << mother->getGeneratorStatus() <<  " "
                                              << this_Kid->vertexIsNotEndpointOfParent() << std::endl;
                        streamlog_out( DEBUG1 ) << " starts in tracker " <<  std::endl;

          } else if ( has_pi0 != 0 || has_bs != 0 || oma_in_calo != 0 ) { // clear-cut case of this_Kid starting in the calo.
                                                                  // We do know that the mother
                                                                  // is a generator particle and/or the "Eve" of the cascade,
                                                                  // so we should attribute hits to the
                                                                  // mother and save it.

                        remap_as_you_go[mcp]=mother;
                        mcp=mother;   // mother started at ip, and reached the calo, and is the one
                                      // the hit should be attributed to.

                        streamlog_out( DEBUG3 ) << "     Hit " << simHit->id() << " / "
                                  << caloHit->id() << " attributed to mother " << mcp.id() <<
                                  " because it's origin could be deduced to be in tracker : originator case 6 " <<std::endl;
            streamlog_out( DEBUG2 ) << "        Case 6 details: kid starts in calo " << " "
                  << this_Kid->getVertex()[1] << " " << this_Kid->getVertex()[2] << " " << rdist <<  " "
                  << has_pi0  << " " << has_bs << " " <<  oma_in_calo << std::endl;


          } else { // un-clear case: no pi0 nor back-scatteres among the sisters to help to decide.
                              // Use distance this_Kid-startpoint to mother-endpoint. We know that the latter is
                              // in the calo, so if this is small, guess that the start point of this_Kid is
                              // also in the calo. Calos are dense, so typically in the case the "non-destructive interaction"
                              // is in the calo, one would guess  that the distance is small, ie. large distance ->
                              // unlikely that it was in the calo.

            if ( rdist > 200. ) { // guess "non-destructive interaction" not in calo -> this_Kid is originator

                          remap_as_you_go[mcp]=this_Kid;
                          mcp=this_Kid; // this_Kid started before the calo, and is the one
                                      // the hit should be attributed to
                          streamlog_out( DEBUG3 ) << "     Hit " << simHit->id() << " / "
                                  << caloHit->id() << " attributed to kid " << mcp.id() <<
                                  " because it's origin is guessed to be in tracker : originator case 7 " <<std::endl;
              streamlog_out( DEBUG2 ) << "        Case 7 details: guess kid starts in tracker " << " "
                << this_Kid->getVertex()[1] << " " << this_Kid->getVertex()[2] << " " << rdist <<  " "
                                                    << has_pi0  << " " << std::endl;

                        } else { // guess in calo -> mother is originator
                          remap_as_you_go[mcp]=mother;
                          mcp=mother;   // mother started at ip, and reached the calo, and is the one
                                      // the hit should be attributed to.
                          streamlog_out( DEBUG3 ) << "     Hit " << simHit->id() << " / "
                                  << caloHit->id() << " attributed to mother " << mcp.id() <<
                                  " because it's origin is guessed in tracker : originator case 8 " <<std::endl;
              streamlog_out( DEBUG2 ) << "        Case 8 details:  guess kid starts in calo " << " "
                << this_Kid->getVertex()[1] << " " << this_Kid->getVertex()[2] << " " << rdist <<  " "
                                                    << has_pi0  << " " << std::endl;
                          }

                      }
                    }
                  }
                }

              } else {
                // first case: the hit generating mcp itself already fulfills the firest
                // criterium, ie. it is a generator particle
                streamlog_out( DEBUG3 ) << "     Hit " << simHit->id() << " / "
                                  << caloHit->id() << " attributed to " << mcp.id() <<
                          " because it's origin is a generator particle : originator case 1 " <<std::endl;

                remap_as_you_go[mcp]=mcp;
        } // genstat if - then - else
            }
        */
            }

            debug() << "    Final assignment for contribution " << k << " to " << simHit.id() << " / "
                    << caloHit.id() << " : " << mcp.id() << endmsg;
            debug() << "      gs "<<mcp.getGeneratorStatus()<<
                      " dint " <<mcp.isDecayedInTracker()<<
                      "  bs "<< mcp.isBackscatter()<<
                      " npar "<<mcp.getParents().size()<<
                      " pdg "<<mcp.getPDG()<<
                      " " <<mcp.getEndpoint()[0] <<
                      " " <<mcp.getEndpoint()[1] <<
                      " " <<mcp.getEndpoint()[2] << endmsg;

            simHitMapEnergy[ mcp.id().index ] += e;
          } // end loop over contributions
          double sumw=0.0; // for debug
          for (const auto &mcp: *mcparticles) {
            // create calo hit<->mc particle associations
            if (simHitMapEnergy[mcp.id().index]>0) {
              auto link = caloHitMCParticleLinkCollection->create();
              link.setFrom(caloHit);
              link.setTo(mcp);
              double w = simHitMapEnergy[mcp.id().index]/ caloHit.getEnergy();
              sumw += w;
              link.setWeight(w);
            }
          }
          debug() << "Sum of weights for this calo hit = " << sumw << endmsg;
        } // end if (assoc.getTo() == simHit)
      } // end loop over calo hit <-> sim calo hit links
      if (nhits==0) {
        error() << "Sim hit with no associated  calo hit!!!" << endmsg;
        return StatusCode::FAILURE;
      }
    } // end loop over sim calo hits
  } // end loop over sim calo hit collections

  // push the CaloHitCollection to event store
  m_cell_mcparticle_links.put(caloHitMCParticleLinkCollection);
  debug() << "Finished linking calo hits to MC particles" << endmsg;

  return StatusCode::SUCCESS;
}

StatusCode CreateHitTruthLinks::finalize() {
  for (size_t ih = 0; ih < m_hitCollectionHandles.size(); ih++)
    delete m_hitCollectionHandles[ih];

  for (size_t ih = 0; ih < m_cellCollectionHandles.size(); ih++)
    delete m_cellCollectionHandles[ih];

  for (size_t ih = 0; ih < m_cell_hit_linkCollectionHandles.size(); ih++)
    delete m_cell_hit_linkCollectionHandles[ih];

  return Gaudi::Algorithm::finalize();
}
