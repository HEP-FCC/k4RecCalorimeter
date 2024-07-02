#include "FilterTruthParticles.h"

// std
#include <vector>
#include <math.h>



DECLARE_COMPONENT(FilterTruthParticles)

FilterTruthParticles::FilterTruthParticles(const std::string& name, ISvcLocator* svcLoc) : Transformer(name, svcLoc,
                    KeyValue("InputCollection", "MCParticles"),
                    KeyValue("OutputCollection", "FilteredMCParticles")) {

  declareProperty("SaveMuons", m_useMuons, "Bool for if muons are stored");
  declareProperty("MinPt", m_minPt, "Minimum pT for filtered particles");
}



colltype_out  FilterTruthParticles::operator()(const colltype_in& input) const{

  edm4hep::MCParticleCollection filteredParticles = edm4hep::MCParticleCollection();

  for(auto particle: input){

    // Only want to include final state particles (before GEANT simulation)
    if(particle.isCreatedInSimulation()) continue;
    if(particle.getGeneratorStatus() != 1) continue;

    // Sometimes we want to ignore muons for jet reconstruction
    if(m_useMuons && std::abs(particle.getPDG()) == 13) continue; 

    // We may want to impose a minimum pT requirement
    double input_particle_pt = std::sqrt(particle.getMomentum().x * particle.getMomentum().x + particle.getMomentum().y * particle.getMomentum().y);
    if(input_particle_pt < m_minPt) continue;

    edm4hep::MutableMCParticle filteredParticle = particle.clone();
    filteredParticles.push_back(filteredParticle);
  }


  return filteredParticles;
}



