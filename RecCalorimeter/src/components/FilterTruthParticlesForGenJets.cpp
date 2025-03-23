// Gaudi
#include "Gaudi/Property.h"

// k4FWCore
#include "k4FWCore/Transformer.h"

// EDM4hep
#include "edm4hep/MCParticleCollection.h"

// FastJet
#include "fastjet/JetDefinition.hh"

/** @struct FilterTruthParticlesForGenJets
 * k4RecCalorimeter/RecCalorimeter/src/components/FilterTruthParticlesForGenJets.h
 *
 *  Tool to filter truth particles based on their status code.
 *  Only detector stable particles are stored in a new MCParticleCollection.
 *  This can provide inputs to truth jet clustering.
 *
 *  MinPt: the minimum pT of particles that are stored
 *  KeepMuons: true if muons should be included in the output collection, false
 * otherwise.
 *
 *  @author Jennifer Roloff
 *  @date 2024-7
 */

struct FilterTruthParticlesForGenJets final
    : k4FWCore::Transformer<edm4hep::MCParticleCollection(const edm4hep::MCParticleCollection&)> {

  FilterTruthParticlesForGenJets(const std::string& name, ISvcLocator* svcLoc)
      : Transformer(name, svcLoc, {KeyValues("InputMCParticleCollection", {"MCParticles"})},
                    {KeyValues("OutputMCParticleCollection", {"FilteredMCParticles"})}) {}

  edm4hep::MCParticleCollection operator()(const edm4hep::MCParticleCollection& input) const override {

    debug() << "Number of input particles: " << input.size() << endmsg;

    auto filteredParticles = edm4hep::MCParticleCollection();
    for (auto particle : input) {

      // Only want to include final state particles (before GEANT simulation)
      if (particle.isCreatedInSimulation())
        continue;
      if (particle.getGeneratorStatus() != 1)
        continue;

      // Sometimes we want to ignore muons for jet reconstruction
      if (!m_keepMuons && std::abs(particle.getPDG()) == 13)
        continue;

      // We may want to impose a minimum pT requirement
      double input_particle_pt = std::sqrt(particle.getMomentum().x * particle.getMomentum().x +
                                           particle.getMomentum().y * particle.getMomentum().y);
      if (input_particle_pt < m_minPt)
        continue;

      auto filteredParticle = particle.clone();
      filteredParticles.push_back(filteredParticle);
    }

    debug() << "Number of filtered particles: " << filteredParticles.size() << endmsg;

    return filteredParticles;
  }

private:
  Gaudi::Property<bool> m_keepMuons{this, "SaveMuons", 0, "Bool for if muons are kept"};
  Gaudi::Property<double> m_minPt{this, "MinPt", 0., "Minimum pT for filtered particles"};
};

DECLARE_COMPONENT(FilterTruthParticlesForGenJets)
