// Gaudi
#include "Gaudi/Property.h"
#include <GaudiKernel/StatusCode.h>

// k4FWCore
#include "k4FWCore/Transformer.h"

// EDM4hep
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/MCRecoParticleAssociationCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"

// k4RecCalorimeter
#include "RecCaloCommon/ClusterJet.h"

// std
#include <string>
#include <vector>

/** @struct CreateTruthJet
 k4RecCalorimeter/RecCalorimeter/src/components/CreateTruthJet.h
 *
 *  Tool for reconstructing jets from truth particles.
 *  It takes as input an MCParticleCollection, and it outputs a
 ReconstructedParticleCollection with the jets.
 *
 *  JetAlg: A string corresponding to the jet algorithm for clustering
 *  JetRadius: The radius of the jets being clustered
 *  MinPt: The pT threshold below which jets are ignored
 *  isExclusiveClustering: 1 if jets should use an exclusive clustering, 0
 otherwise
 *  outputAssociation: Name of the output association collection to link jets to
 their constituents
 *
 *  @author Jennifer Roloff
 *  @date 2024-7
 */

struct CreateTruthJet final
    : k4FWCore::MultiTransformer<
          std::tuple<edm4hep::ReconstructedParticleCollection,
                     edm4hep::MCRecoParticleAssociationCollection>(
              const edm4hep::MCParticleCollection &)> {

  CreateTruthJet(const std::string &name, ISvcLocator *svcLoc)
      : MultiTransformer(
            name, svcLoc,
            {KeyValues("InputMCParticleCollection", {"MCParticles"})},
            {KeyValues("OutputJetCollection", {"TruthJets"}),
             KeyValues("OutputAssociationsCollection",
                       {"TruthJetParticleAssociations"})}) {}

  StatusCode initialize() override {
    m_clusterer = new k4::recCalo::ClusterJet(m_jetAlg, m_jetRadius,
                                              m_isExclusive, m_minPt);
    if (!m_clusterer->initialize()) {
      return StatusCode::FAILURE;
    }

    return StatusCode::SUCCESS;
  }

  std::tuple<edm4hep::ReconstructedParticleCollection,
             edm4hep::MCRecoParticleAssociationCollection>
  operator()(const edm4hep::MCParticleCollection &input) const override {

    std::vector<fastjet::PseudoJet> clustersPJ;
    int i = 0;
    for (auto particle : input) {
      fastjet::PseudoJet clusterPJ(
          particle.getMomentum().x, particle.getMomentum().y,
          particle.getMomentum().z, particle.getEnergy());
      clusterPJ.set_user_info(new k4::recCalo::ClusterInfo(i));
      clustersPJ.push_back(clusterPJ);
      i++;
    }

    std::vector<fastjet::PseudoJet> inclusiveJets =
        m_clusterer->cluster(clustersPJ);

    auto edmJets = edm4hep::ReconstructedParticleCollection();
    auto assoc = edm4hep::MCRecoParticleAssociationCollection();

    for (auto cjet : inclusiveJets) {
      edm4hep::MutableReconstructedParticle jet;
      jet.setMomentum(edm4hep::Vector3f(cjet.px(), cjet.py(), cjet.pz()));
      jet.setEnergy(cjet.e());
      jet.setMass(cjet.m());

      std::vector<fastjet::PseudoJet> constits = cjet.constituents();

      for (auto constit : constits) {
        int index = constit.user_info<k4::recCalo::ClusterInfo>().index();

        auto association = assoc.create();
        association.setRec(jet);
        association.setSim((input)[index]);
      }

      edmJets.push_back(jet);
    }

    return std::make_tuple(std::move(edmJets), std::move(assoc));
  }

  StatusCode finalize() override {
    delete m_clusterer;

    return StatusCode::SUCCESS;
  }

private:
  Gaudi::Property<std::string> m_jetAlg{this, "JetAlg", "antikt",
                                        "Name of jet clustering algorithm"};
  Gaudi::Property<double> m_jetRadius{this, "JetRadius", 0.4,
                                      "Jet clustering radius"};
  Gaudi::Property<double> m_minPt{this, "MinPt", 1.,
                                  "Minimum pT for saved jets"};
  Gaudi::Property<bool> m_isExclusive{this, "isExclusiveClustering", false,
                                      "true if exclusive, false if inclusive"};
  k4::recCalo::ClusterJet *m_clusterer = nullptr;
};

DECLARE_COMPONENT(CreateTruthJet)
