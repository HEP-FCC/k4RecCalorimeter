// Gaudi
#include "Gaudi/Property.h"

// k4FWCore
#include "k4FWCore/Transformer.h"

// std
#include <math.h>
#include <vector>

// EDM4hep
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"

#include "ClusterJet.h"

/** @class CreateCaloJet
 * k4RecCalorimeter/RecFCCeeCalorimeter/src/components/CreateCaloJet.h
 *
 *  Tool for reconstructing jets from calorimeter clusters.
 *  It takes as input a ClusterCollection, and it outputs a
 * ReconstructedParticleCollection with the jets.
 *
 *  JetAlg: A string corresponding to the jet algorithm for clustering
 *  JetRadius: The radius of the jets being clustered
 *  MinPt: The pT threshold below which jets are ignored
 *  isExclusiveClustering: 1 if jets should use an exclusive clustering, 0
 * otherwise
 *
 *  @author Jennifer Roloff
 *  @date 2024-7
 */

struct CreateCaloJet final
    : k4FWCore::Transformer<edm4hep::ReconstructedParticleCollection(
          const edm4hep::ClusterCollection &)> {
  CreateCaloJet(const std::string &name, ISvcLocator *svcLoc)
      : Transformer(
            name, svcLoc,
            {KeyValues("InputClusterCollection", {"CorrectedCaloClusters"})},
            {KeyValues("OutputJetCollection", {"Jets"})}) {}

  StatusCode initialize() {
    clusterer = new k4::recCalo::ClusterJet(m_jetAlg, m_jetRadius,
                                            m_isExclusive, m_minPt);

    return clusterer->initialize();
  }

  edm4hep::ReconstructedParticleCollection
  operator()(const edm4hep::ClusterCollection &input) const override {
    std::vector<fastjet::PseudoJet> clustersPJ;
    int i = 0;

    // Need to convert cluster position in space to the momentum
    for (auto cluster : input) {
      const edm4hep::Vector3f position = cluster.getPosition();
      float x = position.x;
      float y = position.y;
      float z = position.z;
      double theta = acos(sqrt(z * z / (x * x + y * y + z * z)));
      double eta = -log(tan(theta / 2.));
      double phi = atan2(y, x);
      double pT =
          cluster.getEnergy() * sqrt((x * x + y * y) / (x * x + y * y + z * z));

      // Take clusters to be massless
      fastjet::PseudoJet clusterPJ;
      clusterPJ.reset_momentum_PtYPhiM(pT, eta, phi, 0);
      // Attach the cluster index to the pseudojet
      clusterPJ.set_user_info(new k4::recCalo::ClusterInfo(i));
      clustersPJ.push_back(clusterPJ);
      i++;
    }

    std::vector<fastjet::PseudoJet> inclusiveJets =
        clusterer->cluster(clustersPJ);

    edm4hep::ReconstructedParticleCollection edmJets =
        edm4hep::ReconstructedParticleCollection();
    // Add a reconstructed particle for each jet
    for (auto cjet : inclusiveJets) {
      edm4hep::MutableReconstructedParticle jet;
      jet.setMomentum(edm4hep::Vector3f(cjet.px(), cjet.py(), cjet.pz()));
      jet.setEnergy(cjet.e());
      jet.setMass(cjet.m());

      // Also add the clusters that were used to make the jets
      std::vector<fastjet::PseudoJet> constits = cjet.constituents();
      for (auto constit : constits) {
        int index = constit.user_info<k4::recCalo::ClusterInfo>().index();
        jet.addToClusters((input)[index]);
      }
      edmJets.push_back(jet);
    }

    return edmJets;
  }

private:
  Gaudi::Property<std::string> m_jetAlg{this, "JetAlg", "antikt",
                                        "Name of jet clustering algorithm"};
  Gaudi::Property<double> m_jetRadius{this, "JetRadius", 0.4,
                                      "Jet clustering radius"};
  Gaudi::Property<double> m_minPt{this, "MinPt", 10,
                                  "Minimum pT for saved jets"};
  Gaudi::Property<int> m_isExclusive{this, "IsExclusiveClustering", 0,
                                     "1 if exclusive, 0 if inclusive"};

  k4::recCalo::ClusterJet *clusterer;
};

DECLARE_COMPONENT(CreateCaloJet)
