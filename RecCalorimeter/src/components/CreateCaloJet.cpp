#include "CreateCaloJet.h"

// std
#include <vector>
#include <math.h>

// EDM4hep
#include "edm4hep/ClusterCollection.h"

DECLARE_COMPONENT(CreateCaloJet)

CreateCaloJet::CreateCaloJet(const std::string& name, ISvcLocator* svcLoc) : Transformer(name, svcLoc,
                    KeyValue("InputCollection", "CorrectedCaloClusters"),
                    KeyValue("OutputCollection", "Jets")) {
  declareProperty("JetAlg", m_jetAlg, "Name of jet clustering algorithm");
  declareProperty("JetRadius", m_jetRadius, "Jet clustering radius");
  declareProperty("MinPt", m_minPt, "Minimum pT for saved jets");
  declareProperty("isExclusiveClustering", m_isExclusive, "1 if exclusive, 0 if inclusive");

}


StatusCode CreateCaloJet::initialize() {
  clusterer = new k4::recCalo::ClusterJet(m_jetAlg, m_jetRadius, m_isExclusive, m_minPt);


  return clusterer->initialize();
}

edm4hep::ReconstructedParticleCollection  CreateCaloJet::operator()(const edm4hep::ClusterCollection& input) const{
  std::vector<fastjet::PseudoJet> clustersPJ;
  int i=0;

  // Need to convert cluster position in space to the momentum
  for(auto cluster: input){
    const edm4hep::Vector3f position = cluster.getPosition();
    float x = position.x;
    float y = position.y;
    float z = position.z;
    double theta = acos(sqrt(z*z / (x*x + y*y + z*z)));
    double eta = -log(tan(theta / 2.));
    double phi = atan2(y, x);
    double pT = cluster.getEnergy()* sqrt( (x*x + y*y) / (x*x + y*y + z*z));

    // Take clusters to be massless
    fastjet::PseudoJet clusterPJ;
    clusterPJ.reset_momentum_PtYPhiM(pT, eta, phi, 0);
    // Attach the cluster index to the pseudojet
    clusterPJ.set_user_info(new k4::recCalo::ClusterInfo(i));
    clustersPJ.push_back(clusterPJ);
    i++;
  }


  std::vector<fastjet::PseudoJet> inclusiveJets = clusterer->cluster(clustersPJ);


  edm4hep::ReconstructedParticleCollection edmJets = edm4hep::ReconstructedParticleCollection();
  //Add a reconstructed particle for each jet
  for(auto cjet : inclusiveJets){
    edm4hep::MutableReconstructedParticle jet;
    jet.setMomentum(edm4hep::Vector3f(cjet.px(), cjet.py(), cjet.pz()));
    jet.setEnergy(cjet.e());
    jet.setMass(cjet.m());

    // Also add the clusters that were used to make the jets
    std::vector<fastjet::PseudoJet> constits = cjet.constituents();
    for(auto constit : constits){
      int index = constit.user_info<k4::recCalo::ClusterInfo>().index();
      jet.addToClusters((input)[index]);
    }
    edmJets.push_back(jet);
  }

  return edmJets;
}


