#include "CreateTruthJet.h"

// std
#include <vector>
#include <math.h>



DECLARE_COMPONENT(CreateTruthJet)

CreateTruthJet::CreateTruthJet(const std::string& name, ISvcLocator* svcLoc) : Transformer(name, svcLoc,
                    KeyValue("InputCollection", "MCParticles"),
                    KeyValue("OutputCollection", "TruthJets")) {
  declareProperty("JetAlg", m_jetAlg, "Name of jet clustering algorithm");
  declareProperty("JetRadius", m_jetRadius, "Jet clustering radius");
  declareProperty("MinPt", m_minPt, "Minimum pT for saved jets");
}



StatusCode CreateTruthJet::initialize() {
  clusterer = new ClusterJet(m_jetAlg, m_jetRadius, m_minPt);

  return StatusCode::SUCCESS;
}

colltype_out  CreateTruthJet::operator()(const colltype_in& input) const{
  std::vector<fastjet::PseudoJet> clustersPJ;
  int i=0;
  for(auto particle: input){
    fastjet::PseudoJet clusterPJ(particle.getMomentum().x, particle.getMomentum().y, particle.getMomentum().z, particle.getEnergy());
    clusterPJ.set_user_info(new ClusterInfo(i));
    clustersPJ.push_back(clusterPJ);
    i++;
  }

  std::vector<fastjet::PseudoJet> inclusiveJets = clusterer->cluster(clustersPJ);
  

  edm4hep::ReconstructedParticleCollection edmJets = edm4hep::ReconstructedParticleCollection();

  for(auto cjet : inclusiveJets){
    edm4hep::MutableReconstructedParticle jet;
    jet.setMomentum(edm4hep::Vector3f(cjet.px(), cjet.py(), cjet.pz()));
    jet.setEnergy(cjet.e());
    jet.setMass(cjet.m());

    std::vector<fastjet::PseudoJet> constits = cjet.constituents();
    for(auto constit : constits){
      edm4hep::MutableReconstructedParticle jetInput;
      jetInput.setMomentum(edm4hep::Vector3f(constit.px(), constit.py(), constit.pz()));
      jetInput.setEnergy(constit.e());
      jetInput.setMass(constit.m());

      int index = constit.user_info<ClusterInfo>().index();
      jetInput.setPDG((input)[index].getPDG());

      jet.addToParticles(jetInput);
    }

    edmJets.push_back(jet);
  }

  return edmJets;
}



