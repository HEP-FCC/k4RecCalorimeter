#include "TruthJetFCCee.h"

// std
#include <vector>
#include <math.h>



DECLARE_COMPONENT(TruthJetFCCee)

TruthJetFCCee::TruthJetFCCee(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
  declareProperty("MCParticleInput", m_inputParticles, "Handle for input particle collection");
  declareProperty("JetAlg", m_jetAlg, "Name of jet clustering algorithm");
  declareProperty("JetRadius", m_jetRadius, "Jet clustering radius");
  declareProperty("JetOutput", m_jetCollection, "Handle for output jet collection");
  declareProperty("MinPt", m_minPt, "Minimum pT for saved jets");
}



StatusCode TruthJetFCCee::initialize() {
  if (GaudiAlgorithm::initialize().isFailure()) return StatusCode::FAILURE;

  return StatusCode::SUCCESS;
}

StatusCode TruthJetFCCee::execute() {
  // Create output collections
  auto* edmJets = m_jetCollection.createAndPut();

  const edm4hep::MCParticleCollection *inParticles = m_inputParticles.get();

  fastjet::JetDefinition* jetDef = new fastjet::JetDefinition( m_jetAlgMap[m_jetAlg], m_jetRadius);


  std::vector<fastjet::PseudoJet> clustersPJ;
  int i=0;
  for(auto cluster: *inParticles){
    if(cluster.isCreatedInSimulation()) continue;
    if(cluster.getGeneratorStatus() != 1) continue;
    // TODO: decide if muons should be included in the truth jet clustering
    //if(std::abs(cluster.getPDG()) == 13) continue; // No muons
    fastjet::PseudoJet clusterPJ(cluster.getMomentum().x, cluster.getMomentum().y, cluster.getMomentum().z, cluster.getEnergy());
    clusterPJ.set_user_info(new ClusterInfo(i));
    clustersPJ.push_back(clusterPJ);
    i++;
  }


  fastjet::ClusterSequence clustSeq(clustersPJ, *jetDef);
  std::vector <fastjet::PseudoJet> inclusiveJets = fastjet::sorted_by_pt(clustSeq.inclusive_jets(m_minPt));
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
      jetInput.setPDG((*inParticles)[index].getPDG());

      //jet.addToParticles((*inParticles)[index]);
      jet.addToParticles(jetInput);
    }

    edmJets->push_back(jet);
  }


  return StatusCode::SUCCESS;
}


StatusCode TruthJetFCCee::finalize() { return GaudiAlgorithm::finalize(); }


