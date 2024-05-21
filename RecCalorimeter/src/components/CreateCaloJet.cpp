#include "CreateCaloJet.h"

// std
#include <vector>
#include <math.h>

// EDM4hep
#include "edm4hep/ClusterCollection.h"

DECLARE_COMPONENT(CreateCaloJet)

//CreateCaloJet::CreateCaloJet(const std::string& name, ISvcLocator* svcLoc) : Gaudi::Functional::Transformer(name, svcLoc,
CreateCaloJet::CreateCaloJet(const std::string& name, ISvcLocator* svcLoc) : Transformer(name, svcLoc,
                    KeyValue("InputCollection", "CorrectedCaloClusters"),
                    KeyValue("OutputCollection", "Jets")) {
  declareProperty("JetAlg", m_jetAlg, "Name of jet clustering algorithm");
  declareProperty("JetRadius", m_jetRadius, "Jet clustering radius");
  declareProperty("MinPt", m_minPt, "Minimum pT for saved jets");
}


StatusCode CreateCaloJet::initialize() {
  if (m_jetAlgMap.find(m_jetAlg) == m_jetAlgMap.end()) {
    error() << m_jetAlg << " is not in the list of supported jet algorithms" << endmsg;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

colltype_out  CreateCaloJet::operator()(const colltype_in& input) const{
  edm4hep::ReconstructedParticleCollection edmJets = edm4hep::ReconstructedParticleCollection();


  const fastjet::JetAlgorithm jetAlg = m_jetAlgMap.at(m_jetAlg);
  fastjet::JetDefinition* jetDef = new fastjet::JetDefinition( jetAlg, m_jetRadius);

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
    clusterPJ.set_user_info(new ClusterInfo(i));
    clustersPJ.push_back(clusterPJ);
    i++;
  }


  fastjet::ClusterSequence clustSeq(clustersPJ, *jetDef);
  std::vector <fastjet::PseudoJet> inclusiveJets = fastjet::sorted_by_pt(clustSeq.inclusive_jets(m_minPt));

  //Add a reconstructed particle for each jet
  for(auto cjet : inclusiveJets){
    edm4hep::MutableReconstructedParticle jet;
    jet.setMomentum(edm4hep::Vector3f(cjet.px(), cjet.py(), cjet.pz()));
    jet.setEnergy(cjet.e());
    jet.setMass(cjet.m());

    // Also add the clusters that were used to make the jets
    std::vector<fastjet::PseudoJet> constits = cjet.constituents();
    for(auto constit : constits){
      int index = constit.user_info<ClusterInfo>().index();
      jet.addToClusters((input)[index]);
    }
    edmJets.push_back(jet);
  }

  return edmJets;
 // return StatusCode::SUCCESS;
}


