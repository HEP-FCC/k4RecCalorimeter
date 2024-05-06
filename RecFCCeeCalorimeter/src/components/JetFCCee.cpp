#include "JetFCCee.h"

// std
#include <vector>
#include <math.h>

// EDM4hep
#include "edm4hep/ClusterCollection.h"

DECLARE_COMPONENT(JetFCCee)

JetFCCee::JetFCCee(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
  declareProperty("TopoClusterInput", m_inputClusters, "Handle for input cluster collection");
  declareProperty("JetOutput", m_jetCollection, "Handle for output jet collection");
  declareProperty("JetAlg", m_jetAlg, "Name of jet clustering algorithm");
  declareProperty("JetRadius", m_jetRadius, "Jet clustering radius");
  declareProperty("MinPt", m_minPt, "Minimum pT for saved jets");
}


StatusCode JetFCCee::initialize() {
  if (GaudiAlgorithm::initialize().isFailure()) return StatusCode::FAILURE;

  return StatusCode::SUCCESS;
}

StatusCode JetFCCee::execute() {
  // Create output collections
  auto* edmJets = m_jetCollection.createAndPut();


  // TODO: Make it possible to use multiple input types (PFCs, clusters, etc)
  // Retrieve input clusters
  const edm4hep::ClusterCollection *inClusters = m_inputClusters.get();

  fastjet::JetDefinition* jetDef = new fastjet::JetDefinition( m_jetAlgMap[m_jetAlg], m_jetRadius);

  std::vector<fastjet::PseudoJet> clustersPJ;
  int i=0;
  for(auto cluster: *inClusters){
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
      jet.addToClusters((*inClusters)[index]);
    }
    edmJets->push_back(jet);
  }

  return StatusCode::SUCCESS;
}


StatusCode JetFCCee::finalize() { 
  std::cout << m_jetRadius << "\t" << m_minPt << std::endl;
  return GaudiAlgorithm::finalize(); 
}


