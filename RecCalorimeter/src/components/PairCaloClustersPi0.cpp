#include "PairCaloClustersPi0.h"

// Include the <cmath> header for sqrt, pow
#include <cmath>

DECLARE_COMPONENT(PairCaloClustersPi0)

PairCaloClustersPi0::PairCaloClustersPi0(const std::string& name, ISvcLocator* svcLoc)
    : Gaudi::Algorithm(name, svcLoc) {
  declareProperty("inClusters", m_inClusters, "Input cluster collection");
  declareProperty("reconstructedPi0", m_reconstructedPi0, "Output1: Reconstructed pi0 collection");
  declareProperty("unpairedClusters", m_unpairedClusters, "Output2: Unpaired cluster collection");
  declareProperty("pairedClusters", m_pairedClusters, "Output3: Paired cluster collection");
}

StatusCode PairCaloClustersPi0::initialize() {
  {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure()) {
      return sc;
    }
  }

  // print pi0 mass window
  info() << "pi0 mass window for cluster pairing: [" <<m_massLow<<","<<m_massHigh<<"] GeV, peak= "<<m_massPeak<<" GeV"<< endmsg;

  return StatusCode::SUCCESS;
}

StatusCode PairCaloClustersPi0::execute(const EventContext&) const {
  verbose() << "-------------------------------------------" << endmsg;

  // Get the input collection with clusters
  const edm4hep::ClusterCollection* inClusters = m_inClusters.get();

  // Initialize output clusters
  edm4hep::ClusterCollection* outClusters = ClusterPairing(inClusters, m_massPeak, m_massLow, m_massHigh);

  if (!outClusters) {
    error() << "Something went wrong in initialization of the output cluster collection, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

StatusCode PairCaloClustersPi0::finalize() { return Gaudi::Algorithm::finalize(); }

// calculation of invariant mass
double PairCaloClustersPi0::getInvariantMass(
    double E1, edm4hep::Vector3d momentum1, double E2, edm4hep::Vector3d momentum2) const{
  double invM = pow(E1+E2, 2) - pow(momentum1.x+momentum2.x, 2) - pow(momentum1.y+momentum2.y, 2) - pow(momentum1.z+momentum2.z, 2);
  if (invM>0){
    invM=sqrt(invM);
  }
  else{
    invM=0;
  }
  return invM;
}

// projection: energy -> 3 momentum of a cluster
edm4hep::Vector3d PairCaloClustersPi0::projectMomentum(double energy, edm4hep::Vector3d position, edm4hep::Vector3d origin) const{
  // the direction of a cluster is defined to be that from the inferred origin of the cluster to the barycenter of the cluster
  double dirx = position.x-origin.x;
  double diry = position.y-origin.y;
  double dirz = position.z-origin.z;
  double quadsum_dir = sqrt( abs( pow(dirx,2) + pow(diry,2) + pow(dirz,2) ) );
  double px = energy * dirx / quadsum_dir;
  double py = energy * diry / quadsum_dir;
  double pz = energy * dirz / quadsum_dir;
  return edm4hep::Vector3d(px, py, pz);
}

// cluster pairing
edm4hep::ClusterCollection*
PairCaloClustersPi0::ClusterPairing(const edm4hep::ClusterCollection* inClusters, double masspeak, double masslow, double masshigh) const {

  edm4hep::ReconstructedParticleCollection* reconstructedPi0 = m_reconstructedPi0.createAndPut();
  edm4hep::ClusterCollection* unpairedClusters = m_unpairedClusters.createAndPut();
  edm4hep::ClusterCollection* pairedClusters = m_pairedClusters.createAndPut();

  // ***** Step 1: Get all possible cluster pairs in the mass window, overlap of clusters allowed *****
  verbose() << "We are in cluster pairing, step 1" <<endmsg;
  std::vector<std::pair<size_t, size_t>> vec_AllPossiblePairs;
  for (size_t i=0; i<inClusters->size(); ++i){
    double energy_i=inClusters->at(i).getEnergy();
    edm4hep::Vector3d cluster_i_position3d(inClusters->at(i).getPosition().x, inClusters->at(i).getPosition().y, inClusters->at(i).getPosition().z);
    // For the moment, the cluster direction uses the pointing assumption: from (0,0,0) to the cluster position. Waiting for the update of cluster direction pointing algorithm.
    edm4hep::Vector3d cluster_i_momentum=PairCaloClustersPi0::projectMomentum(energy_i, cluster_i_position3d, edm4hep::Vector3d(0,0,0));
    for(size_t j=i+1; j<inClusters->size(); j++){
      double energy_j=inClusters->at(j).getEnergy();
      edm4hep::Vector3d cluster_j_position3d(inClusters->at(j).getPosition().x, inClusters->at(j).getPosition().y, inClusters->at(j).getPosition().z);
      // For the moment, the cluster direction uses the pointing assumption: from (0,0,0) to the cluster position. Waiting for the update of cluster direction pointing algorithm.
      edm4hep::Vector3d cluster_j_momentum=PairCaloClustersPi0::projectMomentum(energy_j, cluster_j_position3d, edm4hep::Vector3d(0,0,0));
      double invM=PairCaloClustersPi0::getInvariantMass(energy_i, cluster_i_momentum, energy_j, cluster_j_momentum);
      if (invM>masslow && invM<masshigh){
	std::pair<size_t, size_t> this_possible_pair=std::make_pair(i,j);
	vec_AllPossiblePairs.push_back(this_possible_pair);
      }
    }
  }
  
  // ***** Step 2: make all possible combinations of cluster pairing without overlap *****
  verbose() << "We are in cluster pairing step 2" <<endmsg;
  std::vector<std::vector<std::pair<size_t, size_t>>> vec_combi_pairs;
  // initialisation
  for (size_t i=0; i<vec_AllPossiblePairs.size(); i++){
    std::vector<std::pair<size_t, size_t>> this_combi_pairs;
    this_combi_pairs.push_back(vec_AllPossiblePairs[i]);
    vec_combi_pairs.push_back(this_combi_pairs);
  }
  verbose() << "Initialisation with all possible pairs, size =  " << vec_combi_pairs.size() <<endmsg;
  // iterations
  size_t start_of_iter = 0; // starting index of pair combination in this iteration
  size_t end_of_iter = vec_combi_pairs.size(); // end index of pair combination in this iteration
  // run at most (Number of Pairs - 1) iterations to get all combinations
  for (size_t i_iter=0; i_iter+1<vec_AllPossiblePairs.size(); i_iter++){
    // in each iteration, loop over relevant pair combinations to this iteration
    for (size_t i_combi=start_of_iter; i_combi<end_of_iter; i_combi++){
      // try to extend the list of pair combinations by one more pair
      for (size_t i_pair=0; i_pair<vec_AllPossiblePairs.size(); i_pair++){
	bool IsClusterOverlap=false;
	for (size_t j_pair=0; j_pair<vec_combi_pairs[i_combi].size(); j_pair++){
	  if (vec_AllPossiblePairs[i_pair].first >= vec_combi_pairs[i_combi][j_pair].first || // avoid double counting combinations
	      vec_AllPossiblePairs[i_pair].first == vec_combi_pairs[i_combi][j_pair].second ||
	      vec_AllPossiblePairs[i_pair].second == vec_combi_pairs[i_combi][j_pair].first ||
	      vec_AllPossiblePairs[i_pair].second == vec_combi_pairs[i_combi][j_pair].second){
	    IsClusterOverlap=true;
	    break;
	  }
	}
	if (IsClusterOverlap==false){
	  std::vector<std::pair<size_t, size_t>> vec_combi_copy;
	  for (size_t k_pair=0; k_pair<vec_combi_pairs[i_combi].size(); k_pair++){
	    vec_combi_copy.push_back(vec_combi_pairs[i_combi][k_pair]);
	  }
	  vec_combi_copy.push_back(vec_AllPossiblePairs[i_pair]);
	  vec_combi_pairs.push_back(vec_combi_copy);
	}
      }
    }
    start_of_iter=end_of_iter;
    end_of_iter=vec_combi_pairs.size();
  }
  
  // ***** Step 3: choose the combination with the most number of pairs and closest to the pi0 mass
  verbose() << "We are in cluster pairing step 3" <<endmsg;
  // Step 3(a): look for the maximal number of pairs
  size_t max_Npairs=0;
  std::vector<std::vector<std::pair<size_t, size_t>>> vec_Maxcombi_pairs;
  std::vector<std::pair<size_t, size_t>> bestcombi_pairs;
  for (size_t i_combi=0; i_combi<vec_combi_pairs.size(); i_combi++){
    size_t this_Npairs=vec_combi_pairs[i_combi].size();
    if (this_Npairs>max_Npairs){
      max_Npairs=this_Npairs;
      vec_Maxcombi_pairs.clear();
      vec_Maxcombi_pairs.push_back(vec_combi_pairs[i_combi]);
    }
    else if (this_Npairs==max_Npairs){
      vec_Maxcombi_pairs.push_back(vec_combi_pairs[i_combi]);
    }
  }
  // Step 3(b): if there are multiple combinations with the same maximal number of pairs
  // -> look for the smallest sum of mass deviation with respect to the input pi0 mass peak
  if (vec_Maxcombi_pairs.size()>1){
    size_t index_best_combi=-1;
    double best_devM=1000.;
    for (size_t i_combi=0; i_combi<vec_Maxcombi_pairs.size(); i_combi++){
      double this_devM=0.;
      for (size_t i_pair=0; i_pair<vec_Maxcombi_pairs[i_combi].size(); i_pair++){
	auto this_cluster1=inClusters->at(vec_Maxcombi_pairs[i_combi][i_pair].first);
	auto this_cluster2=inClusters->at(vec_Maxcombi_pairs[i_combi][i_pair].second);
	edm4hep::Vector3d position1(this_cluster1.getPosition().x, this_cluster1.getPosition().y, this_cluster1.getPosition().z);
	edm4hep::Vector3d position2(this_cluster2.getPosition().x, this_cluster2.getPosition().y, this_cluster2.getPosition().z);
	edm4hep::Vector3d momentum1=PairCaloClustersPi0::projectMomentum(this_cluster1.getEnergy(), position1, edm4hep::Vector3d(0,0,0));
	edm4hep::Vector3d momentum2=PairCaloClustersPi0::projectMomentum(this_cluster2.getEnergy(), position2, edm4hep::Vector3d(0,0,0));
	double this_invM=PairCaloClustersPi0::getInvariantMass(
	    this_cluster1.getEnergy(), momentum1, this_cluster2.getEnergy(), momentum2);
	this_devM+=pow((this_invM - masspeak),2);
      }
      this_devM/=vec_Maxcombi_pairs[i_combi].size();
      if (this_devM<best_devM){
	index_best_combi=i_combi;
	best_devM=this_devM;
      }
      // very unlikely, but if the sum of mass deviation is the same between two combinations, we have to randomly choose one
      else if (this_devM==best_devM){
	index_best_combi = (rand()%2) ? i_combi : index_best_combi;
      }
    }
    bestcombi_pairs=vec_Maxcombi_pairs[index_best_combi];
  }
  else if (vec_Maxcombi_pairs.size()==1){
    bestcombi_pairs=vec_Maxcombi_pairs[0];
  }
  
  // Step 4: save paired and unpaired clusters, reconstruct pi0 candidates
  verbose() << "We are in cluster pairing step 4" <<endmsg;
  std::vector<size_t> vec_index_paired_clusters;
  // save paired clusters
  for (size_t i=0; i<bestcombi_pairs.size(); i++){
    auto outCluster1 = inClusters->at(bestcombi_pairs[i].first).clone();
    pairedClusters->push_back(outCluster1);
    vec_index_paired_clusters.push_back(bestcombi_pairs[i].first);
    auto outCluster2 = inClusters->at(bestcombi_pairs[i].second).clone();
    pairedClusters->push_back(outCluster2);
    vec_index_paired_clusters.push_back(bestcombi_pairs[i].second);
    // reconstruct pi0 from these two clusters
    edm4hep::Vector3d position1(outCluster1.getPosition().x, outCluster1.getPosition().y, outCluster1.getPosition().z);
    edm4hep::Vector3d position2(outCluster2.getPosition().x, outCluster2.getPosition().y, outCluster2.getPosition().z);
    edm4hep::Vector3d momentum1=PairCaloClustersPi0::projectMomentum( outCluster1.getEnergy(), position1, edm4hep::Vector3d(0,0,0));
    edm4hep::Vector3d momentum2=PairCaloClustersPi0::projectMomentum( outCluster2.getEnergy(), position2, edm4hep::Vector3d(0,0,0));
    double this_pi0_invM=PairCaloClustersPi0::getInvariantMass( outCluster1.getEnergy(), momentum1, outCluster2.getEnergy(), momentum2 );
    verbose() << "Final pairing " << i << " first cluster = " << bestcombi_pairs[i].first << ", second cluster =  " << bestcombi_pairs[i].second << ", invariant mass [GeV] = " <<  this_pi0_invM <<endmsg;
    
    edm4hep::MutableReconstructedParticle this_pi0(
	111,
	outCluster1.getEnergy()+outCluster2.getEnergy(),
	edm4hep::Vector3f(momentum1.x+momentum2.x, momentum1.y+momentum2.y, momentum1.z+momentum2.z),
	edm4hep::Vector3f(0,0,0),
	0.,
	this_pi0_invM,
	0.,
	edm4hep::CovMatrix4f()
	);
    this_pi0.addToClusters(outCluster1);
    this_pi0.addToClusters(outCluster2);
    reconstructedPi0->push_back(this_pi0);
  }
  for (size_t i=0; i<inClusters->size(); ++i){
    bool IsPaired=false;
    for (size_t j=0; j<vec_index_paired_clusters.size(); j++){
      if (i==vec_index_paired_clusters[j]){
	IsPaired=true;
	break;
      }
    }
    // save unpaired clusters
    if (!IsPaired){
      auto outCluster = inClusters->at(i).clone();
      unpairedClusters->push_back(outCluster);
      verbose() << "save unpaired cluster. cluster index = " << i <<endmsg;
    }
  }

  return unpairedClusters;
}
