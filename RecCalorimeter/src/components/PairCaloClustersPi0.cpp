/**
 * @file RecCalorimeter/src/components/PairCaloClustersPi0.cpp
 * @author Zhibo Wu, scott snyder <snyder@bnl.gov>
 * @date Rewritten Sep, 2026
 * @brief Make pi0 candidates from cluster pairs.
 */

#include "PairCaloClustersPi0.h"
// Key4HEP
#include "k4FWCore/GaudiChecks.h"
#include "k4FWCore/MetadataUtils.h"

#include "TLorentzVector.h"
#include "TVector3.h"

#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/maximum_weighted_matching.hpp"

// Include the <cmath> header for sqrt, pow
#include <cmath>
#include <fstream>
#include <ranges>

DECLARE_COMPONENT(PairCaloClustersPi0)

namespace {

/// Graph data structures.  This is the recommended structure for a
/// non-sparse graph not being modified dynamically.
/// The maximum_weighted_matching algorithm requires the edge weight
/// as an internal property.
// clang-format off
using EdgeProps = boost::property<boost::edge_weight_t, double,
                                  boost::property<boost::edge_index_t, int> >;
using Graph = boost::adjacency_list<boost::vecS,
                                    boost::vecS,
                                    boost::undirectedS,
                                    boost::no_property,
                                    EdgeProps>;
using Vertex = boost::graph_traits<Graph>::vertex_descriptor;
using Edge = boost::graph_traits<Graph>::edge_descriptor;
// clang-format on

/// Helper: Turn a pair of iterators into a range.
template <class IT>
auto make_range(const std::pair<IT, IT>& p) {
  return std::ranges::subrange(p.first, p.second);
}

/// Helper: Extract a 4-vector from a cluster (assuming it comes
/// from the origin).
TLorentzVector getTLV(const edm4hep::Cluster& cl) {
  double e = cl.getEnergy();
  TVector3 disp(cl.getPosition().x, cl.getPosition().y, cl.getPosition().z);
  return TLorentzVector(disp * (e / disp.Mag()), e);
}

/**
 * @brief Create graph corresponding to a set of clusters.
 * @param inClusters The input set of clusters.
 * @param paris Filled with the 4-vector for each candidate pair
 *              (corresponding to an edge in the graph).
 * @param minClusterEnergy Minimum energy for a cluster to be considered.
 * @param maxDTheta Maximum opening angle for a cluster pair to be considered.
 * @param masspeak The pi0 mass.
 * @param masslow Low end of mass window to acccept.
 * @param masshigh Upper end of mass window to acccept.
 */
// clang-format off
Graph makeGraph(const edm4hep::ClusterCollection& inClusters,
                std::vector<TLorentzVector>& pairs,
                double minClusterEnergy,
                double maxDTheta,
                double masspeak,
                double masslow,
                double masshigh)
// clang-format on
{
  // New graph; number of vertices corresponds to the number of clusters.
  Graph g(inClusters.size());

  // Sum of weights.
  double wsum = 0;

  // Loop over cluster pairs.
  for (size_t i = 0; i < inClusters.size(); ++i) {
    // Apply energy requirement to first cluster and get 4-momentum.
    const auto& cl_i = inClusters.at(i);
    if (cl_i.getEnergy() < minClusterEnergy)
      continue;
    TLorentzVector tlv_i = getTLV(cl_i);

    for (size_t j = i + 1; j < inClusters.size(); j++) {
      // Apply energy requirement to second cluster and get 4-momentum.
      const auto& cl_j = inClusters.at(j);
      if (cl_j.getEnergy() < minClusterEnergy)
        continue;
      TLorentzVector tlv_j = getTLV(cl_j);

      // Calculate total 4-momentum; apply mass and opening angle requirements.
      TLorentzVector vpair = tlv_i + tlv_j;
      double invM = vpair.M();
      if (invM > masslow && invM < masshigh && std::abs(tlv_i.Angle(tlv_j.Vect())) < maxDTheta) {
        // Good candidate.  Set the weight to the square difference from the
        // pi0 mass (we'll adjust it later).
        double w = std::pow(invM - masspeak, 2);
        wsum += w;

        // Add the edge to the graph, and remember the 4-momentum.
        boost::add_edge(i, j, EdgeProps(w, pairs.size()), g);
        pairs.push_back(vpair);
      }
    }
  }

  // Adjust weights such that finding the maximum weight will actually find
  // the maximum cardinality, minimum weight solution of the orginal weights.
  for (auto e : make_range(boost::edges(g))) {
    double w = boost::get(boost::edge_weight, g, e);
    boost::put(boost::edge_weight, g, e, 2 * wsum - w);
  }

  return g;
}

/// Find the edeges of the graph corresponding to the desired solution.
std::vector<Edge> findEdges(const Graph& g) {
  // Run the algorithm.
  size_t nv = boost::num_vertices(g);
  std::vector<Vertex> mate(nv);
  boost::maximum_weighted_matching(g, mate.data());

  // Read out the edges of the solution.
  std::vector<Edge> out;
  for (Vertex v1 = 0; v1 < nv; ++v1) {
    Vertex v2 = mate[v1];
    if (v2 != boost::graph_traits<Graph>::null_vertex() && v1 < v2) {
      out.push_back(boost::edge(v1, v2, g).first);
    }
  }

  return out;
}

} // anonymous namespace

PairCaloClustersPi0::PairCaloClustersPi0(const std::string& name, ISvcLocator* svcLoc)
    : Gaudi::Algorithm(name, svcLoc) {
  declareProperty("inClusters", m_inClusters, "Input cluster collection");
  declareProperty("reconstructedPi0", m_reconstructedPi0, "Output1: Reconstructed pi0 collection");
  declareProperty("unpairedClusters", m_unpairedClusters, "Output2: Unpaired cluster collection");
  declareProperty("pairedClusters", m_pairedClusters, "Output3: Paired cluster collection");
}

StatusCode PairCaloClustersPi0::initialize() {
  K4_GAUDI_CHECK(Gaudi::Algorithm::initialize());

  // If there are shapeParameters metadata in the input cluster collection, ship them to the output cluster collections
  auto shapeParameterNames = k4FWCore::getCollectionParameter<std::vector<std::string>>(
                                 m_inClusters.objKey(), edm4hep::labels::ShapeParameterNames, this)
                                 .value_or(std::vector<std::string>{});
  if (shapeParameterNames.size() > 0) {
    k4FWCore::putCollectionParameter(m_pairedClusters.objKey(), edm4hep::labels::ShapeParameterNames,
                                     shapeParameterNames, this);
    k4FWCore::putCollectionParameter(m_unpairedClusters.objKey(), edm4hep::labels::ShapeParameterNames,
                                     shapeParameterNames, this);
  }
  // print pi0 mass window
  info() << "pi0 mass window for cluster pairing: [" << m_massLow << "," << m_massHigh << "] GeV, peak= " << m_massPeak
         << " GeV" << endmsg;

  return StatusCode::SUCCESS;
}

StatusCode PairCaloClustersPi0::execute(const EventContext&) const {
  verbose() << "-------------------------------------------" << endmsg;

  // Get the input collection with clusters
  const edm4hep::ClusterCollection* inClusters = m_inClusters.get();

  // Initialize output clusters
  edm4hep::ReconstructedParticleCollection* reconstructedPi0 = m_reconstructedPi0.createAndPut();
  edm4hep::ClusterCollection* unpairedClusters = m_unpairedClusters.createAndPut();
  edm4hep::ClusterCollection* pairedClusters = m_pairedClusters.createAndPut();

  // clang-format off
  K4_GAUDI_CHECK( doPairing (*inClusters,
                             *reconstructedPi0,
                             *pairedClusters,
                             *unpairedClusters) );
  // clang-format on

  return StatusCode::SUCCESS;
}

/// Cluster pairing
StatusCode PairCaloClustersPi0::doPairing(const edm4hep::ClusterCollection& inClusters,
                                          edm4hep::ReconstructedParticleCollection& reconstructedPi0s,
                                          edm4hep::ClusterCollection& pairedClusters,
                                          edm4hep::ClusterCollection& unpairedClusters) const {
  // Make the graph and find the matching.
  size_t nclust = inClusters.size();
  std::vector<TLorentzVector> pairs;
  // clang-format off
  Graph g = makeGraph (inClusters,
                       pairs,
                       m_minClusterEnergy,
                       m_maxDTheta,
                       m_massPeak,
                       m_massLow,
                       m_massHigh);
  // clang-format on
  std::vector<Edge> edges = findEdges(g);

  // Sort in order of descending energy.
  std::ranges::sort(edges, [&](const Edge& e1, const Edge& e2) {
    int iedge1 = boost::get(boost::edge_index, g, e1);
    int iedge2 = boost::get(boost::edge_index, g, e2);
    return pairs[iedge1].E() > pairs[iedge2].E();
  });

  // Make the pi0 candidates.  We keep the cluster->pair map in used_clusts.
  std::vector<int> used_clusts(nclust, -1);
  for (const Edge& e : edges) {
    int iedge = boost::get(boost::edge_index, g, e);
    const TLorentzVector tlv_pi = pairs[iedge];
    // clang-format off
    edm4hep::MutableReconstructedParticle this_pi0
      (111, tlv_pi.E(),
       edm4hep::Vector3f(tlv_pi.Px(), tlv_pi.Py(), tlv_pi.Pz()),
       edm4hep::Vector3f(0, 0, 0), 0., tlv_pi.M(), 0.,
       edm4hep::CovMatrix4f());
    // clang-format on
    used_clusts[boost::source(e, g)] = reconstructedPi0s.size();
    used_clusts[boost::target(e, g)] = reconstructedPi0s.size();
    reconstructedPi0s.push_back(this_pi0);
  }

  // Go through each input cluster.  If it was used, add it to pairedClusters,
  // and also associate with the pair.  Otherwise, add it to unpairedClusters.
  for (size_t i = 0; i < nclust; ++i) {
    auto cl = inClusters.at(i).clone();
    int ipair = used_clusts[i];
    if (ipair >= 0) {
      reconstructedPi0s[ipair].addToClusters(cl);
      pairedClusters.push_back(cl);
    } else {
      unpairedClusters.push_back(cl);
    }
  }

  return StatusCode::SUCCESS;
}
