/**
 * @file RecCalorimeter/src/components/PairCaloClustersPi0.h
 * @author Zhibo Wu, scott snyder <snyder@bnl.gov>
 * @date Rewritten Sep, 2026
 * @brief Make pi0 candidates from cluster pairs.
 */

#ifndef RECCALORIMETER_PAIRCALOCLUSTERSPI0_H
#define RECCALORIMETER_PAIRCALOCLUSTERSPI0_H

// Key4HEP
#include "k4FWCore/DataHandle.h"

// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolHandle.h"

// our edm
#include "edm4hep/Cluster.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/ReconstructedParticle.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/Vector3d.h"

// class IRndmGenSvc;

// EDM4HEP
namespace edm4hep {
class Cluster;
class MutableCluster;
class ClusterCollection;
class ReconstructedParticle;
class MutableReconstructedParticle;
class ReconstructedParticleCollection;
class Vector3d;
} // namespace edm4hep

/** @class PairCaloClustersPi0Zhibo Wu
 *
 *  Make pi0 candidate (reconstructed particle) from cluster pairs, according to the definition of a pi0 mass window.
 *
 * We find all pairs of clusters with invariant mass within a mass window,
 * such that no cluster is used in more than one pair.  When doing this,
 * we keep as many pairs as possible.  If there are ties, we keep the
 * combination that minimizes the sum of square differences of each
 * pair mass from the pi0 mass.
 *
 * This combinatoric problem can be written as a graph problem.
 * Each cluster corresponds to a vertex and each candidate pair to an edge.
 * Each edge has a weight given by the square of the difference between
 * the pair invariant mass and the pi0 mass.
 * We then want to remove the minimum number of edges such that no vertex
 * has more than one edge; in case of ties, we choose then configuration that
 * minimizes the sumb of weights.
 *
 * This is an example of what is called a `matching' problem and there are
 * standard algorithms to solve it that run in polynomial time
 * and linear space.  Here we use the maximum_weighted_matching algorithm
 * from boost.graph.  This finds the matching (a configuration with no
 * more than one edge for any vertex) with the maximum sum of edge weights.
 * This sounds like it's not really what we want, but if we transform
 * our weights according to wi' -> C - wi, where C is larger than any wi,
 * than maximizing the sum of wi' will give the matching with maximum
 * cardinality with the minimum sum of weights.
 *
 * The algorithm implemented by boost.graph (Galil, https://doi.org/10.1145/6462.6502)
 * has N^3 complexity.
 * The best known is ~ N^2 log N, but the boost.graph version
 * seems to be good enough.
 *
 * Outputs:
 *
 *  reconstructedPi0: A list of reconstructed particles, with energy, momentum, and pointers to a pair of clusters
 *  unpairedClusters: The rest of clusters not involved in the reconstruction of pi0 candidate through the pairing.
 *  pairedClusters: Clusters used in the reconstruction of pi0 candidate.
 *
 *  @author Zhibo Wu
 */

class PairCaloClustersPi0 : public Gaudi::Algorithm {

public:
  PairCaloClustersPi0(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize();

  StatusCode execute(const EventContext&) const;

private:
  /**
   * Cluster pairing algorithm
   *
   * @param[in] inClusters  Pointer to the input cluster collection.
   * @param[out] reconstructedPi0s Container for reconstructed pi0s.
   * @param[out] pairedClusters Container for clusters used for a pi0.
   * @param[out] unpairedClusters Container for clusters not used for a pi0.
   */
  StatusCode doPairing(const edm4hep::ClusterCollection& inClusters,
                       edm4hep::ReconstructedParticleCollection& reconstructedPi0s,
                       edm4hep::ClusterCollection& pairedClusters, edm4hep::ClusterCollection& unpairedClusters) const;

  /// Handle for input calorimeter clusters collection
  mutable k4FWCore::DataHandle<edm4hep::ClusterCollection> m_inClusters{"inClusters", Gaudi::DataHandle::Reader, this};
  /// Handle for reconstructed pi0 particles (output1) collection
  mutable k4FWCore::DataHandle<edm4hep::ReconstructedParticleCollection> m_reconstructedPi0{
      "reconstructedPi0", Gaudi::DataHandle::Writer, this};
  /// Handle for unpaired (output2) calorimeter clusters collection
  mutable k4FWCore::DataHandle<edm4hep::ClusterCollection> m_unpairedClusters{"unpairedClusters",
                                                                              Gaudi::DataHandle::Writer, this};
  /// Handle for paired (output3) calorimeter clusters collection
  mutable k4FWCore::DataHandle<edm4hep::ClusterCollection> m_pairedClusters{"pairedClusters", Gaudi::DataHandle::Writer,
                                                                            this};

  // pi0 mass window
  Gaudi::Property<double> m_massPeak{this, "massPeak", 0.135, "pi0 mass peak [GeV]"};
  Gaudi::Property<double> m_massLow{this, "massLow", 0.0, "lower boundary of pi0 mass window [GeV]"};
  Gaudi::Property<double> m_massHigh{this, "massHigh", 0.27, "upper boundary of pi0 mass window [GeV]"};

  Gaudi::Property<double> m_minClusterEnergy{this, "minClusterEnergy", 0.0, "minimum cluster energy [GeV]"};
  Gaudi::Property<double> m_maxDTheta{this, "maxDTheta", 999, "maximum opening angle of a pair"};
};

#endif /* RECCALORIMETER_PAIRCALOCLUSTERSPI0_H */
