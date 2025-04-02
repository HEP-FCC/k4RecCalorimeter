#ifndef RECCALORIMETER_PAIRCALOCLUSTERSPI0_H
#define RECCALORIMETER_PAIRCALOCLUSTERSPI0_H

// Key4HEP
#include "k4FWCore/DataHandle.h"

// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
//#include "GaudiKernel/IRndmGenSvc.h"
//#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/ToolHandle.h"

// our edm
#include "edm4hep/Cluster.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/ReconstructedParticle.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/Vector3d.h"

//class IRndmGenSvc;

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

/** @class PairCaloClustersPi0
 *
 *  Make pi0 candidate (reconstructed particle) from cluster pairs, according tothe definition of a pi0 mass window
 *  Output1: A list of reconstructed particles, with energy, momentum, and pointers to a pair of clusters
 *  Output2: The rest of clusters not involved in the reconstructioni of pi0 candidate through the pairing.
 *  
 *  @author Zhibo Wu
 */

class PairCaloClustersPi0 : public Gaudi::Algorithm {

public:
  PairCaloClustersPi0(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize();

  StatusCode execute(const EventContext&) const;

  StatusCode finalize();

private:
  /**
   * Cluster pairing algorithm
   *
   * @param[in] inClusters  Pointer to the input cluster collection.
   * @param[in] masspeak    pi0 mass peak.
   * @param[in] masslow     Lower boundary of the pi0 mass window.
   * @param[in] masshigh    Upper boundary of the pi0 mass window.
   *
   * @return                Pointer to the output cluster collection.
   */
  edm4hep::ClusterCollection* ClusterPairing(const edm4hep::ClusterCollection* inClusters, double masspeak, double masslow, double masshigh) const;

  /**
   * Project the energy of a cluster in the pointing direction of the cluster
   *
   * @param[in]  energy     Energy of the cluster.
   * @param[in]  position   Barycenter of the cluster.
   * @param[in]  origin     Origin of the cluster (inferred from cluster pointing)
   *
   * @return                Effective momentum of cluster.
   */
  edm4hep::Vector3d projectMomentum(double energy, edm4hep::Vector3d position, edm4hep::Vector3d origin) const;

  /**
   * Calculate invariant mass of two input clusters
   *
   * @param[in]  E1          Energy of the first cluster.
   * @param[in]  momentum1   Momentum of the first cluster.
   * @param[in]  E2          Energy of the second cluster.
   * @param[in]  momentum2   Momentum of the second cluster.
   *
   * @return                 Invariant mass.
   */
  double getInvariantMass(double E1, edm4hep::Vector3d momentum1 ,double E2, edm4hep::Vector3d momentum2) const;

  /// Handle for input calorimeter clusters collection
  mutable DataHandle<edm4hep::ClusterCollection> m_inClusters{"inClusters", Gaudi::DataHandle::Reader, this};
  /// Handle for paired (output1) calorimeter clusters collection
  mutable DataHandle<edm4hep::ClusterCollection> m_pairedClusters{"pairedClusters", Gaudi::DataHandle::Writer, this};
  /// Handle for unpaired (output2) calorimeter clusters collection
  mutable DataHandle<edm4hep::ClusterCollection> m_unpairedClusters{"unpairedClusters", Gaudi::DataHandle::Writer, this};

  mutable DataHandle<edm4hep::ReconstructedParticleCollection> m_reconstructedPi0{"reconstructedPi0", Gaudi::DataHandle::Writer, this};
  
  // pi0 mass window
  Gaudi::Property<double> m_massPeak{this, "massPeak", 0.135, "pi0 mass peak [GeV]"};
  Gaudi::Property<double> m_massLow{this, "massLow", 0.0, "lower boundary of pi0 mass window [GeV]"};
  Gaudi::Property<double> m_massHigh{this, "massHigh", 0.27, "upper boundary of pi0 mass window [GeV]"};

  /*
  /// Random Number Service
  SmartIF<IRndmGenSvc> m_randSvc;
  /// Randomly choose between 0 and 1, used in cluster pairing
  Rndm::Numbers m_bit;
  */

};

#endif /* RECCALORIMETER_PAIRCALOCLUSTERSPI0_H */
