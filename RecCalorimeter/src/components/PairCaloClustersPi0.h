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
#include "edm4hep/Vector3d.h"

//class IRndmGenSvc;

// EDM4HEP
namespace edm4hep {
class Cluster;
class MutableCluster;
class ClusterCollection;
class Vector3d;
} // namespace edm4hep

/** @class CorrectCaloClusters
 *
 *  Apply corrections to the clusters reconstructed in ECAL.
 *  * Upstream energy correction corrects for the energy lost by the particles in the material before they enter
 *    active volume of the calorimeter. The correction is parametrized in one (cluster energy) or two variables (cluster
 *    energy, cluster angle)
 *  * Downstream energy correction corrects for the energy lost due to punch trough. This energy is deposited in the
 *    instrumentation behind the active volume of the calorimeter. It is parametrized in one (cluster energy) or two
 *    variables (cluster energy, cluster angle)
 * * Benchmark calibration should be used for the combined simulation of ECal and HCal when using charged pions.
 *    At the input level, the ECal should be calibrated to EM scale and HCal should be calibrated to HAD scale.
 *    The aim of the benchmark calibration is to bring ECal to HAD scale and also to take into account
 *    the energy loss between the ECal and HCal (e.g. in cryostat) - for this, the energy from the last ECal layer and
 * the first HCal layer is used. While, downstream correction is part of the benchmark method (energy lost between ECal
 * and HCal), as well as the upstream correction for ECal For standalone ECal, include upstream and downstream
 * correction; for ECal+HCal simulation apply only benchmark correction should be applied To obtain the actual
 * parameters run RecCalorimeter/tests/options/fcc_ee_caloBenchmarkCalibration.py which calls CalibrateBenchmarkMethod
 *
 *  Based on similar corrections by Jana Faltova and Anna Zaborowska.
 *
 *  @author Juraj Smiesko, benchmark calibration added by Michaela Mlynarikova
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
