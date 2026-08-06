// k4FWCore / k4Interface
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/MetadataUtils.h"
#include "k4FWCore/Transformer.h"

// Gaudi
#include "GaudiKernel/ToolHandle.h"

// EDM4hep
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/Vector3f.h"

// STL
#include <optional>
#include <string>
#include <vector>

/** @struct ClusterPi0PhotonID
 *
 * Gaudi MultiTransformer that identifies photons and pi0 candidates from a
 * collection of clusters with associated classification scores (e.g. produced by the
 * TRAPPISTTool). For each input cluster, the algorithm retrieves the classification score from the cluster
 * shape parameters using a configurable shape parameter name. If the score is greater than
 * a user-defined threshold, the cluster is identified as a pi0 candidate; otherwise,
 * it is identified as a photon candidate. A corresponding ReconstructedParticle is created
 * and added to the corresponding output collection.
 *
 * input: edm4hep::ClusterCollection
 * output: edm4hep::ReconstructedParticleCollection
 *
 *  @author Jacopo Fanini
 *  @date   2026-07
 *
 */

struct ClusterPi0PhotonID final
    : k4FWCore::MultiTransformer<
          std::tuple<edm4hep::ReconstructedParticleCollection, edm4hep::ReconstructedParticleCollection>(
              const edm4hep::ClusterCollection&)> {
public:
  ClusterPi0PhotonID(const std::string& name, ISvcLocator* svcLoc)
      : MultiTransformer(name, svcLoc, {KeyValues("inClusters", {"clustersWithScore"})},
                         {KeyValues("outPi0Particles", {"IdentifiedPi0Particles"}),
                          KeyValues("outPhotonParticles", {"IdentifiedPhotonParticles"})}) {}

  StatusCode initialize() override {
    // Retrieve shape parameter names for input cluster collection
    auto inputKey = this->inputLocations("inClusters")[0];
    auto shapeParameterNames =
        k4FWCore::getCollectionParameter<std::vector<std::string>>(inputKey, edm4hep::labels::ShapeParameterNames, this)
            .value_or(std::vector<std::string>{});
    debug() << "Input cluster has " << shapeParameterNames.size() << " names in metadata" << endmsg;
    // Find and save the index of the shape parameter to be used
    auto it = std::find(shapeParameterNames.begin(), shapeParameterNames.end(), m_shapeParameterScoreName.value());
    if (it != shapeParameterNames.end()) {
      m_scoreIndex = std::distance(shapeParameterNames.begin(), it);
      debug() << "Feature " << m_shapeParameterScoreName.value() << " found in position " << m_scoreIndex.value()
              << " of collection metadata" << endmsg;
    } else {
      throw std::runtime_error("Feature " + m_shapeParameterScoreName.value() + " not found in collection metadata");
    }

    return StatusCode::SUCCESS;
  }

  std::tuple<edm4hep::ReconstructedParticleCollection, edm4hep::ReconstructedParticleCollection>
  operator()(const edm4hep::ClusterCollection& clustersWithScore) const override {
    edm4hep::ReconstructedParticleCollection outputPi0Particles;
    edm4hep::ReconstructedParticleCollection outputPhotonParticles;
    for (const auto& cluster : clustersWithScore) {
      // Retrieve the score from the shape parameters and if it exceeds the threshold
      // create and fill a new ReconstructedParticle
      float clusterScore = cluster.getShapeParameters(m_scoreIndex.value());
      debug() << "Cluster score: " << clusterScore << ", threshold: " << m_threshold << endmsg;
      if (clusterScore > m_threshold) {
        buildParticle(outputPi0Particles, cluster, 111, 0.1349768, clusterScore);
      } else {
        buildParticle(outputPhotonParticles, cluster, 22, 0.0, (1. - clusterScore));
      }
    }
    return std::make_tuple(std::move(outputPi0Particles), std::move(outputPhotonParticles));
  }

  StatusCode finalize() override { return StatusCode::SUCCESS; }

private:
  std::optional<std::size_t> m_scoreIndex = std::nullopt;

  Gaudi::Property<std::string> m_shapeParameterScoreName{
      this, "ShapeParameterName", "ClusterTRAPPISTScore",
      "Name of the shape parameter to be used for particle identification"};

  Gaudi::Property<float> m_threshold{this, "Threshold", 0.5,
                                     "Threshold for particle identification based on the shape parameter value"};

  edm4hep::Vector3f calculateMomentum(double energy, const edm4hep::Vector3f& position,
                                      const edm4hep::Vector3f& origin) const {
    double dirx = position.x - origin.x;
    double diry = position.y - origin.y;
    double dirz = position.z - origin.z;
    double quadsum_dir = std::sqrt(dirx * dirx + diry * diry + dirz * dirz);
    double px = energy * dirx / quadsum_dir;
    double py = energy * diry / quadsum_dir;
    double pz = energy * dirz / quadsum_dir;
    return edm4hep::Vector3f(px, py, pz);
  }

  void buildParticle(edm4hep::ReconstructedParticleCollection& collection, const edm4hep::Cluster& cluster, int pdg,
                     float mass, float pidScore) const {

    auto particle = collection.create();

    particle.addToClusters(cluster);

    const double energy = cluster.getEnergy();
    particle.setEnergy(energy);

    particle.setPDG(pdg);
    particle.setCharge(0);
    particle.setMass(mass);
    particle.setGoodnessOfPID(pidScore);

    edm4hep::Vector3f position{cluster.getPosition().x, cluster.getPosition().y, cluster.getPosition().z};
    // Assuming particle is coming from the IP, will change this once the tool to retrieve cluster direction is
    // implemented
    particle.setMomentum(calculateMomentum(energy, position, edm4hep::Vector3f{0, 0, 0}));
  }
};

DECLARE_COMPONENT(ClusterPi0PhotonID)