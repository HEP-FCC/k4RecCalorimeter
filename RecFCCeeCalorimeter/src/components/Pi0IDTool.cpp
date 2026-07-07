// k4FWCore / k4Interface
#include "k4FWCore/Transformer.h"
#include "k4FWCore/MetadataUtils.h"
#include "k4FWCore/DataHandle.h"

// Gaudi
#include "GaudiKernel/ToolHandle.h"

// EDM4hep
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/Vector3f.h"

// STL
#include <vector>
#include <string>
#include <optional>

/** @struct Pi0IDTool
 *
 * Gaudi MultiTransformer that identifies pi0 particles from a collection of clusters
 * with an associated classification score (e.g., from the TRAPPISTTool). The tool retrieves 
 * the classification score from the shape parameters of each cluster based on a specified
 * shape parameter name, and if the score exceeds a user defined threshold, a new 
 * ReconstructedParticle is created and added to the output collection.
 *
 * input: edm4hep::ClusterCollection
 * output: edm4hep::ReconstructedParticleCollection
 *
 *  @author Jacopo Fanini
 *  @date   2026-07
 *
 */

struct Pi0IDTool final 
    : k4FWCore::MultiTransformer<
        std::tuple<
            edm4hep::ReconstructedParticleCollection
        >(
            const edm4hep::ClusterCollection&
        )
    >
{
public:

    Pi0IDTool(const std::string& name, ISvcLocator* svcLoc)
        : MultiTransformer(
            name,
            svcLoc,
            {   
                KeyValues("inClusters", {"clustersWithScore"})
            },
            {   
                KeyValues("outParticles", {"IdentifiedParticles"})
            }
        )
        {}

    StatusCode initialize() override {
        // Retrieve shape parameter names for input cluster collection
        auto inputKey  = this->inputLocations("inClusters")[0];
        auto shapeParameterNames = k4FWCore::getCollectionParameter<std::vector<std::string>>(
                             inputKey, edm4hep::labels::ShapeParameterNames, this)
                             .value_or(std::vector<std::string>{});
        debug() << "Input cluster has " << shapeParameterNames.size() << " names in metadata" << endmsg;
        // Find and save the index of the shape parameter to be used 
        auto it = std::find(shapeParameterNames.begin(), shapeParameterNames.end(), m_shapeParameterScoreName.value());
        if (it != shapeParameterNames.end()) {
            m_scoreIndex = std::distance(shapeParameterNames.begin(), it);
            debug() << "Feature " << m_shapeParameterScoreName.value() << " found in position " << m_scoreIndex.value() << " of collection metadata" << endmsg;
        } else {
            error() << "Feature " << m_shapeParameterScoreName.value() << " not found, aborting..." << endmsg;
            throw std::runtime_error("Feature " + m_shapeParameterScoreName.value() + " not found in collection metadata");
        }

        return StatusCode::SUCCESS;
    }

    std::tuple <edm4hep::ReconstructedParticleCollection> 
    operator()(const edm4hep::ClusterCollection& clustersWithScore) const override {
        edm4hep::ReconstructedParticleCollection outputParticles;
        for (const auto& cluster : clustersWithScore) {
            // Retrieve the score from the shape parameters and if it exceeds the threshold
            // create and fill a new ReconstructedParticle
            float clusterScore = cluster.getShapeParameters(m_scoreIndex.value());
            debug() << "Cluster score: " << clusterScore << endmsg;
            if (clusterScore > m_threshold) {
                auto particle = outputParticles.create();
                particle.addToClusters(cluster);
                double energy = cluster.getEnergy();
                edm4hep::Vector3f position = edm4hep::Vector3f(cluster.getPosition().x, cluster.getPosition().y, cluster.getPosition().z);
                particle.setPDG(111);
                particle.setCharge(0);
                particle.setGoodnessOfPID(clusterScore);
                particle.setEnergy(energy);
                particle.setMass(0.1349768); // From PDG table
                edm4hep::Vector3f momentum = calculateMomentum(energy, position, edm4hep::Vector3f(0, 0, 0));
                particle.setMomentum(momentum);

            }
        }

        return std::make_tuple(std::move(outputParticles));
    }

    StatusCode finalize() override {
        return StatusCode::SUCCESS;
    }


private:

    std::optional<std::size_t> m_scoreIndex = std::nullopt;

    Gaudi::Property<std::string> m_shapeParameterScoreName{this, "ShapeParameterName", "ClusterScore", 
                                 "Name of the shape parameter to be used for particle identification"};

    Gaudi::Property<float> m_threshold{this, "Threshold", 0.5, 
                           "Threshold for particle identification based on the shape parameter value"};

    edm4hep::Vector3f calculateMomentum(double energy, 
                                        const edm4hep::Vector3f& position,
                                        const edm4hep::Vector3f& origin) const {
      double dirx = position.x - origin.x;
      double diry = position.y - origin.y;
      double dirz = position.z - origin.z;
      double quadsum_dir = sqrt(pow(dirx, 2) + pow(diry, 2) + pow(dirz, 2));
      double px = energy * dirx / quadsum_dir;
      double py = energy * diry / quadsum_dir;
      double pz = energy * dirz / quadsum_dir;
      return edm4hep::Vector3f(px, py, pz);
    }
    
};

DECLARE_COMPONENT(Pi0IDTool)