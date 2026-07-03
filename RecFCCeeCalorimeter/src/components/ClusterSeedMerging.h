#ifndef ClusterSeedMerging_h
#define ClusterSeedMerging_h 1

/*
 * ClusterSeedMerging
 *
 * Gaudi MultiTransformer that merges cluster seeds from the calorimeter-driven
 * (Type A/B) and track-driven (Type C) seeding algorithms into a single set of
 * merged clusters.
 *
 * Merging criterion
 * -----------------
 *   Two seeds i and j are spatially adjacent when the opening angle between
 *   their position vectors is less than arctan(MergeDistance / |pos_i|):
 *
 *     alpha(i,j) < arctan(MergeDistance / |pos_i|)
 *
 *   This defines a fixed angular cone of half-angle arctan(d/r_i) around
 *   seed i.  Adjacency is propagated transitively via BFS. The constraint
 *   that at most one Type-C seed may belong to any merged group is enforced
 *   during BFS (additional C seeds found while growing a component that
 *   already contains one C seed are deferred to their own component).
 *
 * Outputs
 * -------
 *   OutputMergedSeeds  (edm4hep::ClusterCollection)
 *       One cluster per merged group. The cluster position is set to the
 *       energy-weighted barycenter of the constituent seed hits. All
 *       constituent seed hits are attached.
 *
 * Inputs
 * ------
 *   CaloDrivenSeeds   (vector of edm4hep::ClusterCollection)
 *       Output collections from CaloDrivenClusterSeeding (Type A/B).
 *   TrackDrivenSeeds  (vector of edm4hep::ClusterCollection)
 *       Output collections from TrackDrivenClusterSeeding (Type C).
 */

#include "edm4hep/ClusterCollection.h"

#include "k4FWCore/Transformer.h"
#include "Gaudi/Property.h"

#include "ClusterSeedingBase.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

// helper functions
namespace ClusterSeeding {
  enum class SeedType {
    CaloDrivenA = 0x1,
    CaloDrivenB = 0x2,
    TrackDrivenC = 0x4,
    Unseeded = 0x8
  };

  // ---- Cluster::type bit-field layout ----
  //   bits  0..15 : PDG ID  (NH=130, photon=22, ch=211, sat=0)    -- written by downstream PID algorithms
  //   bits 16..31 : SeedType bitmask                              -- written by seeding components
  inline int encodeType(int pdg, SeedType seed) {
    return (pdg & 0xFFFF) | ((static_cast<int>(seed) & 0xFFFF) << 16);
  }
  inline int decodePdg(int t) {
    return static_cast<int16_t>(t & 0xFFFF);
  }
  inline SeedType decodeSeed(int t) {
    return static_cast<SeedType>((t >> 16) & 0xFFFF);
  }

  // ---- Geometry helpers (free functions) ----
  // Convert (x,y,z) to (theta, phi).  Returns (0,0) if |r| ~ 0.
  inline std::pair<float, float> toThetaPhi(float x, float y, float z) {
    const float r = std::sqrt(x * x + y * y + z * z);
    if (r < std::numeric_limits<float>::epsilon()) return {0.f, 0.f};
    return {std::acos(std::max(-1.f, std::min(1.f, z / r))), std::atan2(y, x)};
  }

  // Phi difference ph1 - ph2 wrapped into [-pi, pi].
  inline float deltaPhi(float ph1, float ph2) {
    constexpr float pi = static_cast<float>(M_PI);
    float d = ph1 - ph2;
    if (d > pi)  d -= 2.f * pi;
    if (d < -pi) d += 2.f * pi;
    return d;
  }

  // Angular distance sqrt(dTheta^2 + dPhi_wrapped^2).
  inline float angularDist(float th1, float ph1, float th2, float ph2) {
    const float dth = th1 - th2;
    const float dph = deltaPhi(ph1, ph2);
    return std::sqrt(dth * dth + dph * dph);
  }
} // namespace ClusterSeeding

struct ClusterSeedMerging final
    : ClusterSeedingBase<
          std::tuple<edm4hep::ClusterCollection>( // merged output
               const std::vector<const edm4hep::ClusterCollection*>&, // calo-driven seeds
               const std::vector<const edm4hep::ClusterCollection*>&)> // track-driven seeds
{
  ClusterSeedMerging(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode finalize() override { return StatusCode::SUCCESS; }

  std::tuple<edm4hep::ClusterCollection>
  operator()(const std::vector<const edm4hep::ClusterCollection*>& caloSeeds,
             const std::vector<const edm4hep::ClusterCollection*>& trackSeeds) const override;

  // ---- Merging criterion (configurable) ----
  Gaudi::Property<float> m_mergeDistance{
      this, "MergeDistance", 44.f,
      "d in arctan(d / |pos_i|): half-angle cone size for adjacency (mm)"};
};

#endif // ClusterSeedMerging_h
