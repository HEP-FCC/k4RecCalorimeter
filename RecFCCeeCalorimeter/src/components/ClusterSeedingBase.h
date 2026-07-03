#ifndef ClusterSeedingBase_h
#define ClusterSeedingBase_h 1

/*
 * ClusterSeedingBase<Signature>
 *
 * CRTP-free templated base class for the calorimeter clustering functionals
 *
 * Each daughter class inherits from
 *   ClusterSeedingBase<ReturnType(InputTypes...)>
 * which itself inherits from
 *   k4FWCore::MultiTransformer<ReturnType(InputTypes...)>
 * so the Gaudi framework sees a fully-specified MultiTransformer with exactly
 * the right signature.  The daughter only needs to declare its own
 * constructor/operator() and any extra algorithm-specific properties or
 * methods; everything shared lives here.
 *
 * Shared protected interface
 * --------------------------
 *   initialize() override
 *       Retrieves GeoSvc, resolves the segmentation from ReadoutName, and
 *       builds the internal m_decoder pointer.  Daughters that need additional
 *       initialisation should call ClusterSeedingBase::initialize() first.
 *
 *   passSelection(cellID)
 *       Returns true if the cell passes all (FieldStringsToFilter/Values) and
 *       all (FieldStringsToInclude/Values) cuts.
 *
 *   vonNeumannNeighbors(cellID, half)
 *       Returns the set of cellIDs that are reachable from cellID within Von
 *       Neumann distance *half* using the segmentation's neighbour finder.
 *
 *   isConnected(cells)
 *       Returns true if all cells in the set form one connected component
 *       under VN-d1 adjacency.
 *
 *   calcBarycenter(hits)  /  calcBarycenter(hitMap)
 *       Overloaded log-weighted barycenter helper returning a ClusterState.
 *       w_i = max(0, W0 + ln(E_i / E_tot)).
 *
 * Shared protected properties (all configurable via Gaudi::Property)
 * ------------------------------------------------------------------
 *   ReadoutName            (std::string)
 *   FieldStringsToFilter   (std::vector<std::string>)
 *   FieldValuesToFilter    (std::vector<int>)
 *   W0                     (float)   log-weight offset
 *
 * Shared protected data members
 * ------------------------------
 *   m_geoSvc       SmartIF<IGeoSvc>
 *   m_segmentation dd4hep::DDSegmentation::Segmentation*
 *   m_decoder      const dd4hep::DDSegmentation::BitFieldCoder*
 */

#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/MutableCluster.h"

#include "Gaudi/Property.h"
#include "k4FWCore/Transformer.h"
#include "k4Interface/IGeoSvc.h"

#include "DD4hep/Detector.h"
#include "DD4hep/Readout.h"
#include "DDSegmentation/Segmentation.h"

#include <cmath>
#include <map>
#include <queue>
#include <set>
#include <unordered_map>
#include <vector>

// forward declarations
namespace dd4hep {
namespace DDSegmentation {
  class Segmentation;
  class BitFieldCoder;
} // namespace DDSegmentation
} // namespace dd4hep

// ============================================================
//  ClusterSeedingBase
// ============================================================

template <typename Signature>
struct ClusterSeedingBase : k4FWCore::MultiTransformer<Signature> {

  // ---- Type aliases ----
  // Convenience alias for a map of cellID -> CalorimeterHit, used as the
  // standard hit-pool type throughout the clustering chain.
  using Hitmap = std::unordered_map<uint64_t, edm4hep::CalorimeterHit>;

  // ---- Lightweight cluster state used by barycenter ----
  struct ClusterState {
    float x{0.f}, y{0.f}, z{0.f}, energy{0.f};
  };

protected:
  // ---- Constructor forwarding to MultiTransformer ----
  // Daughters call:
  //   ClusterSeedingBase(name, svcLoc, inputs, outputs)
  // which forwards directly to MultiTransformer.
  using Base = k4FWCore::MultiTransformer<Signature>;
  using Base::Base;

  // ---- Shared initialisation ----
  // Retrieves GeoSvc, resolves segmentation from m_readoutName, and
  // caches m_segmentation / m_decoder.  Daughters that override initialize()
  // should call this first:
  //   if (ClusterSeedingBase::initialize().isFailure()) return StatusCode::FAILURE;
  StatusCode initialize() override {
    m_geoSvc = this->service("GeoSvc");
    if (!m_geoSvc) {
      this->error() << "Unable to locate GeoSvc." << endmsg;
      return StatusCode::FAILURE;
    }

    if (m_geoSvc->getDetector()->readouts().find(m_readoutName) == m_geoSvc->getDetector()->readouts().end()) {
      this->error() << "Readout <<" << m_readoutName << ">> does not exist." << endmsg;
      return StatusCode::FAILURE;
    }

    m_segmentation = m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation();
    m_decoder = m_geoSvc->getDetector()->readout(m_readoutName).segmentation().decoder();

    if (m_fieldStringsToFilter.size() != m_fieldValuesToFilter.size()) {
      this->error() << "FieldStringsToFilter and FieldValuesToFilter must have the same number of entries "
                    << "(got " << m_fieldStringsToFilter.size() << " vs " << m_fieldValuesToFilter.size() << ")."
                    << endmsg;
      return StatusCode::FAILURE;
    }

    return StatusCode::SUCCESS;
  } // initialize

  // ---- Cell-ID selection ----
  // Returns false if:
  //   for each (field, value) in FieldStringsToFilter:  decoder[field] != value
  //   (i.e. the cell is filtered out)
  bool passSelection(uint64_t cellID) const {
    if (!m_decoder)
      return true;

    for (size_t i = 0; i < m_fieldStringsToFilter.size(); ++i) {
      const int64_t value = m_decoder->get(cellID, m_fieldStringsToFilter[i]);
      if (value != static_cast<int64_t>(m_fieldValuesToFilter[i]))
        return false;
    }

    return true;
  } // passSelection

  // ---- Von Neumann neighbourhood ----
  // Returns all cellIDs reachable from cellID within VN distance *half*
  // using the segmentation's neighbour finder (BFS, *half* iterations).
  std::set<uint64_t> vonNeumannNeighbors(uint64_t cellID, unsigned int half) const {
    std::set<uint64_t> frontier{cellID};
    std::set<uint64_t> visited{cellID};
    for (unsigned int d = 0; d < half; ++d) {
      std::set<uint64_t> nextFrontier;
      for (uint64_t cid : frontier) {
        std::set<dd4hep::DDSegmentation::CellID> temp;
        m_segmentation->neighbours(cid, temp);

        // filter out neighbors that fail passSelection
        std::set<dd4hep::DDSegmentation::CellID> nbrs;
        for (const auto& nb : temp) {
          if (passSelection(nb))
            nbrs.insert(nb);
        }

        // form next frontier with unvisited neighbors
        for (auto nb : nbrs) {
          if (visited.insert(nb).second)
            nextFrontier.insert(nb);
        }
      }

      frontier = std::move(nextFrontier);
    } // for each layer

    return visited;
  } // vonNeumannNeighbors

  // ---- Connectivity check (VN-d1) ----
  // Returns true if all cellIDs in *cells* form a single connected component
  // under VN distance-1 adjacency (shared-face neighbours only).
  bool isConnected(const std::set<uint64_t>& cells) const {
    if (cells.size() <= 1)
      return true;

    const uint64_t start = *cells.begin();
    std::set<uint64_t> visited{start};
    std::queue<uint64_t> q;
    q.push(start);

    while (!q.empty()) {
      const uint64_t current = q.front();
      q.pop();
      std::set<uint64_t> nbrs;
      m_segmentation->neighbours(current, nbrs);

      std::set<uint64_t> filteredNbrs;
      for (const auto& nb : nbrs) {
        if (passSelection(nb))
          filteredNbrs.insert(nb);
      }

      for (const uint64_t nb : filteredNbrs) {
        if (cells.count(nb) && !visited.count(nb)) {
          visited.insert(nb);
          q.push(nb);
        }
      }
    }

    return visited.size() == cells.size();
  } // isConnected

  // ---- Log-weighted barycenter (vector<hit> overload) ----
  ClusterState calcBarycenter(const std::vector<edm4hep::CalorimeterHit>& hits) const {
    float totE = 0.f;
    for (const auto& h : hits)
      totE += h.getEnergy();

    ClusterState s;
    s.energy = totE;

    // Guard against zero/negative total energy to avoid NaNs
    if (totE <= 0.f)
      return s;

    float totW = 0.f;
    for (const auto& h : hits) {
      const float w = std::max(0.f, m_w0.value() + std::log(h.getEnergy() / totE));
      totW += w;
      const auto& p = h.getPosition();
      s.x += w * p.x;
      s.y += w * p.y;
      s.z += w * p.z;
    }
    if (totW > 0.f) {
      s.x /= totW;
      s.y /= totW;
      s.z /= totW;
    }

    return s;
  } // calcBarycenter

  // ---- Log-weighted barycenter (Hitmap overload) ----
  ClusterState calcBarycenter(const Hitmap& hitMap) const {
    float totE = 0.f;
    for (const auto& [id, h] : hitMap)
      totE += h.getEnergy();

    ClusterState s;
    s.energy = totE;

    // Guard against zero/negative total energy to avoid NaNs
    if (totE <= 0.f)
      return s;

    float totW = 0.f;
    for (const auto& [id, h] : hitMap) {
      const float w = std::max(0.f, m_w0.value() + std::log(h.getEnergy() / totE));
      totW += w;
      const auto& p = h.getPosition();
      s.x += w * p.x;
      s.y += w * p.y;
      s.z += w * p.z;
    }
    if (totW > 0.f) {
      s.x /= totW;
      s.y /= totW;
      s.z /= totW;
    }

    return s;
  } // calcBarycenter

  // ---- Opening-angle distance ----
  // Returns 1 - cos(angle) between the cluster barycenter direction and (x,y,z).
  // Returns std::numeric_limits<float>::max() if either vector is at the origin.
  float openingAngleDist(const ClusterState& cs, float x, float y, float z) const {
    const float rc = std::sqrt(cs.x * cs.x + cs.y * cs.y + cs.z * cs.z);
    const float rh = std::sqrt(x * x + y * y + z * z);
    if (rc < std::numeric_limits<float>::epsilon() || rh < std::numeric_limits<float>::epsilon()) {
      this->error() << "ClusterSeedingBase::openingAngleDist: zero-length vector! "
                    << "Cluster (x,y,z,E) = (" << cs.x << "," << cs.y << "," << cs.z << "," << cs.energy << ") "
                    << "Point (x,y,z) = (" << x << "," << y << "," << z << ")" << endmsg;
      return std::numeric_limits<float>::max();
    }
    float cosAngle = (cs.x * x + cs.y * y + cs.z * z) / (rc * rh);
    cosAngle = std::max(-1.f, std::min(1.f, cosAngle));

    return 1.f - cosAngle;
  } // openingAngleDist

  // ---- Shared protected properties ----

  // Readout / cell-ID decoding
  Gaudi::Property<std::string> m_readoutName{this, "ReadoutName", "",
                                             "Name of the calorimeter readout (used to retrieve the segmentation)"};

  // Exclude cells matching any (field, value) pair
  Gaudi::Property<std::vector<std::string>> m_fieldStringsToFilter{
      this, "FieldStringsToFilter", {}, "BitField names used to filter (exclude) hits, e.g. 'cherenkov', 'layer'"};
  Gaudi::Property<std::vector<int>> m_fieldValuesToFilter{
      this, "FieldValuesToFilter", {}, "Values corresponding to FieldStringsToFilter"};

  // Barycenter / distance parameters
  Gaudi::Property<float> m_w0{this, "W0", 4.6f,
                              "Log-weight offset in log-weighted barycenter: w_i = max(0, W0 + ln(E_i/E_tot))"};

  // ---- Shared protected data members ----
  SmartIF<IGeoSvc> m_geoSvc;
  dd4hep::DDSegmentation::Segmentation* m_segmentation{nullptr};
  const dd4hep::DDSegmentation::BitFieldCoder* m_decoder{nullptr};
};

#endif // ClusterSeedingBase_h
