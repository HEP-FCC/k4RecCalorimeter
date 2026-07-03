#ifndef ClusterSeedGrower_h
#define ClusterSeedGrower_h 1

/*
 * ClusterSeedGrower
 *
 * Gaudi MultiTransformer that grows pre-formed cluster seeds topologically
 * into full clusters by expanding outward via segmentation neighbours.
 *
 * Growth strategy
 * ---------------
 *   Layer-by-layer BFS: in each layer all seeds simultaneously claim their
 *   reachable neighbours within Von Neumann distance VNDistSeeded.
 *   A neighbour is eligible for growth only when its energy exceeds
 *   GrowThreshold; hits below this threshold block further propagation.
 *   Hits below HardThreshold are discarded entirely.
 *
 *   Contested hits (claimed by two or more seeds at the same BFS layer) are
 *   resolved with the opening angle distance:
 *       d = 1 - cos(angle)
 *   where angle is the opening angle between the cluster log-weighted
 *   barycenter and the hit position.  The cluster with the smaller distance wins.
 *   Cluster energy / barycenter are updated at the end of each layer (not
 *   mid-layer) so that contest resolution within one layer is reproducible.
 *
 *   After seeded BFS, remaining unassigned hitPool cells are used to form
 *   "unseeded" cluster candidates.  A candidate requires at least
 *   MinUnseededHits connected hits all above UnseededThreshold within VN
 *   distance VNDistUnseeded.  Unseeded clusters are then grown topologically
 *   (E >= GrowThreshold, VN distance VNDistUnseeded); contested hits between
 *   two unseeded clusters cause a full merge of the two clusters.
 *
 * Inputs
 * ------
 *   InputSeeds  (std::vector<const edm4hep::ClusterCollection*>)
 *       Seed clusters with pre-attached hits (e.g. from ClusterSeedMerging).
 *       These hits form the initial BFS frontier.
 *
 *   InputHits   (std::vector<const edm4hep::CalorimeterHitCollection*>)
 *       All digitised calorimeter hit collections to grow into.
 *       Treated as a single unified pool; filtered by FieldStrings/FieldValues.
 *
 * Output
 * ------
 *   OutputClusters  (edm4hep::ClusterCollection)
 *       Grown seeded clusters followed by unseeded clusters.
 *       Position = log-weighted barycenter; energy = total energy of all
 *       attached hits.
 */

#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ClusterCollection.h"

#include "k4FWCore/Transformer.h"
#include "Gaudi/Property.h"

#include "ClusterSeedingBase.h"

#include <tuple>
#include <vector>

struct ClusterSeedGrower final
    : ClusterSeedingBase<
          std::tuple<edm4hep::ClusterCollection>(
              const std::vector<const edm4hep::ClusterCollection*>&,
              const std::vector<const edm4hep::CalorimeterHitCollection*>&)>
{
  ClusterSeedGrower(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode finalize() override { return StatusCode::SUCCESS; }

  std::tuple<edm4hep::ClusterCollection>
  operator()(const std::vector<const edm4hep::ClusterCollection*>& seedColls,
             const std::vector<const edm4hep::CalorimeterHitCollection*>& hitColls) const override;

private:
  // Contest resolution strategy for grow().
  enum class ContestStrategy {
    ResolveByDist,   // closest cluster by opening-angle distance keeps hit; others get nothing
    MergeClusters    // contested hit fully merges the two clusters (union-find)
  };

  // Find all connected components in *pool* where every member has
  // E >= threshold, using VN neighbourhood distance vnDist.
  // Returns one Hitmap per component; components with fewer than
  // minHits members are already filtered out before returning.
  std::vector<ClusterSeedingBase::Hitmap> connectedSeeds(const ClusterSeedingBase::Hitmap& pool,
                                                         unsigned int vnDist,
                                                         float threshold,
                                                         unsigned int minHits) const;

  // Layer-by-layer BFS growth of *seeds* into *pool* hits above *growThreshold*.
  // Seeds are moved in; the fully-grown Hitmaps are returned.
  //   ResolveByDist:  contested hits are awarded to the closest cluster by
  //                   opening-angle distance; states are updated after each layer.
  //   MergeClusters:  contested hits cause a full union-find merge of the
  //                   two competing clusters; output contains one Hitmap
  //                   per surviving (root) cluster.
  std::vector<ClusterSeedingBase::Hitmap>
  grow(ContestStrategy                         strategy,
       const ClusterSeedingBase::Hitmap&       pool,
       const std::vector<ClusterSeedingBase::Hitmap>& seeds,
       unsigned int                            vnDist,
       float                                   growThreshold) const;

  // Attach isolated hits above HardThreshold to the nearest cluster by opening-angle distance
  void attachIsolatedHits(const ClusterSeedingBase::Hitmap& pool,
                          std::vector<ClusterSeedingBase::Hitmap>& seededClusters,
                          std::vector<ClusterSeedingBase::Hitmap>& unseededClusters) const;

  // ---- Cell-ID inclusion (daughter-specific) ----
  Gaudi::Property<std::vector<std::string>> m_fieldStringsToInclude{
      this, "FieldStringsToInclude", {},
      "BitField names used to include hits (e.g. 'cherenkov', 'layer')"};
  Gaudi::Property<std::vector<int>> m_fieldValuesToInclude{
      this, "FieldValuesToInclude", {},
      "Values corresponding to FieldStringsToInclude for hit selection"};

  // ---- Energy thresholds ----
  Gaudi::Property<float> m_growThreshold{
      this, "GrowThreshold", 0.01f,
      "Minimum hit energy [GeV] for BFS propagation. "
      "Hits below this threshold block further topological growth."};
  Gaudi::Property<float> m_hardThreshold{
      this, "HardThreshold", 0.01f,
      "Absolute minimum hit energy [GeV]. Hits below this are never attached "
      "to any cluster (HardThreshold <= GrowThreshold is expected)."};
  Gaudi::Property<float> m_unseededThreshold{
      this, "UnseededThreshold", 0.04f,
      "Minimum hit energy [GeV] required to seed an unseeded cluster candidate. "
      "A group of >= MinUnseededHits connected hits all above this threshold "
      "forms one unseeded cluster seed."};

  // ---- Neighbourhood distances ----
  Gaudi::Property<unsigned int> m_vnDistSeeded{
      this, "VNDistSeeded", 2u,
      "Von Neumann neighbour distance used for BFS growth of seeded clusters. "
      "Distance 1 means only direct face-sharing neighbours; distance 2 allows "
      "one jump between the frontier and the candidate hit."};
  Gaudi::Property<unsigned int> m_vnDistUnseeded{
      this, "VNDistUnseeded", 2u,
      "Von Neumann neighbour distance used both for seeding and BFS growth of "
      "unseeded clusters."};
  Gaudi::Property<unsigned int> m_minUnseededHits{
      this, "MinUnseededHits", 2u,
      "Minimum number of connected above-UnseededThreshold hits required to "
      "form an unseeded cluster seed."};

  // ---- Barycenter / distance parameters ----
  Gaudi::Property<float> m_maxOpeningAngle{
      this, "MaxOpeningAngle", 0.3f,
      "Maximum opening angle [rad] between a cluster barycenter and an unowned hit "
      "for the hit to be considered for attachment in Step 5. "
      "Only clusters within this angle are compared by opening-angle distance."};
  Gaudi::Property<bool> m_attachIsolatedHits{
      this, "AttachIsolatedHits", false,
      "Whether to attach remaining unassigned hitPool cells to the nearest cluster "
      "after BFS growth. Only hits above HardThreshold are eligible for attachment."};
};

#endif // ClusterSeedGrower_h
