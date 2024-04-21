#ifndef RECFCCEECALORIMETER_AUGMENTCLUSTERSFCCEE_H
#define RECFCCEECALORIMETER_AUGMENTCLUSTERSFCCEE_H

// Key4HEP
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/MetaDataHandle.h"
class IGeoSvc;

// Gaudi
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"

// EDM4HEP
namespace edm4hep
{
  class ClusterCollection;
}

namespace dd4hep
{
  namespace DDSegmentation
  {
    class BitFieldCoder;
  }
}

/** @class AugmentClustersFCCee
 *
 *  Add to the cluster shape parameters the sum of the cluster cells energy and barycenter theta/phi coordinates
 *  per layer. The theta position is calculated with a log(E) weighting
 *
 *  @author Alexis Maloizel
 *  @author Giovanni Marchiori
 */

class AugmentClustersFCCee : public Gaudi::Algorithm
{

public:
  AugmentClustersFCCee(const std::string &name, ISvcLocator *svcLoc);
  StatusCode initialize();
  StatusCode execute(const EventContext &evtCtx) const;
  StatusCode finalize();

private:
  /// Handle for input clusters
  mutable DataHandle<edm4hep::ClusterCollection> m_inClusters{"inClusters", Gaudi::DataHandle::Reader, this};
  /// Handle for output clusters
  mutable DataHandle<edm4hep::ClusterCollection> m_outClusters{"outClusters", Gaudi::DataHandle::Writer, this};
  /// Handle for the cluster shape metadata to read and to write
  MetaDataHandle<std::vector<std::string>> m_inShapeParameterHandle{
    m_inClusters,
    edm4hep::shapeParameterNames,
    Gaudi::DataHandle::Reader};
  MetaDataHandle<std::vector<std::string>> m_showerShapeHandle{
    m_outClusters,
    edm4hep::shapeParameterNames,
    Gaudi::DataHandle::Writer};

  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;
  /// ID of the detectors
  Gaudi::Property<std::vector<int>> m_systemIDs{
      this, "systemIDs", {4}, "IDs of systems"};
  /// Name of the detectors (for the metadata)  
  Gaudi::Property<std::vector<std::string>> m_detectorNames{
      this, "systemNames", {"EMB"}, "Names of the detectors, corresponding to systemIDs"};
  /// Numbers of layers of the detectors
  Gaudi::Property<std::vector<size_t>> m_numLayers{
      this, "numLayers", {12}, "Numbers of layers of the systems"};
  /// Weights for each detector layer for theta position log-weighting
  Gaudi::Property<std::vector<std::vector<double>>> m_thetaRecalcLayerWeights{
      this,
      "thetaRecalcWeights",
      {{-1, 3.0, 3.0, 3.0, 4.25, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0}},
      "Weights for each detector layer for theta position log-weighting. If negative use linear weight."};
  /// Name of the detector readouts, corresponding to system IDs in m_systemIDs
  Gaudi::Property<std::vector<std::string>> m_readoutNames{
      this, "readoutNames", {"ECalBarrelModuleThetaMerged"}, "Names of the detector readouts, corresponding to systemIDs"};
  /// Name of the layer/cell field
  Gaudi::Property<std::vector<std::string>> m_layerFieldNames{
      this, "layerFieldNames", {"layer"}, "Identifiers of layers, corresponding to systemIDs"};
};

#endif /* RECFCCEECALORIMETER_AUGMENTCLUSTERSFCCEE_H */
