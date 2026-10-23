#ifndef RECFCCEECALORIMETER_AUGMENTCLUSTERSFCCEE_H
#define RECFCCEECALORIMETER_AUGMENTCLUSTERSFCCEE_H

// Key4HEP
#include "k4FWCore/DataHandle.h"
class IGeoSvc;

// Gaudi
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"

// DD4HEP
// #include "DDSegmentation/Segmentation.h"

// EDM4HEP
namespace edm4hep {
class ClusterCollection;
}

// DD4HEP
namespace dd4hep {
namespace DDSegmentation {
  class BitFieldCoder;
  class Segmentation;
} // namespace DDSegmentation
} // namespace dd4hep

/** @class AugmentClustersFCCee
 *
 *  Add to the cluster shape parameters the sum of the cluster cells energy and barycenter theta/phi coordinates
 *  per layer. The theta position is calculated with a log(E) weighting
 *
 *  Optionally (thetaLayerCorrectionFile), the theta barycentres of the ECAL barrel layers are corrected for the
 *  S-curve (bias towards the cell centre due to the finite granularity), with a parametrisation read from a JSON
 *  file; the corrected values are stored in theta_<detector>_layer_<i> in place of the uncorrected ones.
 *
 *  @author Alexis Maloizel
 *  @author Giovanni Marchiori
 *  @author Tong Li
 */

class AugmentClustersFCCee : public Gaudi::Algorithm {

public:
  AugmentClustersFCCee(const std::string& name, ISvcLocator* svcLoc);
  StatusCode initialize();
  StatusCode execute(const EventContext& evtCtx) const;
  StatusCode finalize();

private:
  /// Handle for input clusters
  mutable k4FWCore::DataHandle<edm4hep::ClusterCollection> m_inClusters{"inClusters", Gaudi::DataHandle::Reader, this};
  /// Handle for output clusters
  mutable k4FWCore::DataHandle<edm4hep::ClusterCollection> m_outClusters{"outClusters", Gaudi::DataHandle::Writer,
                                                                         this};

  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;
  /// ID of the detectors
  Gaudi::Property<std::vector<int>> m_systemIDs{this, "systemIDs", {4}, "IDs of systems"};
  /// Name of the detectors (for the metadata)
  Gaudi::Property<std::vector<std::string>> m_detectorNames{
      this, "systemNames", {"EMB"}, "Names of the detectors, corresponding to systemIDs"};
  /// systemID of EMB, should be retrieved from mydetector.constantAsDouble("DetID_ECAL_Barrel")
  int systemID_EMB;
  /// Numbers of layers of the detectors
  Gaudi::Property<std::vector<size_t>> m_numLayers{this, "numLayers", {11}, "Numbers of layers of the systems"};
  /// Weights for each detector layer for theta position log-weighting
  Gaudi::Property<std::vector<std::vector<double>>> m_thetaRecalcLayerWeights{
      this,
      "thetaRecalcWeights",
      {{-1, 3.0, 3.0, 3.0, 4.25, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0}},
      "Weights for each detector layer for theta position log-weighting. If negative use linear weight."};
  /// Name of the detector readouts, corresponding to system IDs in m_systemIDs
  Gaudi::Property<std::vector<std::string>> m_readoutNames{
      this,
      "readoutNames",
      {"ECalBarrelModuleThetaMerged"},
      "Names of the detector readouts, corresponding to systemIDs"};
  /// Name of the layer/cell field
  Gaudi::Property<std::vector<std::string>> m_layerFieldNames{
      this, "layerFieldNames", {"layer"}, "Identifiers of layers, corresponding to systemIDs"};
  Gaudi::Property<std::vector<std::string>> m_thetaFieldNames{
      this, "thetaFieldNames", {"theta"}, "Identifiers of theta, corresponding to systemIDs"};
  Gaudi::Property<std::vector<std::string>> m_moduleFieldNames{
      this, "moduleFieldNames", {"module"}, "Identifiers of module, corresponding to systemIDs"};
  Gaudi::Property<bool> m_do_photon_shapeVar{
      this, "do_photon_shapeVar", false, "Calculate shape variables for pi0/photon separation: E_ratio, Delta_E etc."};
  Gaudi::Property<bool> m_do_widthTheta_logE_weights{this, "do_widthTheta_logE_weights", false,
                                                     "Calculate width in theta using logE weights in shape variables"};

  // the number of grouped theta and phi cells
  std::vector<int> nMergedThetaCells;
  std::vector<int> nMergedModules;
  // the number of modules / phi cells
  std::vector<int> nModules;

  /// JSON file with the per-layer S-curve correction of the theta barycentres ("": no correction)
  Gaudi::Property<std::string> m_thetaLayerCorrectionFile{
      this, "thetaLayerCorrectionFile", "",
      "JSON file with the parametrisation of the S-curve of the theta barycentre of each ECAL barrel layer. If set, "
      "theta_<detector>_layer_<i> contains the corrected barycentres. Requires linear weights "
      "(thetaRecalcWeights < 0) for the corrected layers."};

  /// One term of the S-curve parametrisation of a layer:
  /// f(2 pi n nu) * sum_i c_i x_i, with f = sin or cos and x_i = hat_k(log E) [* |cos theta|] or |cos theta| or 1
  struct SCurveTerm {
    bool isSin;
    int harmonic;
    std::vector<int> energyNode;    // node of the hat function of each variable, -1 if none
    std::vector<bool> timesAbsCos;  // whether each variable is multiplied by |cos theta|
    std::vector<double> coefficients;
  };
  /// S-curve terms per ECAL barrel layer (empty: layer not corrected)
  std::vector<std::vector<SCurveTerm>> m_thetaCorrTerms;
  /// log(E/GeV) of the energy nodes of the hat functions
  std::vector<double> m_thetaCorrLogENodes;
  /// theta cell width per ECAL barrel layer and lower edge of the theta cells
  std::vector<double> m_thetaCorrCellWidth;
  double m_thetaCorrCellEdge = 0.;
  /// index in m_readoutNames of the ECAL barrel, -1 if absent or no correction
  int m_thetaCorrSystemIndex = -1;
  /// read the correction file and set up the members above
  StatusCode initThetaLayerCorrection();
  /// S-curve correction (in units of the layer cell width) for a theta barycentre
  double thetaLayerCorrection(unsigned layer, double theta, double energy) const;

  /// Limit of debug printing
  mutable uint debugIter = 0;
  Gaudi::Property<uint> m_maxDebugPrint{this, "maxDebugPrint", 10,
                                        "maximum number of debug/warning messages from execute()"};

  void PrintDebugMessage(MsgStream stream, const std::string& text) const;
};

#endif /* RECFCCEECALORIMETER_AUGMENTCLUSTERSFCCEE_H */
