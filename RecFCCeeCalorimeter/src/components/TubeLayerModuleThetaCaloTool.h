#ifndef RECFCCEECALORIMETER_TUBELAYERMODULETHETAMERGEDCALOTOOL_H
#define RECFCCEECALORIMETER_TUBELAYERMODULETHETAMERGEDCALOTOOL_H

// from Gaudi
#include "GaudiKernel/AlgTool.h"

// k4FWCore
#include "k4Interface/ICalorimeterTool.h"

class IGeoSvc;

/** @class TubeLayerModuleThetaCaloTool k4RecCalorimeter/RecFCCeeCalorimeter/src/components/TubeLayerModuleThetaCaloTool.h
 *  TubeLayerModuleThetaCaloTool.h
 *
 *  Tool for geometry-dependent settings of the digitisation.
 *  It assumes cylindrical geometry (layers) and phi-theta segmentation.
 *
 *  @author Anna Zaborowska
 *  @author Zhibo Wu
 */

class TubeLayerModuleThetaCaloTool : public AlgTool, virtual public ICalorimeterTool {
public:
  TubeLayerModuleThetaCaloTool(const std::string& type, const std::string& name, const IInterface* parent);
  virtual ~TubeLayerModuleThetaCaloTool() = default;
  virtual StatusCode initialize() final;
  virtual StatusCode finalize() final;
  /** Prepare a map of all existing cells in current geometry.
   *   Active layers (cylindrical tubes) are looked in the geometry manager by name ('\b activeVolumeName').
   *   Corresponding bitfield name is given in '\b activeFieldName'.
   *   If users wants to limit the number of active layers, it is possible by setting '\b activeVolumesNumber'.
   *   The total number of cells N = n_layer * n_theta * n_phi, where
   *   n_layer is number of layers (taken from geometry if activeVolumesNumber not set),
   *   n_theta is number of eta bins in that layer,
   *   n_phi is number of phi bins (the same for each layer).
   *   For more explanation please [see reconstruction documentation](@ref md_reconstruction_doc_reccalorimeter).
   *   @param[out] aCells map of existing cells (and deposited energy, set to 0)
   *   return Status code.
   */
  virtual StatusCode prepareEmptyCells(std::unordered_map<uint64_t, double>& aCells) final;

private:
  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;
  /// Name of the detector readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", ""};
  /// Name of active volumes
  Gaudi::Property<std::string> m_activeVolumeName{this, "activeVolumeName", "LAr_sensitive"};
  /// Name of active layers for sampling calorimeter
  Gaudi::Property<std::string> m_activeFieldName{this, "activeFieldName", "active_layer"};
  /// Name of the fields describing the segmented volume
  Gaudi::Property<std::vector<std::string>> m_fieldNames{this, "fieldNames"};
  /// Values of the fields describing the segmented volume
  Gaudi::Property<std::vector<int>> m_fieldValues{this, "fieldValues"};
  /// Number of layers
  Gaudi::Property<unsigned int> m_activeVolumesNumber{this, "activeVolumesNumber", 0};
};

#endif /* RECFCCEECALORIMETER_TUBELAYERMODULETHETAMERGEDCALOTOOL_H */
