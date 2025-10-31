#ifndef RECCALORIMETER_TUBELAYERPHIETACALOTOOL_H
#define RECCALORIMETER_TUBELAYERPHIETACALOTOOL_H

#include "RecCaloCommon/CalorimeterToolBase.h"

/** @class TubeLayerPhiEtaCaloTool Reconstruction/RecCalorimeter/src/components/TubeLayerPhiEtaCaloTool.h
 *TubeLayerPhiEtaCaloTool.h
 *
 *  Tool for geometry-dependent settings of the digitisation.
 *  It assumes cylindrical geometry (layers) and phi-eta segmentation.
 *
 *  @author Anna Zaborowska
 */

class TubeLayerPhiEtaCaloTool : public CalorimeterToolBase {
public:
  using CalorimeterToolBase::CalorimeterToolBase;
  virtual ~TubeLayerPhiEtaCaloTool() = default;


protected:
  /** Fill vector with all existing cells for this geometry.
   */
  virtual StatusCode collectCells(std::vector<uint64_t>& cells) const final;


private:
  /// Name of active volumes
  Gaudi::Property<std::string> m_activeVolumeName{this, "activeVolumeName", "LAr_sensitive"};
  /// Name of active layers for sampling calorimeter
  Gaudi::Property<std::string> m_activeFieldName{this, "activeFieldName", "active_layer"};
  /// Name of the fields describing the segmented volume
  Gaudi::Property<std::vector<std::string>> m_fieldNames{this, "fieldNames"};
  /// Values of the fields describing the segmented volume
  Gaudi::Property<std::vector<int>> m_fieldValues{this, "fieldValues"};
  /// Temporary: for use with MergeLayer tool
  Gaudi::Property<unsigned int> m_activeVolumesNumber{this, "activeVolumesNumber", 0};
};

#endif /* RECCALORIMETER_TUBELAYERPHIETACALOTOOL_H */
