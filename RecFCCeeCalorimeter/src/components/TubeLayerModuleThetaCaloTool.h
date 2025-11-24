#ifndef RECFCCEECALORIMETER_TUBELAYERMODULETHETACALOTOOL_H
#define RECFCCEECALORIMETER_TUBELAYERMODULETHETACALOTOOL_H

#include "RecCaloCommon/CalorimeterToolBase.h"

/** @class TubeLayerModuleThetaCaloTool
 * k4RecCalorimeter/RecFCCeeCalorimeter/src/components/TubeLayerModuleThetaCaloTool.h TubeLayerModuleThetaCaloTool.h
 *
 *  Tool for geometry-dependent settings of the digitisation.
 *  It assumes cylindrical geometry (layers) and phi-theta segmentation.
 *
 *  @author Anna Zaborowska
 *  @author Zhibo Wu
 */

class TubeLayerModuleThetaCaloTool : public CalorimeterToolBase {
public:
  using CalorimeterToolBase::CalorimeterToolBase;
  virtual ~TubeLayerModuleThetaCaloTool() = default;


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
  /// Number of layers
  Gaudi::Property<unsigned int> m_activeVolumesNumber{this, "activeVolumesNumber", 0};
};

#endif /* RECFCCEECALORIMETER_TUBELAYERMODULETHETACALOTOOL_H */
