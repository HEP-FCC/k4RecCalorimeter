#ifndef RECCALORIMETER_LAYERPHIETACALOTOOL_H
#define RECCALORIMETER_LAYERPHIETACALOTOOL_H

#include "RecCaloCommon/CalorimeterToolBase.h"

/** @class LayerPhiEtaCaloTool Reconstruction/RecCalorimeter/src/components/LayerPhiEtaCaloTool.h
 *LayerPhiEtaCaloTool.h
 *
 *  Tool for geometry-dependent settings of the digitisation.
 *  Needs the eta range per layer as input to calculate the total number of calo cells, layers and phi-eta segmentation
 *is taken from the geometry.
 *
 *  @author Anna Zaborowska, Coralie Neubueser
 */

class LayerPhiEtaCaloTool : public CalorimeterToolBase {
public:
  using CalorimeterToolBase::CalorimeterToolBase;
  virtual ~LayerPhiEtaCaloTool() = default;


protected:
  /** Fill vector with all existing cells for this geometry.
   */
  virtual StatusCode collectCells(std::vector<uint64_t>& cells) const final;


private:
  /// Name of active volumes
  Gaudi::Property<std::string> m_activeVolumeName{this, "activeVolumeName", "layerVolume"};
  /// Name of active layers for sampling calorimeter
  Gaudi::Property<std::string> m_activeFieldName{this, "activeFieldName", "layer"};
  /// Name of the fields describing the segmented volume
  Gaudi::Property<std::vector<std::string>> m_fieldNames{this, "fieldNames", {"system"}};
  /// Values of the fields describing the segmented volume
  Gaudi::Property<std::vector<int>> m_fieldValues{this, "fieldValues", 8};
  /// Temporary: for use with MergeLayer tool
  Gaudi::Property<unsigned int> m_activeVolumesNumber{this, "activeVolumesNumber", 10};
  /// Temporary: for use with Tile Calo
  Gaudi::Property<std::vector<double>> m_activeVolumesEta{
      this, "activeVolumesEta", {1.2524, 1.2234, 1.1956, 1.15609, 1.1189, 1.08397, 1.0509, 0.9999, 0.9534, 0.91072}};
};

#endif /* RECCALORIMETER_LAYERPHIETACALOTOOL_H */
