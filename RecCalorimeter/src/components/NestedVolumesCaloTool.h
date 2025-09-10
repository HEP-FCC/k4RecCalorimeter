#ifndef RECCALORIMETER_NESTEDVOLUMESCALOTOOL_H
#define RECCALORIMETER_NESTEDVOLUMESCALOTOOL_H

#include "RecCaloCommon/CalorimeterToolBase.h"

/** @class NestedVolumesCaloTool Reconstruction/RecCalorimeter/src/components/NestedVolumesCaloTool.h
 *NestedVolumesCaloTool.h
 *
 *  Tool for geometry-dependent settings of the digitisation.
 *  It assumes no segmentation is used. It may be used for nested volumes.
 *   For more explanation please [see reconstruction documentation](@ref md_reconstruction_doc_reccalorimeter).
 *
 *  @author Anna Zaborowska
 */

class NestedVolumesCaloTool : public CalorimeterToolBase {
public:
  using CalorimeterToolBase::CalorimeterToolBase;
  virtual ~NestedVolumesCaloTool() = default;


protected:
  /** Fill vector with all existing cells for this geometry.
   */
  virtual StatusCode collectCells(std::vector<uint64_t>& cells) const final;


private:
  /// Name of active volumes (if different than all)
  Gaudi::Property<std::vector<std::string>> m_activeVolumeName{
      this, "activeVolumeName", {"LAr_sensitive"}, "Name of active volumes (if different than all)"};
  /// Name of active layers for sampling calorimeter
  Gaudi::Property<std::vector<std::string>> m_activeFieldName{
      this, "activeFieldName", {"active_layer"}, "Name of active layers for sampling calorimeter"};
  /// Name of the fields describing the segmented volume
  Gaudi::Property<std::vector<std::string>> m_fieldNames{
      this, "fieldNames", {}, "Name of the fields describing the segmented volume"};
  /// Values of the fields describing the segmented volume
  Gaudi::Property<std::vector<int>> m_fieldValues{
      this, "fieldValues", {}, "Values of the fields describing the segmented volume"};
};

#endif /* RECCALORIMETER_NESTEDVOLUMESCALOTOOL_H */
