// This file's extension implies that it's C, but it's really -*- C++ -*-.
/**
 * @file k4RecCalorimeter/RecCalorimeter/src/components/CaloCellPositionsTool.h
 * @author scott snyder <snyder@bnl.gov>
 * @date Apr, 2025
 * @brief Generic tool to find positions of calorimeter cells.
 */


#ifndef RECCALORIMETER_CALOCELLPOSITIONSTOOL_H
#define RECCALORIMETER_CALOCELLPOSITIONSTOOL_H


#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ServiceHandle.h"
#include "RecCaloInterface/ICellPositionsTool.h"
#include "k4Interface/IGeoSvc.h"
#include "DD4hep/Segmentations.h"
#include "DD4hep/Volumes.h"


/** Generic tool to find positions of calorimeter cells.
 */
class CaloCellPositionsTool
  : public extends<AlgTool, ICellPositionsTool>
{
public:
  using base_class::base_class;

  /** Standard Gaudi initialize method.
   */
  virtual StatusCode initialize() override;

  /** Return the cartesian global coordinates of a cell center
   */
  virtual dd4hep::Position xyzPosition(const uint64_t& aCellId) const final override;

  /** Copy cells from aCells to outputColl filling in cell positions.
   */
  virtual void getPositions(const edm4hep::CalorimeterHitCollection& aCells,
                            edm4hep::CalorimeterHitCollection&       outputColl) const final override;

  /** Return the layer number of a cell.
   */
  virtual int              layerId(const uint64_t& aCellId)     const final override;


private:
  Gaudi::Property<std::string> m_readoutName
  { this, "readoutName", "", "Name of the readout for this detector" };
  Gaudi::Property<std::string> m_layerFieldName
  { this, "layerFieldName", "layer", "Name of the decoder field for layer" };

  ServiceHandle<IGeoSvc> m_geoSvc
  { this, "GeoSvc", "GeoSvc", "Geometry service" };

  // DD4hep volume manager.
  dd4hep::VolumeManager m_volman;

  // Segmentation for this geometry.
  dd4hep::Segmentation m_segmentation;

  // Decoder for this segmentation.
  const dd4hep::BitFieldCoder* m_decoder = nullptr;

  // Index of the layer field.
  size_t m_layerFieldIdx = 0;
};


#endif // not RECCALORIMETER_CALOCELLPOSITIONSTOOL_H
