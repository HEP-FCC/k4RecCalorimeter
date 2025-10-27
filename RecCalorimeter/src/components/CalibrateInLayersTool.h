#ifndef RECCALORIMETER_CALIBRATEINLAYERSTOOL_H
#define RECCALORIMETER_CALIBRATEINLAYERSTOOL_H

// from Gaudi
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ServiceHandle.h"

#include "DDSegmentation/BitFieldCoder.h"

// k4FWCore
#include "k4Interface/ICalibrateCaloHitsTool.h"
class IGeoSvc;

#include <vector>

/** @class CalibrateInLayersTool Reconstruction/RecCalorimeter/src/components/CalibrateInLayersTool.h
 * CalibrateInLayersTool.h
 *
 *  Tool for energy calibration to the electromagnetic scale.
 *  Geant4 energy deposits calibration in sampling calorimeters.
 *
 *  Sampling fraction (energy deposited in the active material to the total deposited energy) is defined for layers of
 *  the calorimeter.
 *  It is used for calorimeters with varying sampling fraction in radial distance.
 *  Sampling fraction for each layer is used for the calibration of the deposits within this layer.
 *  Sampling fractions for layers may be obtained using SamplingFractionInLayers algorithm.
 *
 *  @author Anna Zaborowska
 */

class CalibrateInLayersTool : public extends<AlgTool, ICalibrateCaloHitsTool> {
public:
  using base_class::base_class;
  ~CalibrateInLayersTool() = default;
  /**  Initialize.
   *   @return status code
   */
  virtual StatusCode initialize() final;

  /** @brief  Calibrate Geant4 hit energy to EM scale
   */
  virtual void calibrate(std::unordered_map<uint64_t, double>& aHits) const final;
  virtual void calibrate(std::unordered_map<uint64_t, double>& aHits) final
  { const auto* cthis = this;  cthis->calibrate(aHits); }
  virtual void calibrate(std::vector<std::pair<uint64_t, double> >& aHits) const final;

private:
  /// Do calibration for a single cell.
  void calibrateCell (uint64_t cID, double& energy) const;

  /// Decoder associated with this readout.
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder = nullptr;

  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc"};
  /// Name of the detector readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "", "Name of the detector readout"};
  /// Name of the cells/layer field
  Gaudi::Property<std::string> m_layerFieldName{this, "layerFieldName", "", "Identifier of layers"};
  /// Id of the first layer (current design starts layer ids at 1)
  Gaudi::Property<uint> m_firstLayerId{this, "firstLayerId", 0, "ID of first layer"};
  /// Values of sampling fraction
  Gaudi::Property<std::vector<double>> m_samplingFraction{
      this, "samplingFraction", {}, "Values of sampling fraction per layer"};
};
#endif /* RECCALORIMETER_CALIBRATEINLAYERSTOOL_H */
