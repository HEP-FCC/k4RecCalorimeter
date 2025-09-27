#ifndef RECCALORIMETER_CONSTNOISETOOL_H
#define RECCALORIMETER_CONSTNOISETOOL_H

// from Gaudi
#include "GaudiKernel/AlgTool.h"
class IRndmGenSvc;

// DD4HEP
#include "DDSegmentation/BitFieldCoder.h"

// k4FWCore
#include "k4Interface/INoiseConstTool.h"
class IGeoSvc;

#include <map>

/** @class ConstNoiseTool
 *
 *  Tool to set noise offset and RMS for all cells in Calorimeters.
 *  By default, set to rough estimates from ATLAS, can be changed in arguments for each Calo sub-system.
 *
 *  @author Coralie Neubueser
 *  @date   2018-01
 *
 *  @author Giovanni Marchiori
 *  @date   2024-07
 *
 */

class ConstNoiseTool : public AlgTool, virtual public INoiseConstTool {
public:
  ConstNoiseTool(const std::string& type, const std::string& name, const IInterface* parent);
  virtual ~ConstNoiseTool() = default;

  virtual StatusCode initialize() final;

  virtual StatusCode finalize() final;

  /// Find the appropriate noise constant from the histogram
  virtual double getNoiseRMSPerCell(uint64_t aCellID) const final;
  virtual double getNoiseRMSPerCell(uint64_t aCellID) final
  { const auto* cthis = this;  return cthis->getNoiseRMSPerCell(aCellID); }
  virtual double getNoiseOffsetPerCell(uint64_t aCellID) const final;
  virtual double getNoiseOffsetPerCell(uint64_t aCellID) final
  { const auto* cthis = this;  return cthis->getNoiseOffsetPerCell(aCellID); }

private:
  std::map<uint, double> m_systemNoiseRMSMap;
  std::map<uint, double> m_systemNoiseOffsetMap;

  /// List of subdetector names (they must match what is defined in DectDimension)
  Gaudi::Property<std::vector<std::string>> m_detectors{
      this,
      "detectors",
      {"ECAL_Barrel", "ECAL_Endcap", "HCal_Endcap", "HCalBarrel", "HCalExtBarrel", "ECalFwd", "HCalFwd"},
      "names of the calorimeters"};

  /// Cell noise RMS (=1Sigma threshold) in GeV
  /// default values estimates for the LAr and TileCal from ATLAS Barrel: ATL-LARG-PUB-2008_002
  /// effective seed thresholds in topo-clustering of 7.5 and 11.5MeV in LAr and TileCal
  Gaudi::Property<std::vector<double>> m_detectorsNoiseRMS{
      this,
      "detectorsNoiseRMS",
      {0.0075 / 4, 0.0075 / 4, 0.0115 / 4, 0.0115 / 4, 0.0115 / 4, 0.0075 / 4, 0.0075 / 4},
      "Cell noise RMS in GeV"};

  /// Set noise offset or not (false by default)
  Gaudi::Property<bool> m_setNoiseOffset{this, "setNoiseOffset", false, "Set a noise offset per cell"};

  /// Cell noise offset in GeV. 0 by default
  Gaudi::Property<std::vector<double>> m_detectorsNoiseOffset{
      this, "detectorsNoiseOffset", {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, "Cell noise offset in GeV"};

  /// System encoding string
  Gaudi::Property<std::string> m_systemEncoding{this, "systemEncoding", "system:4", "System encoding string"};

  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;

  /// Decoder
  std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_decoder;
};

#endif /* RECCALORIMETER_CONSTNOISETOOL_H */
