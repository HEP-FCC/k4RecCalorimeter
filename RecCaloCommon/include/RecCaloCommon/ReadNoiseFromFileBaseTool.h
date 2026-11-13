// This file's extension implies that it's C, but it's really -*- C++ -*-.
/**
 * @file RecCaloCommon/ReadNoiseFromFileBaseTool.h
 * @author scott snyder <snyder@bnl.gov>, Jana Faltova, Coralie Neubueser, Giovanni Marchior
 * @date Jul, 2026
 * @brief Base class for reading noise data taken from ROOT histograms
 *        read from a file.
 */

#ifndef RECCALOCOMMON_READNOISEFROMFILEBASETOOL_H
#define RECCALOCOMMON_READNOISEFROMFILEBASETOOL_H

// from Gaudi
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolHandle.h"

// k4FWCore
#include "k4FWCore/DataHandle.h"

// Interfaces
#include "RecCaloCommon/ICaloCellConstantsSvc.h"
#include "RecCaloCommon/ICaloCellIndexerSvc.h"
#include "RecCaloCommon/INoiseConstTool.h"

#include <memory>

class IGeoSvc;

// Root
class TH1;

/** @class ReadNoiseFromFileBaseTool
 *
 *  Tool to read the stored noise constant per cell in the calorimeters
 *  Access noise constants from TH1F histogram (noise vs. |eta| or |theta|)
 *
 *  @author Jana Faltova, Coralie Neubueser
 *  @date   2018-01
 *
 * Returns noise RMS/offset per cell.
 *
 * Noise information is taken from histograms read from the ROOT file
 * given by noiseFileName.  There are up to to four sets of histograms,
 * for electronic noise RMS, pileup noise RMS, electronic noise offset,
 * and pileup noise offset.  Offset histograms are present if setNoiseOffset
 * is true, and pileup histograms are present if addPileup is true.
 * Total rms is calculated as the quadrature sum of electronic
 * and pileup contributions, multiplied by m_scaleFactor; total offset
 * is the linear sum of the same.
 *
 * Each set of histograms is provided as an indexed list.  The names are
 * of the form HNAME1, HNAME2, etc., where the 1 in the name corresponds
 * to index 0.  The index used is taken from the CellID field
 * activeFieldName.  This would usually be a layer index, but not
 * necessarily.  The number of histograms within each set is given
 * by numHistograms.
 *
 * All the histograms are held a NoiseData structure, registered with
 * CaloCellConstantsSvc.
 *
 * The NoiseData structure also contains a mapping from CellIDs to histogram
 * bin indices (separately for noise and offset histograms).
 * This is filled in by the initBinning virtual method.
 */
class ReadNoiseFromFileBaseTool : public extends<AlgTool, k4::recCalo::INoiseConstTool> {
public:
  using base_class::base_class;
  virtual ~ReadNoiseFromFileBaseTool() = default;

  virtual StatusCode initialize() override;

  /// Return noise RMS for a given CellID.
  virtual double getNoiseRMSPerCell(CellID aCellID) const override final;

  /// Return noise offset for a given CellID.
  virtual double getNoiseOffsetPerCell(CellID aCellID) const override final;

  /// Return both noise RMS (first) and offset (second)  for a given CellID.
  virtual std::pair<double, double> getNoisePerCell(CellID aCellID) const override final;

protected:
  // clang-format off
  /// Name of the detector readout
  Gaudi::Property<std::string> m_readoutName
    {this, "readoutName", "", "Name of the detector readout"};
  /// Noise offset, if false, mean is set to 0
  Gaudi::Property<bool> m_setNoiseOffset
    {this, "setNoiseOffset", true, "Set a noise offset per cell"};
  /// Name of electronics noise histogram
  Gaudi::Property<std::string> m_elecNoiseHistoName
    {this, "elecNoiseHistoName", "h_elecNoise_layer", "Name of electronics noise histogram"};
  /// Handle to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc
    {this, "GeoSvc", "GeoSvc"};
  // clang-format on

  /// Noise data read from the file.
  struct NoiseData {
    NoiseData() = default;
    NoiseData(const NoiseData&);
    NoiseData(NoiseData&&) = default;
    NoiseData& operator=(const NoiseData&) = delete;
    NoiseData& operator=(NoiseData&&) = default;
    /// Histograms with pileup RMS (index in array - radial layer)
    std::vector<std::unique_ptr<TH1>> m_histoPileupNoiseRMS;
    /// Histograms with electronics noise RMS (index in array - radial layer)
    std::vector<std::unique_ptr<TH1>> m_histoElecNoiseRMS;

    /// Histograms with pileup offset (index in array - radial layer)
    std::vector<std::unique_ptr<TH1>> m_histoPileupNoiseOffset;
    /// Histograms with electronics noise offset (index in array - radial layer)
    std::vector<std::unique_ptr<TH1>> m_histoElecNoiseOffset;

    // Histogram bin indices (for RMS and offset) for each CellID.
    std::vector<std::pair<unsigned, unsigned>> m_bins;
  };

  /// Initialize the m_bins member of NoiseData.
  virtual StatusCode initBinning(NoiseData& data, const k4::recCalo::ICaloIndexer& indexer) const = 0;

  // Decoder
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  int m_index_activeField = -1;

private:
  /// Return noise RMS for a given bin/layer.
  double getNoiseRMSPerCell(int ibin, unsigned cellLayer) const;

  /// Return noise offset for a given bin/layer.
  double getNoiseOffsetPerCell(int ibin, unsigned cellLayer) const;

  /// Open file and read noise histograms into a NoiseData structure.
  StatusCode initNoiseFromFile(NoiseData& data) const;

  // clang-format off
  /// Add pileup contribution to the electronics noise? (only if read from file)
  Gaudi::Property<bool> m_addPileup
    {this, "addPileup", true, "Add pileup contribution to the electronics noise? (only if read from file)"};

  /// Name of the file with noise constants
  Gaudi::Property<std::string> m_noiseFileName
    {this, "noiseFileName", "", "Name of the file with noise constants"};
  /// Name of active layers for sampling calorimeter
  Gaudi::Property<std::string> m_activeFieldName
    {this, "activeFieldName", "active_layer", "Name of active layers for sampling calorimeter"};
  /// Name of pileup histogram
  Gaudi::Property<std::string> m_pileupHistoName
    {this, "pileupHistoName", "h_pileup_layer", "Name of pileup histogram"};
  /// Name of electronics noise offset histogram
  Gaudi::Property<std::string> m_elecNoiseOffsetHistoName
    {this, "elecNoiseOffsetHistoName", "h_mean_pileup_layer", "Name of electronics noise offset histogram"};
  /// Name of pileup offset histogram
  Gaudi::Property<std::string> m_pileupOffsetHistoName
    {this, "pileupOffsetHistoName", "h_pileup_layer", "Name of pileup offset histogram"};

  /// For config compatibility.
  Gaudi::Property<uint> m_numRadialLayers
    {this, "numRadialLayers", 0, "Number of radial layers.  Deprecated: use numHistograms instead."};
  /// Number of radial layers/wheels
  Gaudi::Property<uint> m_numHistograms
    {this, "numHistograms", 3, "Number of histograms"};

  /// Factor to apply to the noise values to get them in GeV if e.g. they were produced in MeV
  Gaudi::Property<float> m_scaleFactor
    {this, "scaleFactor", 1, "Factor to apply to the noise values"};

  const NoiseData* m_data = nullptr;

  ServiceHandle<k4::recCalo::ICaloCellIndexerSvc> m_indexerSvc
    { this, "CaloCellIndexerSvc", "k4::recCalo::CaloCellIndexerSvc", "" };
  ServiceHandle<k4::recCalo::ICaloCellConstantsSvc> m_constantsSvc
    { this, "CaloCellConstantsSvc", "k4::recCalo::CaloCellConstantsSvc", "" };
  // clang-format on

  const k4::recCalo::ICaloIndexer* m_indexer = nullptr;
};

#endif // not RECCALOCOMMON_READNOISEFROMFILEBASETOOL_H
