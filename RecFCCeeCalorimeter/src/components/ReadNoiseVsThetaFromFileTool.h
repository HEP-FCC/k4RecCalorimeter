#ifndef RECFCCEECALORIMETER_READNOISEVSTHETAFROMFILETOOL_H
#define RECFCCEECALORIMETER_READNOISEVSTHETAFROMFILETOOL_H

// Gaudi
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/ToolHandle.h"

// k4FWCore
#include "k4Interface/INoiseConstTool.h"
#include "k4Interface/ICellPositionsTool.h"

// k4geo
#include "detectorSegmentations/FCCSWGridModuleThetaMerged_k4geo.h"

class IGeoSvc;

// Root
class TH1F;

/** @class ReadNoiseVsThetaFromFileTool
 *
 *  Tool to read the stored noise constants per cell in the calorimeters
 *  Access noise constants from TH1F histogram (noise vs. theta)
 *  Cell positioning can be done via segmentation or via a cell 
 *  positioning tool
 *
 *  @author Tong Li, Giovanni Marchiori
 *  @date   2023-09
 *
 */

class ReadNoiseVsThetaFromFileTool : public GaudiTool, virtual public INoiseConstTool {
public:
  ReadNoiseVsThetaFromFileTool(const std::string& type, const std::string& name, const IInterface* parent);
  virtual ~ReadNoiseVsThetaFromFileTool() = default;

  virtual StatusCode initialize() final;

  virtual StatusCode finalize() final;

  /// Open file and read noise histograms in the memory
  StatusCode initNoiseFromFile();
  /// Find the appropriate noise constant from the histogram
  double getNoiseRMSPerCell(uint64_t aCellID);
  double getNoiseOffsetPerCell(uint64_t aCellID);

private:
  /// Handle for tool to get cell positions (optional - not needed if segmentation class is used)
  /// Provided as a way to use the tool for other theta-based segmentations (using the positioning
  /// tool instead of the segmentation to compute the theta of a cell)
  /// (maybe it suffices to use a GridTheta segmentation base class?)
  ToolHandle<ICellPositionsTool> m_cellPositionsTool;
  /// use segmentation (theta-module only so far!) in case no cell position tool is defined
  Gaudi::Property<bool> m_useSeg{this, "useSegmentation", true, "Specify if segmentation is used to determine cell position"};
  /// Add pileup contribution to the electronics noise? (only if read from file)
  Gaudi::Property<bool> m_addPileup{this, "addPileup", true,
                                    "Add pileup contribution to the electronics noise? (only if read from file)"};
  /// Noise offset, if false, mean is set to 0
  Gaudi::Property<bool> m_setNoiseOffset{this, "setNoiseOffset", true, "Set a noise offset per cell"};
  /// Name of the file with noise constants
  Gaudi::Property<std::string> m_noiseFileName{this, "noiseFileName", "", "Name of the file with noise constants"};
  /// Name of the detector readout (needed to get the decoder to extract e.g. layer information)
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "ECalBarrelThetaModuleMerged", "Name of the detector readout"};
  /// Name of active layers for sampling calorimeter
  Gaudi::Property<std::string> m_activeFieldName{this, "activeFieldName", "active_layer",
                                                 "Name of active layers for sampling calorimeter"};
  /// Name of pileup histogram
  Gaudi::Property<std::string> m_pileupNoiseHistoName{this, "pileupHistoName", "h_pileup_layer", "Name of pileup histogram"};
  /// Name of electronics noise histogram
  Gaudi::Property<std::string> m_elecNoiseHistoName{this, "elecNoiseHistoName", "h_elecNoise_layer",
                                                    "Name of electronics noise histogram"};
  /// Name of electronics noise offset histogram
  Gaudi::Property<std::string> m_elecNoiseOffsetHistoName{this, "elecNoiseOffsetHistoName", "h_mean_pileup_layer",
                                                    "Name of electronics noise offset histogram"};
  /// Name of pileup offset histogram
  Gaudi::Property<std::string> m_pileupNoiseOffsetHistoName{this, "pileupOffsetHistoName", "h_pileup_layer", "Name of pileup offset histogram"};
  /// Number of radial layers
  Gaudi::Property<uint> m_numRadialLayers{this, "numRadialLayers", 3, "Number of radial layers"};
  /// Factor to apply to the noise values to get them in GeV if e.g. they were produced in MeV
  Gaudi::Property<float> m_scaleFactor{this, "scaleFactor", 1, "Factor to apply to the noise values"};
  /// Histograms with pileup constants (index in array - radial layer)
  std::vector<TH1F> m_histoPileupNoiseRMS;
  /// Histograms with electronics noise constants (index in array - radial layer)
  std::vector<TH1F> m_histoElecNoiseRMS;
  /// Histograms with pileup offset (index in array - radial layer)
  std::vector<TH1F> m_histoPileupNoiseOffset;
  /// Histograms with electronics noise offset (index in array - radial layer)
  std::vector<TH1F> m_histoElecNoiseOffset;
  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;
  /// Theta-Module segmentation
  dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo* m_segmentation;
  // Decoder for ECal layers
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

};

#endif /* RECFCCEECALORIMETER_READNOISEVSTHETAFROMFILETOOL_H */
