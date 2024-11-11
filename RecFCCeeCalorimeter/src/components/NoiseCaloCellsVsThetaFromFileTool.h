#ifndef RECFCCEECALORIMETER_NOISECALOCELLSVSTHETAFROMFILETOOL_H
#define RECFCCEECALORIMETER_NOISECALOCELLSVSTHETAFROMFILETOOL_H

// from Gaudi
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/IRndmGenSvc.h"

// k4geo
// #include "detectorSegmentations/FCCSWGridPhiEta_k4geo.h"

// k4FWCore
#include "k4Interface/INoiseConstTool.h"
#include "k4Interface/INoiseCaloCellsTool.h"
#include "k4Interface/ICellPositionsTool.h"

class IGeoSvc;

// Root
class TH1F;

/** @class NoiseCaloCellsVsThetaFromFileTool
 *
 *  Tool for calorimeter noise
 *  Access noise constants from TH1F histogram (noise vs. theta)
 *  createRandomCellNoise: Create random CaloHits (gaussian distribution) for the vector of cells
 *  filterCellNoise: remove cells with energy below threshold*sigma from the vector of cells
 * The tool needs a cell positioning tool to translate cellID to cell theta.
 * In alternative, the tool could be rewritten to use a specific segmentation class for the cellID->theta
 * translation, but it would be coupled to a specific readout.
 * To avoid all these calls to the positioning tool, one could either
 * - save directly the noise histograms as histos of noise vs thetaID
 * - or, keep histos of noise vs theta, but change the interfaces and the tool to accept
 *   cells rather than cellIDs as input. One would then get theta from the cells.
 *
 *  @author Giovanni Marchiori
 *  @date   2024-07
 *
 */

class NoiseCaloCellsVsThetaFromFileTool : public AlgTool, virtual public INoiseCaloCellsTool, virtual public INoiseConstTool {
public:
  NoiseCaloCellsVsThetaFromFileTool(const std::string& type, const std::string& name, const IInterface* parent);
  virtual ~NoiseCaloCellsVsThetaFromFileTool() = default;
  virtual StatusCode initialize() final;
  virtual StatusCode finalize() final;

  /** @brief Create random CaloHits (gaussian distribution) for the vector of cells (aCells).
   * Vector of cells must contain all cells in the calorimeter with their cellIDs.
   */
  virtual void addRandomCellNoise(std::unordered_map<uint64_t, double>& aCells) final;
  /** @brief Remove cells with energy below threshold*sigma from the vector of cells
   */
  virtual void filterCellNoise(std::unordered_map<uint64_t, double>& aCells) final;

  /// Open file and read noise histograms in the memory
  StatusCode initNoiseFromFile();
  /// Find the appropriate noise RMS from the histogram
  double getNoiseRMSPerCell(uint64_t aCellID);
  double getNoiseOffsetPerCell(uint64_t aCellID);

private:
  /// Handle for tool to get cell positions
  ToolHandle<ICellPositionsTool> m_cellPositionsTool{"CellPositionsDummyTool", this};

  /// Add pileup contribution to the electronics noise? (only if read from file)
  Gaudi::Property<bool> m_addPileup{this, "addPileup", true,
                                    "Add pileup contribution to the electronics noise? (only if read from file)"};
  /// Read noise offset from histograms in files or not (if false, offset is set to 0)
  Gaudi::Property<bool> m_setNoiseOffset{this, "setNoiseOffset", false, "Set a noise offset per cell"};
  /// Name of the file with noise constants
  Gaudi::Property<std::string> m_noiseFileName{this, "noiseFileName", "", "Name of the file with noise constants"};
  /// Factor to apply to the noise values to get them in GeV if e.g. they were produced in MeV
  Gaudi::Property<float> m_scaleFactor{this, "scaleFactor", 1, "Factor to apply to the noise values"};
  /// Name of the detector readout (only needed to retrieve the decoder)
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "ECalBarrelThetaModuleMerged",
                                             "Name of the detector readout"};
  /// Name of active layers for sampling calorimeter
  Gaudi::Property<std::string> m_activeFieldName{this, "activeFieldName", "active_layer",
                                                 "Name of active layers for sampling calorimeter"};
  /// Name of pileup noise RMS histogram
  Gaudi::Property<std::string> m_pileupNoiseRMSHistoName{this, "pileupNoiseRMSHistoName", "h_pileup_layer",
                                                         "Name of pileup noise RMS histogram"};
  /// Name of electronics noise RMS histogram
  Gaudi::Property<std::string> m_elecNoiseRMSHistoName{this, "elecNoiseRMSHistoName", "h_elecNoise_layer",
                                                       "Name of electronics noise RMS histogram"};
  /// Name of electronics noise offset histogram
  Gaudi::Property<std::string> m_elecNoiseOffsetHistoName{this, "elecNoiseOffsetHistoName", "h_mean_pileup_layer",
                                                          "Name of electronics noise offset histogram"};
  /// Name of pileup noise offset histogram
  Gaudi::Property<std::string> m_pileupNoiseOffsetHistoName{this, "pileupNoiseOffsetHistoName", "h_pileup_layer",
                                                            "Name of pileup noise offset histogram"};

  /// Energy threshold (cells with Ecell < filterThreshold*m_cellNoise removed)
  Gaudi::Property<double> m_filterThreshold{this, "filterNoiseThreshold", 3,
                                            "Energy threshold (cells with Ecell < offset + filterThreshold*m_cellNoise removed)"};
  /// Change the cell filter condition to remove only cells with abs(Ecell) < offset + filterThreshold*m_cellNoise removed)
  /// This avoids to keep only 'one side'  of the noise fluctuations and prevents biasing cluster energy towards higher energies
  Gaudi::Property<bool> m_useAbsInFilter{this, "useAbsInFilter", false,
                                         "If set, cell filtering condition becomes: drop cell if abs(Ecell-offset) < filterThreshold*m_cellNoise"};
  /// Number of radial layers
  Gaudi::Property<uint> m_numRadialLayers{this, "numRadialLayers", 3,
                                          "Number of radial layers"};

  /// Histograms with pileup noise RMS (index in array - radial layer)
  std::vector<TH1F> m_histoPileupNoiseRMS;
  /// Histograms with electronics noise RMS (index in array - radial layer)
  std::vector<TH1F> m_histoElecNoiseRMS;
  /// Histograms with pileup offset (index in array - radial layer)
  std::vector<TH1F> m_histoPileupNoiseOffset;
  /// Histograms with electronics noise offset (index in array - radial layer)
  std::vector<TH1F> m_histoElecNoiseOffset;

  /// Random Number Service
  SmartIF<IRndmGenSvc> m_randSvc;
  /// Gaussian random number generator used for the generation of random noise hits
  Rndm::Numbers m_gauss;

  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;

  /// Decoder for ECal layers
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  int m_index_activeField;
};

#endif /* RECFCCEECALORIMETER_NOISECALOCELLSVSTHETAFROMFILETOOL_H */
