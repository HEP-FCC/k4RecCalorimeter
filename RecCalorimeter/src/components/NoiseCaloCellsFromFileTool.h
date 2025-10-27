#ifndef RECCALORIMETER_NOISECALOCELLSFROMFILETOOL_H
#define RECCALORIMETER_NOISECALOCELLSFROMFILETOOL_H

// from Gaudi
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"

// k4geo
#include "detectorSegmentations/FCCSWGridPhiEta_k4geo.h"

// k4FWCore
#include "k4Interface/ICellPositionsTool.h"
#include "k4Interface/INoiseCaloCellsTool.h"
class IGeoSvc;

// DD4hep
#include "DDSegmentation/MultiSegmentation.h"

// Root
class TH1F;

/** @class NoiseCaloCellsFromFileTool
 *
 *  Tool for calorimeter noise
 *  Access noise constants from TH1F histogram (noise vs. |eta|)
 *  createRandomCellNoise: Create random CaloHits (gaussian distribution) for the vector of cells
 *  filterCellNoise: remove cells with energy bellow threshold*sigma from the vector of cells
 *
 *  @author Jana Faltova
 *  @date   2016-09
 *
 */

class NoiseCaloCellsFromFileTool : public extends<AlgTool, INoiseCaloCellsTool> {
public:
  using base_class::base_class;
  virtual ~NoiseCaloCellsFromFileTool() = default;
  virtual StatusCode initialize() final;

  /** @brief Create random CaloHits (gaussian distribution) for the vector of cells (aCells).
   * Vector of cells must contain all cells in the calorimeter with their cellIDs.
   */
  virtual void addRandomCellNoise(std::unordered_map<uint64_t, double>& aCells) const final;
  virtual void addRandomCellNoise(std::unordered_map<uint64_t, double>& aCells) final
  { const auto* cthis = this;  cthis->addRandomCellNoise(aCells); }

  /** @brief Create random CaloHits (gaussian distribution) for the vector of cells (aCells).
   * Vector of cells must contain all cells in the calorimeter with their cellIDs.
   */
  virtual void addRandomCellNoise(std::vector<std::pair<uint64_t, double> >& aCells) const final;

  /** @brief Remove cells with energy bellow threshold*sigma from the vector of cells
   */
  virtual void filterCellNoise(std::unordered_map<uint64_t, double>& aCells) const final;
  virtual void filterCellNoise(std::unordered_map<uint64_t, double>& aCells) final
  { const auto* cthis = this;  cthis->filterCellNoise(aCells); }

  /** @brief Remove cells with energy below threshold*sigma from the vector of cells
   */
  virtual void filterCellNoise(std::vector<std::pair<uint64_t, double> >& aCells)    const final;

  /// Open file and read noise histograms in the memory
  StatusCode initNoiseFromFile();
  /// Find the appropriate noise RMS from the histogram
  double getNoiseRMSPerCell(uint64_t aCellID) const;

private:
  template <typename C>
  void addRandomCellNoiseT (C& aCells) const;
  template <typename C>
  void filterCellNoiseT (C& aCells) const;

  /// Handle for tool to get cell positions
  ToolHandle<ICellPositionsTool> m_cellPositionsTool
  { this, "cellPositionsTool", "CellPositionsDummyTool", "Handle for tool to retrieve cell positions" };

  /// Add pileup contribution to the electronics noise? (only if read from file)
  Gaudi::Property<bool> m_addPileup{this, "addPileup", true,
                                    "Add pileup contribution to the electronics noise? (only if read from file)"};
  /// use segmentation in case no cell psotion tool defined.
  Gaudi::Property<bool> m_useSeg{this, "useSegmentation", true,
                                 "Specify if segmentation is used to determine cell position."};
  /// Name of the file with noise constants
  Gaudi::Property<std::string> m_noiseFileName{this, "noiseFileName", "", "Name of the file with noise constants"};
  /// Factor to apply to the noise values to get them in GeV if e.g. they were produced in MeV
  Gaudi::Property<float> m_scaleFactor{this, "scaleFactor", 1, "Factor to apply to the noise values"};
  /// Name of the detector readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "ECalHitsPhiEta", "Name of the detector readout"};
  /// Name of active layers for sampling calorimeter
  Gaudi::Property<std::string> m_activeFieldName{this, "activeFieldName", "active_layer",
                                                 "Name of active layers for sampling calorimeter"};
  /// Name of pileup histogram
  Gaudi::Property<std::string> m_pileupHistoName{this, "pileupHistoName", "h_pileup_layer", "Name of pileup histogram"};
  /// Name of electronics noise histogram
  Gaudi::Property<std::string> m_elecNoiseHistoName{this, "elecNoiseHistoName", "h_elecNoise_layer",
                                                    "Name of electronics noise histogram"};
  /// Energy threshold (cells with Ecell < filterThreshold*m_cellNoise removed)
  Gaudi::Property<double> m_filterThreshold{
      this, "filterNoiseThreshold", 3, "Energy threshold (cells with Ecell < filterThreshold*m_cellNoise removed)"};
  /// Change the cell filter condition to remove only cells with abs(Ecell) < filterThreshold*m_cellNoise removed)
  /// This avoids to keep only 'one side'  of the noise fluctuations and prevents biasing cluster energy towards higher
  /// energies
  Gaudi::Property<bool> m_useAbsInFilter{
      this, "useAbsInFilter", false,
      "Cell filtering condition becomes: drop cell if abs(Ecell) < filterThreshold*m_cellNoise"};
  /// Number of radial layers
  Gaudi::Property<uint> m_numRadialLayers{this, "numRadialLayers", 3, "Number of radial layers"};
  /// Histograms with pileup RMS (index in array - radial layer)
  std::vector<TH1F> m_histoPileupNoiseRMS;
  /// Histograms with electronics noise RMS (index in array - radial layer)
  std::vector<TH1F> m_histoElecNoiseRMS;

  /// Random Number Service
  SmartIF<IRndmGenSvc> m_randSvc;
  /// Gaussian random number generator used for the generation of random noise hits
  Rndm::Numbers m_gauss;

  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc { this, "GeoSvc", "GeoSvc" };
  /// PhiEta segmentation
  dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo* m_segmentationPhiEta;
  /// Multi segmentation
  dd4hep::DDSegmentation::MultiSegmentation* m_segmentationMulti;

  /// Decoder for ECal layers
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  int m_index_activeField;
};

#endif /* RECCALORIMETER_NOISECALOCELLSFROMFILETOOL_H */
