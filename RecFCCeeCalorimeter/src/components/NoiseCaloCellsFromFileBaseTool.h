#ifndef RECFCCEECALORIMETER_NOISECALOCELLSFROMFILEBASETOOL_H
#define RECFCCEECALORIMETER_NOISECALOCELLSFROMFILEBASETOOL_H

// from Gaudi
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"

// Interfaces
#include "RecCaloCommon/ICellPositionsTool.h"
#include "RecCaloCommon/INoiseCaloCellsTool.h"
#include "RecCaloCommon/INoiseConstTool.h"

// k4FWCore
#include "k4Interface/IGeoSvc.h"
// class IGeoSvc;

// Root
#include "TH1F.h"
class TFile;


/** @class NoiseCaloCellsFromFileBaseTool
 *
 *  Common base class for calorimeter noise-from-file tools
 *  Access noise constants from histograms saved in ROOT file
 *  createRandomCellNoise: Create random CaloHits (gaussian distribution) for the vector of cells
 *  filterCellNoise: remove cells with energy below threshold*sigma from the vector of cells
 *  The tool needs a cell positioning tool to translate cellID to position
 *  Geometry-specific maaping is delegated to derived classes
 *
 *  @author Giovanni Marchiori
 *  @date   2026-06
 *
 */

class NoiseCaloCellsFromFileBaseTool : public extends<AlgTool, k4::recCalo::INoiseCaloCellsTool,  k4::recCalo::INoiseConstTool> {
public:
  using CellID = k4::recCalo::INoiseCaloCellsTool::CellID;

  using base_class::base_class;
  virtual ~NoiseCaloCellsFromFileBaseTool() = default;
  StatusCode initialize() override;

  /** @brief Create random CaloHits (gaussian distribution) for the vector of cells (aCells).
   * Vector of cells must contain all cells in the calorimeter with their cellIDs.
   */
  virtual void addRandomCellNoise(std::unordered_map<CellID, double>& aCells) const override final;

  /** @brief Create random CaloHits (gaussian distribution) for the vector of cells (aCells).
   * Vector of cells must contain all cells in the calorimeter with their cellIDs.
   */
  virtual void addRandomCellNoise(std::vector<std::pair<CellID, double> >& aCells) const override final;

  /** @brief Remove cells with energy below threshold*sigma from the vector of cells
   */
  virtual void filterCellNoise(std::unordered_map<CellID, double>& aCells) const override final;

  /** @brief Remove cells with energy below threshold*sigma from the vector of cells
   */
  virtual void filterCellNoise(std::vector<std::pair<CellID, double> >& aCells) const override final;

  /// Open file and read noise histograms in the memory
  StatusCode initNoiseFromFile();
  StatusCode readHistograms(TFile* noiseFile);
  /// Find the appropriate noise RMS from the histogram
  virtual double getNoiseRMSPerCell(CellID aCellID) const override final;
  virtual double getNoiseOffsetPerCell(CellID aCellID) const override final;
  virtual std::pair<double, double> getNoisePerCell(CellID aCellID) const override final;

protected:
  /// Handle for tool to get cell positions - available also to derived classes
  ToolHandle<k4::recCalo::ICellPositionsTool> m_cellPositionsTool{this, "cellPositionsTool", "CellPositionsDummyTool",
      "Handle for tool to retrieve cell positions"};
  // Decoder for readout
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  /// index of layer or wheel
  int m_index_activeField;
  /// Histograms with pileup noise RMS (index in array - radial layer)
  std::vector<TH1*> m_histoPileupNoiseRMS;
  /// Histograms with electronics noise RMS (index in array - radial layer)
  std::vector<TH1*> m_histoElecNoiseRMS;
  /// Histograms with pileup offset (index in array - radial layer)
  std::vector<TH1*> m_histoPileupNoiseOffset;
  /// Histograms with electronics noise offset (index in array - radial layer)
  std::vector<TH1*> m_histoElecNoiseOffset;

private:
  template <typename C>
  void addRandomCellNoiseT (C& aCells) const;
  template <typename C>
  void filterCellNoiseT (C& aCells) const;

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
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "ECalBarrelModuleThetaMerged",
                                             "Name of the detector readout"};
  /// Name of active layers/wheels for barrel/endcap sampling calorimeter
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
  Gaudi::Property<double> m_filterThreshold{
      this, "filterNoiseThreshold", 3,
      "Energy threshold (cells with Ecell < offset + filterThreshold*m_cellNoise removed)"};
  /// Change the cell filter condition to remove only cells with abs(Ecell) < offset + filterThreshold*m_cellNoise
  /// removed) This avoids to keep only 'one side'  of the noise fluctuations and prevents biasing cluster energy
  /// towards higher energies
  Gaudi::Property<bool> m_useAbsInFilter{
      this, "useAbsInFilter", false,
      "If set, cell filtering condition becomes: drop cell if abs(Ecell-offset) < filterThreshold*m_cellNoise"};
  /// Number of radial layers/wheels
  Gaudi::Property<uint> m_numHistograms{this, "numHistograms", 0, "Number of histograms"};

  /// Random Number Service
  SmartIF<IRndmGenSvc> m_randSvc;
  /// Gaussian random number generator used for the generation of random noise hits
  Rndm::Numbers m_gauss;

  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc { this, "GeoSvc", "GeoSvc" };

  /// get index of histogram for given cellID
  unsigned getIndexHistogram(CellID aCellId) const;

  /// get bin in histogram for given cellID
  virtual unsigned getBin(CellID aCellId) const = 0;
};

#endif /* RECFCCEECALORIMETER_NOISECALOCELLSVSTHETAFROMFILEBASETOOL_H */
