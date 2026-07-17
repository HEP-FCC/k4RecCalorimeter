#ifndef RECFCCEECALORIMETER_NOISECALOCELLSFROMFILEBASETOOL_H
#define RECFCCEECALORIMETER_NOISECALOCELLSFROMFILEBASETOOL_H

// from Gaudi
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"

// Interfaces
#include "RecCaloCommon/ICalorimeterTool.h"
#include "RecCaloCommon/INoiseCaloCellsTool.h"
#include "RecCaloCommon/INoiseConstTool.h"

#include "RecCaloCommon/ReadNoiseFromFileBaseTool.h"

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

class NoiseCaloCellsFromFileBaseTool : public extends<ReadNoiseFromFileBaseTool, k4::recCalo::INoiseCaloCellsTool> {
public:
  using CellID = k4::recCalo::INoiseCaloCellsTool::CellID;

  using base_class::base_class;
  virtual ~NoiseCaloCellsFromFileBaseTool() = default;
  virtual StatusCode initialize() override;

  /** @brief Create random CaloHits (gaussian distribution) for the vector of cells (aCells).
   * Vector of cells must contain all cells in the calorimeter with their cellIDs.
   */
  virtual void addRandomCellNoise(std::unordered_map<CellID, double>& aCells) const override final;

  /** @brief Create random CaloHits (gaussian distribution) for the vector of cells (aCells).
   * Vector of cells must contain all cells in the calorimeter with their cellIDs.
   */
  virtual void addRandomCellNoise(std::vector<std::pair<CellID, double>>& aCells) const override final;

  /** @brief Remove cells with energy below threshold*sigma from the vector of cells
   */
  virtual void filterCellNoise(std::unordered_map<CellID, double>& aCells) const override final;

  /** @brief Remove cells with energy below threshold*sigma from the vector of cells
   */
  virtual void filterCellNoise(std::vector<std::pair<CellID, double>>& aCells) const override final;

private:
  template <typename C>
  void addRandomCellNoiseT(C& aCells) const;
  template <typename C>
  void filterCellNoiseT(C& aCells) const;

  /// Name of electronics noise RMS histogram
  Gaudi::Property<std::string> m_elecNoiseRMSHistoName{this, "elecNoiseRMSHistoName", "h_elecNoise_layer",
                                                       "Name of electronics noise RMS histogram"};

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

  /// Random Number Service
  SmartIF<IRndmGenSvc> m_randSvc;
  /// Gaussian random number generator used for the generation of random noise hits
  Rndm::Numbers m_gauss;
};

#endif /* RECFCCEECALORIMETER_NOISECALOCELLSVSTHETAFROMFILEBASETOOL_H */
