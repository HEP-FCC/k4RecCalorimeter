#ifndef RECCALORIMETER_NOISECALOCELLSFLATTOOL_H
#define RECCALORIMETER_NOISECALOCELLSFLATTOOL_H

// Gaudi
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"

// k4FWCore
#include "k4Interface/INoiseCaloCellsTool.h"

/** @class NoiseCaloCellsFlatTool
 *
 *  Very simple tool for calorimeter noise using a single noise value for all cells
 *  createRandomCellNoise: Create random CaloHits (gaussian distribution) for the vector of cells
 *  filterCellNoise: remove cells with energy below threshold*sigma from the vector of cells
 *
 *  @author Jana Faltova
 *  @date   2016-09
 *
 *  @author Giovanni Marchiori
 *  @date   2024-07
 */

class NoiseCaloCellsFlatTool : public extends<AlgTool,  INoiseCaloCellsTool> {
public:
  using base_class::base_class;
  virtual ~NoiseCaloCellsFlatTool() = default;
  virtual StatusCode initialize() override final;

  /** @brief Create random CaloHits (gaussian distribution) for the vector of cells (aCells).
   * Vector of cells must contain all cells in the calorimeter with their cellIDs.
   */
  virtual void addRandomCellNoise(std::unordered_map<uint64_t, double>& aCells) const override final;

  /** @brief Create random CaloHits (gaussian distribution) for the vector of cells (aCells).
   * Vector of cells must contain all cells in the calorimeter with their cellIDs.
   */
  virtual void addRandomCellNoise(std::vector<std::pair<uint64_t, double> >& aCells) const override final;

  /** @brief Remove cells with energy below threshold*sigma from the vector of cells
   */
  virtual void filterCellNoise(std::unordered_map<uint64_t, double>& aCells) const override final;

  /** @brief Remove cells with energy below threshold*sigma from the vector of cells
   */
  virtual void filterCellNoise(std::vector<std::pair<uint64_t, double> >& aCells) const override final;


private:
  template <typename C>
  void addRandomCellNoiseT (C& aCells) const;
  template <typename C>
  void filterCellNoiseT (C& aCells) const;

  /// RMS of noise -- uniform RMS per cell in GeV
  Gaudi::Property<double> m_cellNoiseRMS{this, "cellNoiseRMS", 0.003, "uniform noise RMS per cell in GeV"};
  /// Offset of noise -- uniform offset per cell in GeV
  Gaudi::Property<double> m_cellNoiseOffset{this, "cellNoiseOffset", 0.0, "uniform noise offset per cell in GeV"};
  /// Energy threshold (Ecell < m_cellNoiseOffset + filterThreshold*m_cellNoiseRMS removed)
  Gaudi::Property<double> m_filterThreshold{this, "filterNoiseThreshold", 3,
                                            "remove cells with energy below offset + threshold * noise RMS"};
  /// Random Number Service
  SmartIF<IRndmGenSvc> m_randSvc;
  /// Gaussian random number generator used for smearing with a constant resolution (m_sigma)
  Rndm::Numbers m_gauss;
};

#endif /* RECCALORIMETER_NOISECALOCELLSFLATTOOL_H */
