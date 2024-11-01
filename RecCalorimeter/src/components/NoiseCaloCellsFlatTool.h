#ifndef RECCALORIMETER_NOISECALOCELLSFLATTOOL_H
#define RECCALORIMETER_NOISECALOCELLSFLATTOOL_H

// Gaudi
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/IRndmGenSvc.h"

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

class NoiseCaloCellsFlatTool : public AlgTool, virtual public INoiseCaloCellsTool {
public:
  NoiseCaloCellsFlatTool(const std::string& type, const std::string& name, const IInterface* parent);
  virtual ~NoiseCaloCellsFlatTool() = default;
  virtual StatusCode initialize() final;
  virtual StatusCode finalize() final;

  /** @brief Create random CaloHits (gaussian distribution) for the vector of cells (aCells).
   * Vector of cells must contain all cells in the calorimeter with their cellIDs.
   */
  virtual void addRandomCellNoise(std::unordered_map<uint64_t, double>& aCells) final;
  /** @brief Remove cells with energy below threshold*sigma from the vector of cells
   */
  virtual void filterCellNoise(std::unordered_map<uint64_t, double>& aCells) final;

private:
  /// RMS of noise -- uniform RMS per cell in GeV
  Gaudi::Property<double> m_cellNoiseRMS{this, "cellNoiseRMS", 0.003, "uniform noise RMS per cell in GeV"};
  /// Offset of noise -- uniform offset per cell in GeV
  Gaudi::Property<double> m_cellNoiseOffset{this, "cellNoiseOffset", 0.0, "uniform noise offset per cell in GeV"};
  /// Energy threshold (Ecell < m_cellNoiseOffset + filterThreshold*m_cellNoiseRMS removed)
  Gaudi::Property<double> m_filterThreshold{
      this, "filterNoiseThreshold", 3,
      "remove cells with energy below offset + threshold * noise RMS"};
  /// Random Number Service
  SmartIF<IRndmGenSvc> m_randSvc;
  /// Gaussian random number generator used for smearing with a constant resolution (m_sigma)
  Rndm::Numbers m_gauss;
};

#endif /* RECCALORIMETER_NOISECALOCELLSFLATTOOL_H */
