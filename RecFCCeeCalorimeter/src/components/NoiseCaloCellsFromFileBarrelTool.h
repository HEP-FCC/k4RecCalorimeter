#ifndef RECFCCEECALORIMETER_NOISECALOCELLSFROMFILEBARRELTOOL_H
#define RECFCCEECALORIMETER_NOISECALOCELLSFROMFILEBARRELTOOL_H

#include "NoiseCaloCellsFromFileBaseTool.h"

/** @class NoiseCaloCellsFromFileBarrelTool
 *
 *  Tool for calorimeter noise - in barrel, with noise histograms per layer,
 *  implemented as 1D hists vs theta.
 *  Inherits from common base tool and defines how to retrieve proper bin in
 *  noise histograms for cell with given cellID
 *
 *  @author Giovanni Marchiori
 *  @date   2026-06
 *
 */

class NoiseCaloCellsFromFileBarrelTool : public NoiseCaloCellsFromFileBaseTool {
public:
  using NoiseCaloCellsFromFileBaseTool::NoiseCaloCellsFromFileBaseTool;

private:
  /// get bin in histogram for given cellID
  unsigned getBin(k4::recCalo::INoiseCaloCellsTool::CellID aCellId) const override final;
};

#endif /* RECFCCEECALORIMETER_NOISECALOCELLFROMFILEBARRELTOOL_H */
