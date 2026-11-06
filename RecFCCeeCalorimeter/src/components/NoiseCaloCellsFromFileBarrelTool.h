#ifndef RECFCCEECALORIMETER_NOISECALOCELLSFROMFILEBARRELTOOL_H
#define RECFCCEECALORIMETER_NOISECALOCELLSFROMFILEBARRELTOOL_H

#include "NoiseCaloCellsFromFileBaseTool.h"
#include "RecCaloCommon/ICellPositionsTool.h"

#include "GaudiKernel/ToolHandle.h"

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
  NoiseCaloCellsFromFileBarrelTool(const std::string& type, const std::string& name, const IInterface* parent)
      : NoiseCaloCellsFromFileBaseTool(type, name, parent) {
    // Override some property defaults from the base class.
    m_readoutName = "ECalBarrelThetaModuleMerged";
  }

  virtual StatusCode initialize() override;

protected:
  virtual StatusCode initBinning(NoiseData& data, const k4::recCalo::ICaloIndexer& indexer) const override;

private:
  /// Handle for tool to get cell positions
  ToolHandle<k4::recCalo::ICellPositionsTool> m_cellPositionsTool{this, "cellPositionsTool", "",
                                                                  "Handle for tool to retrieve cell positions"};
};

#endif /* RECFCCEECALORIMETER_NOISECALOCELLFROMFILEBARRELTOOL_H */
