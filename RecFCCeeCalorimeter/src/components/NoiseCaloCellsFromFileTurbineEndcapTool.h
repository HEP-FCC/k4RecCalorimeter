#ifndef RECFCCEECALORIMETER_NOISECALOCELLSFROMFILETURBINEENDCAPTOOL_H
#define RECFCCEECALORIMETER_NOISECALOCELLSFROMFILETURBINEENDCAPTOOL_H

#include "NoiseCaloCellsFromFileBaseTool.h"
#include "RecCaloCommon/ICellPositionsTool.h"

#include "GaudiKernel/ToolHandle.h"

/** @class NoiseCaloCellsFromFileTurbineEndcapTool
 *
 *  Tool for calorimeter noise - in endcap, with noise histograms per wheel,
 *  implemented as 2D hists vs rho, z.
 *  Inherits from common base tool and defines how to retrieve proper bin in
 *  noise histograms for cell with given cellID
 *
 *  @author Erich Varnes
 *  @author Giovanni Marchiori
 *  @date   2026-06
 *
 */

class NoiseCaloCellsFromFileTurbineEndcapTool : public NoiseCaloCellsFromFileBaseTool {
public:
  NoiseCaloCellsFromFileTurbineEndcapTool(const std::string& type, const std::string& name, const IInterface* parent)
      : NoiseCaloCellsFromFileBaseTool(type, name, parent) {
    // Override some property defaults from the base class.
    m_readoutName = "ECalEndcapTurbine";
  }

protected:
  virtual StatusCode initBinning(NoiseData& data, const k4::recCalo::ICaloIndexer& indexer) const override;

private:
  unsigned getBin(const char* what, const TH1& h, unsigned iRho, unsigned iZ) const;

  /// Unused, but temporarily here for config compatibility
  ToolHandle<k4::recCalo::ICellPositionsTool> m_cellPositionsTool{this, "cellPositionsTool", "",
                                                                  "Handle for tool to retrieve cell positions"};
};

#endif /* RECFCCEECALORIMETER_NOISECALOCELLFROMFILETURBINEENDCAPTOOL_H */
