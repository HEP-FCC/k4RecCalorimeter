#ifndef RECFCCEECALORIMETER_NOISECALOCELLSFROMFILETURBINEENDCAPTOOL_H
#define RECFCCEECALORIMETER_NOISECALOCELLSFROMFILETURBINEENDCAPTOOL_H

#include "NoiseCaloCellsFromFileBaseTool.h"

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
  using NoiseCaloCellsFromFileBaseTool::NoiseCaloCellsFromFileBaseTool;

private:
  /// get bin in histogram for given cellID
  unsigned getBin(CellID aCellId) const override final;
};

#endif /* RECFCCEECALORIMETER_NOISECALOCELLFROMFILETURBINEENDCAPTOOL_H */
