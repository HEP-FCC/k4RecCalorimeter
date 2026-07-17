#ifndef RECCALORIMETER_READNOISEFROMFILETOOL_H
#define RECCALORIMETER_READNOISEFROMFILETOOL_H

#include "RecCaloCommon/ReadNoiseFromFileBaseTool.h"

/** @class ReadNoiseFromFileTool
 *
 *  Tool to read the stored noise constant per cell in the calorimeters
 *  Access noise constants from TH1F histogram (noise vs. |eta|)
 *
 *  @author Jana Faltova, Coralie Neubueser
 *  @date   2018-01
 *
 */

class ReadNoiseFromFileTool : public ReadNoiseFromFileBaseTool {
public:
  using ReadNoiseFromFileBaseTool::ReadNoiseFromFileBaseTool;
  virtual ~ReadNoiseFromFileTool() = default;

protected:
  virtual StatusCode initBinning(NoiseData& data, const k4::recCalo::ICaloIndexer& indexer) const;
};

#endif /* RECCALORIMETER_READNOISEFROMFILETOOL_H */
