#ifndef RECCALORIMETER_TOPOCALONOISYCELLS_H
#define RECCALORIMETER_TOPOCALONOISYCELLS_H

// from Gaudi
#include "GaudiKernel/AlgTool.h"

// Interfaces
#include "RecCaloCommon/ICaloCellConstantsSvc.h"
#include "RecCaloCommon/ICaloCellIndexerSvc.h"
#include "RecCaloCommon/INoiseConstTool.h"
#include <memory>
#include <utility>

// DD4HEP
namespace dd4hep::DDSegmentation {
class BitFieldCoder;
} // namespace dd4hep::DDSegmentation

class IGeoSvc;
class TFile;

/** @class TopoCaloNoisyCells Reconstruction/RecCalorimeter/src/components/TopoCaloNoisyCells.h
 *TopoCaloNoisyCells.h
 *
 *  Tool that reads a ROOT file containing the TTree with branchs "cellId", "noiseLevel", and "noiseOffset".
 *  This tool reads the tree, creates a map, and allows a lookup of noise level and mean noise of a cell, by its cellID.
 *
 *  @author Coralie Neubueser
 */

class TopoCaloNoisyCells : public extends<AlgTool, k4::recCalo::INoiseConstTool> {
public:
  using base_class::base_class;
  virtual ~TopoCaloNoisyCells() = default;
  /** Read a root file and the stored TTree of cellIDs to noise values.
   * return StatusCode
   */
  virtual StatusCode initialize() override final;

  /** Expected noise per cell in terms of sigma of Gaussian distibution.
   *   @param[in] aCellId of the cell of interest.
   *   return double.
   */
  virtual double getNoiseRMSPerCell(CellID aCellId) const override final;

  /** Expected noise per cell in terms of mean of distibution.
   *   @param[in] aCellId of the cell of interest.
   *   return double.
   */
  virtual double getNoiseOffsetPerCell(CellID aCellId) const override final;

  /** Expected noise per cell.
   *   @param[in] aCellId of the cell of interest.
   *   return [rms, offset]
   */
  virtual std::pair<double, double> getNoisePerCell(CellID aCellId) const override final;

private:
  /// Name
  Gaudi::Property<std::string> m_fileName{this, "fileName",
                                          "/afs/cern.ch/user/c/cneubuse/public/FCChh/cellNoise_map_segHcal.root"};
  // clang-format off
  ServiceHandle<k4::recCalo::ICaloCellConstantsSvc> m_constantsSvc
  { this, "CaloCellConstantsSvc", "k4::recCalo::CaloCellConstantsSvc", "" };
  ServiceHandle<k4::recCalo::ICaloCellIndexerSvc> m_indexerSvc
  { this, "CaloCellIndexerSvc", "k4::recCalo::CaloCellIndexerSvc", "" };
  // clang-format on

  /// System encoding string
  Gaudi::Property<std::string> m_systemEncoding{this, "systemEncoding", "system:4", "System encoding string"};

  struct NoiseData {
    // rms, offset
    std::vector<std::pair<double, double>> m_noise;
    const k4::recCalo::ICaloIndexer* m_indexer;
  };
  const NoiseData* m_data = nullptr;

  const k4::recCalo::ICaloIndexer* m_indexer = nullptr;
  std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_decoder;
  int m_indexSystem = -1;

  NoiseData readData(TFile& inFile) const;
};

#endif /* RECCALORIMETER_TOPOCALONOISYCELLS_H */
