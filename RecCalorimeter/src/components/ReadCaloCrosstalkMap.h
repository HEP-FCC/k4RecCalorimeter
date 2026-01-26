/**
 * @file RecCalorimeter/src/components/ReadCaloCrosstalkMap.h
 * @author Zhibo Wu, scott snyder <snyder@bnl.gov>
 * @date Updated Jul, 2026
 * @brief Return crosstalk information for a subdetector read from a TTree.
 */

#ifndef RECCALORIMETER_READCALOXTALKMAP_H
#define RECCALORIMETER_READCALOXTALKMAP_H

// from Gaudi
#include "GaudiKernel/AlgTool.h"

// Interface
#include "RecCaloCommon/ICaloCellConstantsSvc.h"
#include "RecCaloCommon/ICaloCellIndexerSvc.h"
#include "RecCaloCommon/ICaloReadCrosstalkMap.h"
#include <span>

class IGeoSvc;
class TFile;

/**
 * @brief Return crosstalk information for a subdetector read from a TTree.
 *
 * Reads a ROOT file containing the TTree with branches "cellId",
 * "list_crosstalk_neighbours" and "list_crosstalks".   The data are
 * read from the tree into a CrosstalkData object, which is stored
 * in the CaloCellConstantsSvc.   Allows returning the neighbours
 * and crosstalk coefficients for a given cell.
 *
 *  @author Zhibo Wu
 */
class ReadCaloCrosstalkMap : public extends<AlgTool, k4::recCalo::ICaloReadCrosstalkMap> {
public:
  using base_class::base_class;
  virtual ~ReadCaloCrosstalkMap() = default;

  /**
   * Standard Gaudi initialization method.
   */
  virtual StatusCode initialize() override final;

  /** Function to be called for the crosstalk neighbours of a cell.
   *   @param[in] aCellId, cellid of the cell of interest.
   *   @return span of cellIDs, corresponding to the crosstalk neighbours.
   */
  virtual std::span<const CellID> getNeighbours(CellID aCellId) const final override;

  /** Function to be called for the crosstalk coefficients between the input cell and its neighbouring cells.
   *   @param[in] aCellId, cellid of the cell of interest.
   *   @return span of crosstalk coefficients.
   */
  virtual std::span<const double> getCrosstalks(CellID aCellId) const final override;

private:
  /// Name of input root file that contains the TTree with
  /// cellID->vec<list_crosstalk_neighboursCellID> and
  /// cellId->vec<list_crosstalksCellID>
  Gaudi::Property<std::string> m_fileName{this, "fileName", "", "Name of the file that contains the crosstalk map."};

  /// Our subdetector ID.  If defaulted, ECAL_Barrel will be used.
  Gaudi::Property<int> m_detID{this, "detID", -1, "Subsystem ID for this detector.  Defaults to ECAL_Barrel."};

  /// Cell constants service.
  ServiceHandle<k4::recCalo::ICaloCellConstantsSvc> m_constantsSvc{
      this, "CaloCellConstantsSvc", "k4::recCalo::CaloCellConstantsSvc", "The cell constants service"};

  /// Cell indexing service.
  ServiceHandle<k4::recCalo::ICaloCellIndexerSvc> m_indexerSvc{
      this, "CaloCellIndexerSvc", "k4::recCalo::CaloCellIndexerSvc", "The cell indexing service."};

  /// Structure storing the crosstalk data.
  // We store both ID and coefficient data in flat vectors.
  struct CrosstalkData {
    /// This vector is indexed by the cell index.  It gives
    /// (INDEX,SIZE) pairs in m_neighbours giving the neighbour IDS
    /// for each cellID.
    std::vector<std::pair<size_t, size_t>> m_neighbourIndices;
    /// Neighbour cellID data.
    std::vector<CellID> m_neighbours;
    /// This vector is indexed by the cell index.  It gives
    /// (INDEX,SIZE) pairs in m_crosstalks giving the crosstalk coefficients
    /// for each cellID.
    std::vector<std::pair<size_t, size_t>> m_crosstalkIndices;
    /// Crosstalk coefficient data.
    std::vector<double> m_crosstalks;

    /// Given a cell index, return a span of neighbouring cell IDs.
    std::span<const CellID> getNeighbours(size_t cellNdx) const {
      const auto& indices = m_neighbourIndices.at(cellNdx);
      return std::span<const CellID>(m_neighbours.data() + indices.first, indices.second);
    }

    /// Given a cell index, return a span of crosstalk coefficients.
    std::span<const double> getCrosstalks(size_t cellNdx) const {
      const auto& indices = m_crosstalkIndices.at(cellNdx);
      return std::span<const double>(m_crosstalks.data() + indices.first, indices.second);
    }
  };

  /// Pointer to the crosstalk data.
  const CrosstalkData* m_data = nullptr;

  /// Pointer to the cell indexer.
  const k4::recCalo::ICaloIndexer* m_indexer = nullptr;

  /// Read crosstalk data from a file.
  CrosstalkData readData(TFile& xtalkFile) const;
};

#endif /* RECCALORIMETER_READCALOXTALKMAP_H */
