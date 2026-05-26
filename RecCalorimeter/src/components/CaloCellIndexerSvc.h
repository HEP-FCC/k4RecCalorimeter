// This file's extension implies that it's C, but it's really -*- C++ -*-.
/**
 * @file RecCalorimeter/src/components/CaloCellIndexerSvc.h
 * @author scott snyder <snyder@bnl.gov>
 * @date Jan, 2026
 * @brief Holder for calorimeter cell indexers.
 */

#ifndef RECCALORIMETER_CALOCELLINDEXERSVC_H
#define RECCALORIMETER_CALOCELLINDEXERSVC_H

// clang-format off: headers best listed in reverse-dependency order.
#include "RecCaloCommon/ICaloCellIndexerSvc.h"
#include "RecCaloCommon/ICaloCellConstantsSvc.h"
#include "RecCaloCommon/ICalorimeterTool.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"
#include <bitset>
#include <climits>
#include <mutex>
#include <unordered_map>
// clang-format on

namespace k4::recCalo {

/**
 * @brief Holder for calorimeter cell indexers.
 *
 * This holds the indexer objects for different subdetectors.
 * We put it in a service to avoid duplicating them across different tools.
 *
 * One might think of having this be part of ICaloCellConstantsSvc,
 * but we run into initialization loops in that case.
 */
class CaloCellIndexerSvc : public extends<Service, ICaloCellIndexerSvc> {
public:
  using base_class::base_class;

  /**
   * @brief Gaudi initialize method.
   */
  virtual StatusCode initialize() override;

  /**
   * @brief Return indexer for a given subdetector.
   * @param detID Subdetector ID for the desired indexer.
   * @param quiet If true, don't print an error if we don't find an indexer.
   *
   * Returns a pointer to the indexer or nullptr if there isn't one defined.
   */
  const ICaloIndexer* indexer(int detID, bool quiet = false) const override;

  /**
   * @brief Return indexer for a given set of subdetectors.
   * @param detIDs Subdetector IDs to index.
   * @param quiet If true, don't print an error if we don't find an indexer.
   *
   * Returns a pointer to the indexer or nullptr if there isn't one defined.
   */
  virtual const ICaloIndexer* indexer(std::span<const int> detIDs, bool quiet = false) override;

private:
  /// List of calorimeter tools that can provide indexers.
  ToolHandleArray<ICalorimeterTool> m_geoTools{this, "GeoTools", {}};

  /// The cell constants service.
  // clang-format off: this gets broken illogically.
  ServiceHandle<k4::recCalo::ICaloCellConstantsSvc> m_constantsSvc
    { this, "CaloCellConstantsSvc", "k4::recCalo::CaloCellConstantsSvc" };
  // clang-format on

  /// Map from bit mask of detector IDs to indexer objects.
  constexpr static int MAX_DETID = 255;
  using mask_t = std::bitset<MAX_DETID + 1>;
  std::unordered_map<mask_t, std::unique_ptr<ICaloIndexer>> m_indexers;

  /// Guard access to the map.
  mutable std::recursive_mutex m_mutex;
};

} // namespace k4::recCalo

#endif // not RECCALORIMETER_CALOCELLINDEXERSVC_H
