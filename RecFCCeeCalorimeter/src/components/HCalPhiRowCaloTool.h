// This file's extension implies that it's C, but it's really -*- C++ -*-.
/**
 * @file RecFCCeeCalorimeter/src/components/HCalPhiRowCaloTool.h
 * @author scott snyder <snyder@bnl.gov>
 * @date May, 2026
 * @brief Calorimeter tool for Allegro HCal (with row indexing).
 */

#ifndef RECFCCEECALORIMETER_HCALPHIROWCALOTOOL_H
#define RECFCCEECALORIMETER_HCALPHIROWCALOTOOL_H

#include "RecCaloCommon/CalorimeterToolBase.h"

/** @class HCalPhiRowCaloTool
 *
 *  Manage IDs for HCal (with row indexing).
 */
class HCalPhiRowCaloTool : public CalorimeterToolBase {
public:
  using CalorimeterToolBase::CalorimeterToolBase;
  virtual ~HCalPhiRowCaloTool() = default;

  /** Return the name of this subdetector, to be used to find the
      subdetector ID.
   */
  virtual std::string detectorName() const override;

  /** Return a new indexer object for this subdetector.
   */
  virtual std::unique_ptr<k4::recCalo::ICaloIndexer> indexer() const override final;

protected:
  /** Fill vector with all existing cells for this geometry.
   */
  virtual StatusCode collectCells(std::vector<uint64_t>& cells) const override final;
};

#endif /* RECFCCEECALORIMETER_HCALPHIROWCALOTOOL_H */
