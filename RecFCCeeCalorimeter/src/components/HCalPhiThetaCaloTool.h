// This file's extension implies that it's C, but it's really -*- C++ -*-.
/**
 * @file RecFCCeeCalorimeter/src/components/HCalPhiThetaCaloTool.h
 * @author scott snyder <snyder@bnl.gov>
 * @date Feb, 2026
 * @brief Calorimeter tool for Allegro HCal.
 */

#ifndef RECFCCEECALORIMETER_HCALPHITHETACALOTOOL_H
#define RECFCCEECALORIMETER_HCALPHITHETACALOTOOL_H

#include "RecCaloCommon/CalorimeterToolBase.h"

/** @class HCalPhiThetaCaloTool
 *
 *  Manage IDs for HCal.
 */
class HCalPhiThetaCaloTool : public CalorimeterToolBase {
public:
  using CalorimeterToolBase::CalorimeterToolBase;
  virtual ~HCalPhiThetaCaloTool() = default;

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

#endif /* RECFCCEECALORIMETER_HCALPHITHETACALOTOOL_H */
