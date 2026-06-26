// This file's extension implies that it's C, but it's really -*- C++ -*-.
/**
 * @file RecFCCeeCalorimeter/src/components/TurbineEndcapCaloTool.h
 * @author Giovanni Marchiori <giovanni.marchiori@cern.ch>
 * @date June, 2026
 * @brief Calorimeter tool for ALLEGRO ECal endcap
 */

#ifndef RECFCCEECALORIMETER_TURBINEENDCAPCALOTOOL_H
#define RECFCCEECALORIMETER_TURBINEENDCAPCALOTOOL_H

#include "RecCaloCommon/CalorimeterToolBase.h"

/** @class TurbineEndcapCaloTool
 *
 *  Manage IDs for ALLEGRO ECal turbine endcap
 */
class TurbineEndcapCaloTool : public CalorimeterToolBase {
public:
  using CalorimeterToolBase::CalorimeterToolBase;
  virtual ~TurbineEndcapCaloTool() = default;

  /** Return a new indexer object for this subdetector.
   */
  virtual std::unique_ptr<k4::recCalo::ICaloIndexer> indexer() const override final;

protected:
  /** Fill vector with all existing cells for this geometry.
   */
  virtual StatusCode collectCells(std::vector<uint64_t>& cells) const override final;
};

#endif /* RECFCCEECALORIMETER_TURBINEENDCAPCALOTOOL_H */
