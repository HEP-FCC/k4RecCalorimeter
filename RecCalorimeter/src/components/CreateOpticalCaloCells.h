#ifndef RECCALORIMETER_CREATEOPTICALCALOCELLS_H
#define RECCALORIMETER_CREATEOPTICALCALOCELLS_H

// k4FWCore
#include "k4FWCore/Transformer.h"

// EDM4hep
#include "edm4hep/CaloHitSimCaloHitLinkCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

#include <cstdint>
#include <string>
#include <tuple>

/** @class CreateOpticalCaloCells
 *
 *  Digitizer for optical-readout calorimeters. Turns the optical-photon SimCalorimeterHits into
 *  CalorimeterHits (energy = signal * calibrationConstant, position taken from
 *  the SimHit) and links each cell to the truth (energy-deposit)
 *  SimCalorimeterHit at the same cellID.
 *
 *  Inputs:
 *   - OpticalHits : optical-photon SimCalorimeterHits, provide the cell signal
 *   - TruthHits   : energy-deposit SimCalorimeterHits, the link targets
 *
 *  The optical->truth match uses the full cellID by default. When the scint and
 *  Cherenkov signals of one physical cell share a cellID that differs only in
 *  the 'cherenkov' field - as in the crystal ECAL, where both readings come from
 *  the same energy deposit - set maskCherenkovForTruthLink=true to drop that
 *  field from the match, so both readings link to the same truth hit. In the
 *  HCAL the scint and Cherenkov fibres are distinct cells, so the
 *  full-cellID default already links them to their own deposits.
 */

class CreateOpticalCaloCells final
    : public k4FWCore::MultiTransformer<std::tuple<edm4hep::CalorimeterHitCollection,
                                                   edm4hep::CaloHitSimCaloHitLinkCollection>(
          const edm4hep::SimCalorimeterHitCollection&, const edm4hep::SimCalorimeterHitCollection&)> {
public:
  CreateOpticalCaloCells(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;

  std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>
  operator()(const edm4hep::SimCalorimeterHitCollection& opticalHits,
             const edm4hep::SimCalorimeterHitCollection& truthHits) const override;

private:
  /// cellID key used to match an optical hit to a truth hit (optionally cherenkov-masked)
  uint64_t truthKey(uint64_t cellID) const;

  Gaudi::Property<float> m_calibConst{this, "calibrationConstant", 1.,
                                      "Multiplicative calibration from optical signal to cell energy"};
  Gaudi::Property<bool> m_maskCherenkov{this, "maskCherenkovForTruthLink", false,
                                        "Ignore the 'cherenkov' cellID field when matching optical hits to truth hits"};
  Gaudi::Property<std::string> m_cherenkovField{
      this, "cherenkovFieldName", "cherenkov",
      "Name of the cellID field flagging the Cherenkov readout (used only when maskCherenkovForTruthLink=true)"};

  /// mask of the cherenkov field, filled in initialize() when masking is enabled
  uint64_t m_cherenkovMask{0};
};

#endif /* RECCALORIMETER_CREATEOPTICALCALOCELLS_H */
