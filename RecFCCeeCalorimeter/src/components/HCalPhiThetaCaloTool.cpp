/**
 * @file RecFCCeeCalorimeter/src/components/HCalPhiThetaCaloTool.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date Feb, 2026
 * @brief Calorimeter tool for Allegro HCal.
 */

#include "HCalPhiThetaCaloTool.h"
#include "RecCaloCommon/IDMapIndexer.h"
#include "detectorSegmentations/FCCSWHCalPhiTheta_k4geo.h"
#include "k4FWCore/GaudiChecks.h"

DECLARE_COMPONENT(HCalPhiThetaCaloTool)

/** Return the name of this subdetector, to be used to find the
    subdetector ID.
*/
std::string HCalPhiThetaCaloTool::detectorName() const {
  if (readoutName() == "HCalBarrelReadout")
    return "HCAL_Barrel";
  else
    return "HCAL_Endcap";
}

/** Fill vector with all existing cells for this geometry.
 */
StatusCode HCalPhiThetaCaloTool::collectCells(std::vector<uint64_t>& cells) const {
  cells.reserve(readoutName() == "HCalBarrelReadout" ? 210944 : 80896);
  const auto* seg =
      dynamic_cast<const dd4hep::DDSegmentation::FCCSWHCalPhiTheta_k4geo*>(readout().segmentation().segmentation());
  if (!seg) {
    error() << "Unable to cast segmentation pointer!!!! Tool only applicable to FCCSWHCalPhiTheta_k4geo "
               "segmentation."
            << endmsg;
    return StatusCode::FAILURE;
  }

  const dd4hep::DDSegmentation::BitFieldCoder& decoder = *readout().idSpec().decoder();

  dd4hep::DDSegmentation::CellID cID = 0;
  decoder.set(cID, "system", this->id());

  size_t layer_id = decoder.index(seg->fieldNameLayer());
  size_t theta_id = decoder.index(seg->fieldNameTheta());
  size_t phi_id = decoder.index(seg->fieldNamePhi());

  int numLayers = 0;
  for (int n : seg->numLayers())
    numLayers += n;
  for (int layer = 0; layer < numLayers; ++layer) {
    const std::vector<int>& thetaBins = seg->thetaBins(layer);
    decoder.set(cID, layer_id, layer);
    for (int theta : thetaBins) {
      decoder.set(cID, theta_id, theta);
      for (int phi = 0; phi < seg->phiBins(); ++phi) {
        decoder.set(cID, phi_id, phi);
        cells.push_back(cID);
      }
    }
  }

  return StatusCode::SUCCESS;
}

/** Return a new indexer object for this subdetector.
 */
std::unique_ptr<k4::recCalo::ICaloIndexer> HCalPhiThetaCaloTool::indexer() const {
  const auto* seg =
      dynamic_cast<const dd4hep::DDSegmentation::FCCSWHCalPhiTheta_k4geo*>(readout().segmentation().segmentation());
  if (!seg) {
    error() << "Unable to cast segmentation pointer!!!! Tool only applicable to FCCSWHCalPhiTheta_k4geo "
               "segmentation."
            << endmsg;
    return nullptr;
  }

  using Indexer_t = k4::recCalo::IDMapIndexer<3>;
  dd4hep::IDDescriptor idSpec = readout().idSpec();
  // Sorry, clang-format insists on making this hard to read.
  std::vector<Indexer_t::FieldDesc_t> fields{Indexer_t::IDMap_t::makeDesc(*idSpec.field(seg->fieldNameLayer())),
                                             Indexer_t::IDMap_t::makeDesc(*idSpec.field(seg->fieldNameTheta())),
                                             Indexer_t::IDMap_t::makeDesc(*idSpec.field(seg->fieldNamePhi()))};

  const dd4hep::BitFieldElement& sysField = *idSpec.field("system");
  if (sysField.offset() != 0) {
    throw std::runtime_error("Bad system field offset; must be zero");
  }

  // Again, clang-format...
  return std::make_unique<Indexer_t>(this->id(), sysField.width(), fields, cellIDs(),
                                     readoutName() == "HCalBarrelReadout" ? 900000 : 700000);
}
