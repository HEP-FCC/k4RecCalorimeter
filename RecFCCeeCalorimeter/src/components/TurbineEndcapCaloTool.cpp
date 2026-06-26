/**
 * @file RecFCCeeCalorimeter/src/components/TurbineEndcapCaloTool.cpp
 * @author Giovanni Marchiori <giovanni.marchiori@cern.ch>
 * @date June, 2026
 * @brief Calorimeter tool for ALLEGRO ECal endcap
 */

#include "TurbineEndcapCaloTool.h"
#include "RecCaloCommon/IDMapIndexer.h"
#include "detectorSegmentations/FCCSWEndcapTurbine_k4geo.h"
#include "k4FWCore/GaudiChecks.h"

DECLARE_COMPONENT(TurbineEndcapCaloTool)

/** Fill vector with all existing cells for this geometry.
 */
StatusCode TurbineEndcapCaloTool::collectCells(std::vector<uint64_t>& cells) const {

  const auto* seg =
      dynamic_cast<const dd4hep::DDSegmentation::FCCSWEndcapTurbine_k4geo*>(readout().segmentation().segmentation());
  if (!seg) {
    error() << "Unable to cast segmentation pointer!!!! Tool only applicable to FCCSWEndcapTurbine_k4geo "
               "segmentation."
            << endmsg;
    return StatusCode::FAILURE;
  }

  const dd4hep::DDSegmentation::BitFieldCoder& decoder = *readout().idSpec().decoder();

  dd4hep::DDSegmentation::CellID cID = 0;
  decoder.set(cID, "system", this->id());
  decoder.set(cID, "cryo", 0);
  decoder.set(cID, "type", 0);
  decoder.set(cID, "subtype", 0);

  size_t side_id = decoder.index(seg->fieldNameSide());
  size_t wheel_id = decoder.index(seg->fieldNameWheel());
  size_t module_id = decoder.index(seg->fieldNameModule());
  size_t layer_id = decoder.index(seg->fieldNameLayer());
  size_t rho_id = decoder.index(seg->fieldNameRho());
  size_t z_id = decoder.index(seg->fieldNameZ());

  const int nsides = 2;
  int sides[nsides] = {-1, 1};
  const int nWheels = seg->numWheels();
  int ncells = 0;
  for (int iWheel = 0; iWheel < nWheels; ++iWheel) {
    ncells += seg->nModules(iWheel) * seg->numCellsRho(iWheel) * seg->numCellsZ(iWheel);
  }
  ncells *= nsides;
  cells.reserve(ncells);
  for (int iside : sides) {
    decoder.set(cID, side_id, iside);
    for (int iWheel = 0; iWheel < nWheels; ++iWheel) {
      decoder.set(cID, wheel_id, iWheel);
      for (int iModule = 0; iModule < seg->nModules(iWheel); iModule++) {
        decoder.set(cID, module_id, iModule);
        for (int iRho = 0; iRho < seg->numCellsRho(iWheel); iRho++) {
          decoder.set(cID, rho_id, iRho);
          for (int iZ = 0; iZ < seg->numCellsZ(iWheel); iZ++) {
            decoder.set(cID, z_id, iZ);
            int iLayer = seg->expLayer(iWheel, iRho, iZ);
            decoder.set(cID, layer_id, iLayer);
            cells.push_back(cID);
          }
        }
      }
    }
  }

  return StatusCode::SUCCESS;
}

/** Return a new indexer object for this subdetector.
 */
std::unique_ptr<k4::recCalo::ICaloIndexer> TurbineEndcapCaloTool::indexer() const {
  const auto* seg =
      dynamic_cast<const dd4hep::DDSegmentation::FCCSWEndcapTurbine_k4geo*>(readout().segmentation().segmentation());
  if (!seg) {
    error() << "Unable to cast segmentation pointer!!!! Tool only applicable to FCCSWEndcapTurbine_k4geo "
               "segmentation."
            << endmsg;
    return nullptr;
  }

  using Indexer_t = k4::recCalo::IDMapIndexer<5>; // 5 fields to use for indexing (including ignored ones): side, wheel,
                                                  // rho, z, module
  dd4hep::IDDescriptor idSpec = readout().idSpec();
  std::vector<Indexer_t::FieldDesc_t> fields{Indexer_t::IDMap_t::makeDesc(*idSpec.field(seg->fieldNameSide())),
                                             Indexer_t::IDMap_t::makeDesc(*idSpec.field(seg->fieldNameWheel())),
                                             Indexer_t::IDMap_t::makeDesc(*idSpec.field(seg->fieldNameRho())),
                                             Indexer_t::IDMap_t::makeDesc(*idSpec.field(seg->fieldNameZ())),
                                             Indexer_t::IDMap_t::makeDesc(*idSpec.field(seg->fieldNameModule()))};

  const dd4hep::BitFieldElement& sysField = *idSpec.field("system");
  if (sysField.offset() != 0) {
    throw std::runtime_error("Bad system field offset; must be zero");
  }

  // ignore layer field (is function of other fields)
  std::vector<Indexer_t::FieldDesc_t> ignoredFields{Indexer_t::IDMap_t::makeDesc(*idSpec.field(seg->fieldNameLayer()))};

  return std::make_unique<Indexer_t>(this->id(), sysField.width(), fields, cellIDs(), 1300000, ignoredFields);
}
