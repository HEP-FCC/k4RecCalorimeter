#include "CreateOpticalCaloCells.h"

// EDM4hep
#include "edm4hep/CalorimeterHit.h"

// k4FWCore
#include <k4FWCore/MetadataUtils.h>

// DD4hep
#include "DDSegmentation/BitFieldCoder.h"

// STL
#include <stdexcept>
#include <unordered_map>

DECLARE_COMPONENT(CreateOpticalCaloCells)

CreateOpticalCaloCells::CreateOpticalCaloCells(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {KeyValue("OpticalHits", "SimCalorimeterHitCollection"),
                        KeyValue("TruthHits", "SimCalorimeterHitCollection")},
                       {KeyValue("OutputCollection", "CalorimeterHitCollection"),
                        KeyValue("OutputLinks", "CaloHitSimCaloHitLinkCollection")}) {}

StatusCode CreateOpticalCaloCells::initialize() {
  info() << name() << ": calibration constant = " << m_calibConst.value() << endmsg;

  if (m_maskCherenkov.value()) {
    // Locate the 'cherenkov' field in the input cellID encoding to build its mask
    const std::string opticalKey = inputLocations("OpticalHits")[0];
    auto encoding = k4FWCore::getCellIDEncoding(opticalKey, this);
    if (!encoding.has_value()) {
      error() << "maskCherenkovForTruthLink is set but the cellID encoding of '" << opticalKey
              << "' could not be retrieved" << endmsg;
      return StatusCode::FAILURE;
    }
    dd4hep::DDSegmentation::BitFieldCoder decoder(encoding.value());
    try {
      m_cherenkovMask = decoder[m_cherenkovField.value()].mask();
    } catch (const std::exception& e) {
      error() << "cellID field '" << m_cherenkovField.value() << "' not found in encoding '" << encoding.value()
              << "': " << e.what() << endmsg;
      return StatusCode::FAILURE;
    }
    info() << "Cherenkov masking enabled: optical->truth match ignores field '" << m_cherenkovField.value()
           << "' (mask 0x" << std::hex << m_cherenkovMask << std::dec << ")" << endmsg;
  }
  return StatusCode::SUCCESS;
}

uint64_t CreateOpticalCaloCells::truthKey(uint64_t cellID) const {
  return m_maskCherenkov.value() ? (cellID & ~m_cherenkovMask) : cellID;
}

std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>
CreateOpticalCaloCells::operator()(const edm4hep::SimCalorimeterHitCollection& opticalHits,
                                   const edm4hep::SimCalorimeterHitCollection& truthHits) const {
  edm4hep::CalorimeterHitCollection outputHits;
  edm4hep::CaloHitSimCaloHitLinkCollection outputLinks;

  // Index the truth (energy-deposit) hits by their (optionally cherenkov-masked) cellID
  std::unordered_map<uint64_t, edm4hep::SimCalorimeterHit> truthHitMap;
  truthHitMap.reserve(truthHits.size());
  for (const auto& truthHit : truthHits) {
    truthHitMap[truthKey(truthHit.getCellID())] = truthHit;
  }

  // Digitize the optical hits and link each cell to its truth deposit
  for (const auto& opticalHit : opticalHits) {
    if (opticalHit.getEnergy() <= 0) {
      continue; // skip cells with no signal
    }
    auto caloHit = outputHits.create();
    caloHit.setCellID(opticalHit.getCellID());
    caloHit.setPosition(opticalHit.getPosition());
    caloHit.setEnergy(opticalHit.getEnergy() * m_calibConst.value());

    auto truthIt = truthHitMap.find(truthKey(opticalHit.getCellID()));
    if (truthIt != truthHitMap.end()) {
      auto hitLink = outputLinks.create();
      hitLink.setFrom(caloHit);
      hitLink.setTo(truthIt->second);
    }
  }

  debug() << "Digitized " << opticalHits.size() << " optical hits -> " << outputHits.size() << " cells with "
          << outputLinks.size() << " truth links" << endmsg;

  return std::make_tuple(std::move(outputHits), std::move(outputLinks));
}
