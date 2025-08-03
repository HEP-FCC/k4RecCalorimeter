#include "DualCrysCalDigi.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/MutableCalorimeterHit.h"
#include "edm4hep/MutableCaloHitContribution.h"
#include "edm4hep/MutableCaloHitSimCaloHitLink.h"

#include <unordered_map>

DECLARE_COMPONENT(DualCrysCalDigi)

DualCrysCalDigi::DualCrysCalDigi(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc) {}

StatusCode DualCrysCalDigi::initialize() {
  if (Service("GeoSvc", m_geoSvc).isFailure()) {
    error() << "Unable to locate Geometry Service." << endmsg;
    return StatusCode::FAILURE;
  }

  if (Service("UniqueIDGenSvc", m_uidSvc).isFailure()) {
    error() << "Unable to locate UniqueIDGenSvc." << endmsg;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>
DualCrysCalDigi::operator()(const edm4hep::SimCalorimeterHitCollection& simHits,
                            const edm4hep::EventHeaderCollection& /*headers*/) const {
  using namespace edm4hep;

  CalorimeterHitCollection outHits;
  CaloHitSimCaloHitLinkCollection outLinks;

  std::unordered_map<uint64_t, MutableCalorimeterHit> hitMap;

  for (const auto& simHit : simHits) {
    float totalEdep = 0.0;
    for (unsigned i = 0; i < simHit.getContributions().size(); ++i) {
      totalEdep += simHit.getEdep(i);
    }

    if (totalEdep < m_thresholdCal.value()) continue;
    if (totalEdep > m_maxHitEnergyCal.value()) continue;

    uint64_t cellID = simHit.getCellID();

    auto it = hitMap.find(cellID);
    if (it == hitMap.end()) {
      MutableCalorimeterHit newHit;
      newHit.setCellID(cellID);
      newHit.setEnergy(totalEdep * m_calibrCoeffCal.value());
      newHit.setType(static_cast<int>(calorimeter::HitType::Cluster));  // or Raw
      newHit.setTime(simHit.getTime());
      newHit.setPosition({simHit.getPosition().x, simHit.getPosition().y, simHit.getPosition().z});
      outHits.push_back(newHit);
      hitMap[cellID] = newHit;
    } else {
      it->second.setEnergy(it->second.getEnergy() + totalEdep * m_calibrCoeffCal.value());
    }

    // Create Sim -> Digi link
    MutableCaloHitSimCaloHitLink link;
    link.setSim(hitMap[cellID]);
    link.setContributions(simHit.getContributions());
    link.setWeight(1.0);  // or fraction of total
    outLinks.push_back(link);
  }

  return {outHits, outLinks};
}
