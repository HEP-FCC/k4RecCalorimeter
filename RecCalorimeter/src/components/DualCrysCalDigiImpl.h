#include "DualCrysCalDigi.h"

#include "DDSegmentation/Segmentation.h"
#include "DD4hep/DD4hepUnits.h"

#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CaloHitContribution.h"
#include "edm4hep/MutableCalorimeterHit.h"
#include "edm4hep/MutableCaloHitContribution.h"
#include "edm4hep/MutableCaloHitSimCaloHitLink.h"

#include <map>
#include <cmath>

DualCrysCalDigi::DualCrysCalDigi(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc) {}

StatusCode DualCrysCalDigi::initialize() {
  if (GaudiAlgorithm::initialize().isFailure())
    return StatusCode::FAILURE;

  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "GeoSvc not found." << endmsg;
    return StatusCode::FAILURE;
  }

  m_uidSvc = service("UniqueIDGenSvc");
  if (!m_uidSvc) {
    error() << "UniqueIDGenSvc not found." << endmsg;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>
DualCrysCalDigi::operator()(const edm4hep::SimCalorimeterHitCollection& simHits,
                            const edm4hep::EventHeaderCollection&) const {

  edm4hep::CalorimeterHitCollection caloOut;
  edm4hep::CaloHitSimCaloHitLinkCollection links;

  std::unordered_map<uint64_t, std::vector<const edm4hep::SimCalorimeterHit*>> cellHitMap;

  // Group hits by cellID
  for (auto& simHit : simHits) {
    if (simHit.getEnergy() > m_maxHitEnergyCal)
      continue;
    cellHitMap[simHit.getCellID()].push_back(&simHit);
  }

  for (const auto& [cellID, hits] : cellHitMap) {
    double totalEnergy = 0.0;
    double timeWeightedSum = 0.0;
    double maxContribution = 0.0;

    for (const auto* simHit : hits) {
      for (auto contrib : simHit->getContributions()) {
        totalEnergy += contrib.getEnergy();
        timeWeightedSum += contrib.getTime() * contrib.getEnergy();
        if (contrib.getEnergy() > maxContribution)
          maxContribution = contrib.getEnergy();
      }
    }

    if (totalEnergy < m_thresholdCal)
      continue;

    edm4hep::MutableCalorimeterHit digiHit;
    digiHit.setCellID(cellID);
    digiHit.setEnergy(totalEnergy * m_calibrCoeffCal);
    digiHit.setTime(totalEnergy > 0 ? timeWeightedSum / totalEnergy : 0);
    digiHit.setType(CaloHitType::ECAL);
    caloOut.push_back(digiHit);

    for (const auto* simHit : hits) {
      edm4hep::MutableCaloHitSimCaloHitLink link;
      link.setRec(digiHit);
      link.setSim(*simHit);
      link.setWeight(simHit->getEnergy() / totalEnergy);
      links.push_back(link);
    }
  }

  return {caloOut, links};
}
