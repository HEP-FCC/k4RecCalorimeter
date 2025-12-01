/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <cctype>
#include <cstdlib>  // abs
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"
#include "GaudiKernel/MsgStream.h"
#include "edm4hep/CalorimeterHit.h"
#include "DualCrysCalorimeterHit.h"
#include "edm4hep/Constants.h"
#include "DualCrysCalDigi.h"

DECLARE_COMPONENT(DualCrysCalDigi)

DualCrysCalDigi::DualCrysCalDigi(const std::string& aName, ISvcLocator* aSvcLoc)
    : MultiTransformer(aName, aSvcLoc,
                       {
                           KeyValues("CALCollection", {"ECalEcalCollection"}),
                           KeyValues("HeaderName", {"EventHeader"}),
                       },
                       {KeyValues("CALOutputCollections", {"CalorimeterHit"}),
                        KeyValues("RelationOutputCollection", {"RelationCalHit"})}) {
  m_uidSvc = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
  if (!m_uidSvc) {
    error() << "Unable to get UniqueIDGenSvc" << endmsg;
  }

  m_geoSvc = serviceLocator()->service("GeoSvc");  // important to initialize m_geoSvc
}

StatusCode DualCrysCalDigi::initialize() {
  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection> DualCrysCalDigi::operator()(
    const edm4hep::SimCalorimeterHitCollection& SimCaloHits, const edm4hep::EventHeaderCollection& headers) const {
  debug() << " process event : " << headers[0].getEventNumber() << " - run  " << headers[0].getRunNumber()
          << endmsg;  // headers[0].getRunNumber(),headers[0].getEventNumber()

  auto calcol    = edm4hep::CalorimeterHitCollection();
  auto calRelcol = edm4hep::CaloHitSimCaloHitLinkCollection();

  std::string initString;

  std::string colName    = m_calCollections;
  CHT::Layout caloLayout = layoutFromString(colName);


  initString = m_geoSvc->constantAsString(m_encodingStringVariable.value());
  dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(initString);  // check!


  for (const auto& hit : SimCaloHits) {
    const int cellID = hit.getCellID();
    float     energy = hit.getEnergy();
    //Get the layer number
    unsigned int layer = bitFieldCoder.get(cellID, "layer");
    //Check if we want to use this later, else go to the next hit
    if (!useLayer(caloLayout, layer))
      continue;
    //Do the digitalization
    float calibr_coeff = 1.;
    calibr_coeff       = m_calibrCoeffCal;
    float hitEnergy    = calibr_coeff * energy;
    if (hitEnergy > m_maxHitEnergyCal) {
      hitEnergy = m_maxHitEnergyCal;
    }
    if (hitEnergy > m_thresholdCal) {
      edm4hep::MutableCalorimeterHit calHit = calcol.create();
      calHit.setCellID(cellID);
      calHit.setEnergy(hitEnergy);
      calHit.setPosition(hit.getPosition());
      calHit.setType(CHT(CHT::muon, CHT::yoke, caloLayout, layer));
      auto muonRel = muonRelcol.create();
      muonRel.setFrom(calHit);
      muonRel.setTo(hit);
    }
  }

  return std::make_tuple(std::move(calcol), std::move(calRelcol));
}

//StatusCode DualCrysCalDigi::finalize() { return StatusCode::SUCCESS; }

//If the vectors are empty, we are keeping everything
bool DualCrysCalDigi::useLayer(CHT::Layout caloLayout, unsigned int layer) const {
  switch (caloLayout) {
    case CHT::ecal:
      if (layer > m_useLayersEcalVec.size() || m_useLayersEcalVec.size() == 0)
        return true;
      return m_useLayersEcalVec[layer];  //break not needed, because of return
    case CHT::hcal:
      if (layer > m_useLayersHcalVec.size() || m_useLayersHcalVec.size() == 0)
        return true;
      return m_useLayersHcalVec[layer];  //break not needed, because of return
      //For all other cases, always keep the hit
    default:
      return true;
  }
}  //useLayer

