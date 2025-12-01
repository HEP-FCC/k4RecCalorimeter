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
#ifndef DualCrysCalDigi_H
#define DualCrysCalDigi_H

#include "Gaudi/Property.h"
#include "edm4hep/CaloHitContributionCollection.h"
#include "edm4hep/CaloHitSimCaloHitLinkCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

#include "CalorimeterHitType.h"
#include "DualCrysCalorimeterHit.h"
#include "DDRec/SurfaceManager.h"
#include "k4FWCore/Transformer.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

#include <random>
#include <string>
#include <vector>


// inside the <> for the multitransformer, the structure is
// <output type ( input type ), baseclass template >
// so CalorimeterHitCollection and CaloHitSimCaloHitLink Collection are the output types
// SimCalorimeterHit Collection& is the base class for these?
// and somehow the baseclass template is not used
struct DualCrysCalDigi final
    : k4FWCore::MultiTransformer<
          std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>(
              const edm4hep::SimCalorimeterHitCollection&, const edm4hep::EventHeaderCollection&)> {
  
  DualCrysCalDigi(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  // StatusCode finalize() override;

  std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection> operator()(
      const edm4hep::SimCalorimeterHitCollection& simCaloHits,
      const edm4hep::EventHeaderCollection&       headers) const override;

private:

private:
  // declare useLayer here:
  bool useLayer(CHT::Layout caloLayout, unsigned int layer) const;

  Gaudi::Property<std::string> m_calCollections{this, "calCollections", "DRCNoSegment",
                                                "The input collection of calorimeters"};
  Gaudi::Property<std::string> outputRelCollection{this, "outputRelCollection", "outputRelCollection",
                                                  "The output collection of relations"};
  Gaudi::Property<std::string> outputCalCollection{this, "outputCalCollection", "outputCalCollection",
                                                   "The output collection of calorimeters"};

  Gaudi::Property<std::string> m_encodingStringVariable{
      this, "EncodingStringParameterName", "GlobalTrackerReadoutID",
      "The name of the DD4hep constant that contains the encoding string for tracking detectors"};

  Gaudi::Property<float> m_thresholdCal{this, "CalThreshold", 0.025, "Threshold for calorimeters"};
  Gaudi::Property<float> m_calibrCoeffCal{this, "calibrationCoeffCal", 120000.0, "Calibration coefficient of calorimeters"};
  Gaudi::Property<float> m_maxHitEnergyCal{this, "maxCalHitEnergy", 2.0, "Threshold for maximum calorimeter hit energy"};

  Gaudi::Property<std::string> m_detectorNameEcal{this, "detectorNameEcal", "DRCrystal", "Name of ECAL"};
  Gaudi::Property<std::string> m_detectorNameHcal{this, "detectorNameHcal", "DRFtubeFiber", "Name of HCAL"};
  Gaudi::Property<std::vector<bool>> m_useLayersEcalVec{this, "useLayersEcal", {}, "Enable/disable ECAL layers"};
  Gaudi::Property<std::vector<bool>> m_useLayersHcalVec{this, "useLayersHcal", {}, "Enable/disable HCAL layers"};

  std::string m_collName;

  SmartIF<IGeoSvc> m_geoSvc;
  SmartIF<IUniqueIDGenSvc> m_uidSvc;
};

DECLARE_COMPONENT(DualCrysCalDigi)
#endif

