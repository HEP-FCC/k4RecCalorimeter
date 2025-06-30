/*
* Copyright (c) 2014-2024 Key4hep-Project.
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

#include "Gaudi/Property.h"

// edm4hep
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/Constants.h"
#include "edm4hep/RawTimeSeriesCollection.h"
#include "edm4hep/TimeSeriesCollection.h"
#include "edm4hep/CaloHitContributionCollection.h"

#include "k4FWCore/Transformer.h"

#include <algorithm>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"

#include <TRandom3.h>
 
 struct CaloAddNoise2Digits final
     : k4FWCore::Transformer<edm4hep::TimeSeriesCollection(const edm4hep::TimeSeriesCollection&)> {
        CaloAddNoise2Digits(const std::string& name, ISvcLocator* svcLoc)
       : Transformer(name, svcLoc, 
        {KeyValues("InputCollection", {"DigitsFloat"})},
        {KeyValues("OutputCollection", {"DigitsWithNoiseFloat"})}) {}
 
    StatusCode initialize() override{
        r3->SetSeed(0); // Set the seed for the random number generator
        return StatusCode::SUCCESS;
    }
   
   // This is the function that will be called to transform the data
   // Note that the function has to be const, as well as all pointers to collections
   // we get from the input
   edm4hep::TimeSeriesCollection operator()(const edm4hep::TimeSeriesCollection& DigitsPulse) const override {
    info() << "Digitized pulse collection size: " << DigitsPulse.size() << endmsg;

    edm4hep::TimeSeriesCollection DigitsWNoiseCollection;

    // Loop over DigitsPulse to extract the pulse amplitudes
    for (const auto& Digit : DigitsPulse) {
        const auto InputPulse = Digit.getAmplitude();

        auto DigitWNoise = DigitsWNoiseCollection.create();
        DigitWNoise.setCellID(Digit.getCellID());

        DigitWNoise.setTime(0.0); // Placeholder for time info
        DigitWNoise.setInterval(Digit.getInterval()); // Set the interval for the digitized pulse in ns

        // Apply noise to the digitized pulse
        auto Out = applyGaussianNoise(InputPulse, m_noiseEnergy, m_noiseWidth);

        for (unsigned int i = 0; i < Out.size(); i++) {
            DigitWNoise.addToAmplitude(Out[i]);
        }
    }
     return DigitsWNoiseCollection;
   }

   std::vector<float> applyGaussianNoise(const podio::RelationRange<float> Digits, float NoiseEnergy, float NoiseWidth) const
   {     
        std::vector<float> OutVector(Digits.size(), 0.0f);

        // Loop over the DigitVector
        for (unsigned int i = 0; i < OutVector.size(); ++i) {
            OutVector[i] = Digits[i] + (NoiseEnergy * r3->Gaus(0, NoiseWidth));
        }
        return OutVector;
   }

   StatusCode finalize() override{return StatusCode::SUCCESS;}

 
 private:
  Gaudi::Property<float> m_noiseEnergy {this, "noiseEnergy", 0.1, "Noise energy to scale Gaussian by - GeV" };
  Gaudi::Property<float> m_noiseWidth {this, "noiseWidth", 1, "Noise width" };
  TRandom *r3 = new TRandom3();

 };
 
 DECLARE_COMPONENT(CaloAddNoise2Digits)