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
 
 struct CaloDigitizerFunc final
     : k4FWCore::Transformer<edm4hep::TimeSeriesCollection(const edm4hep::SimCalorimeterHitCollection&)> {
        CaloDigitizerFunc(const std::string& name, ISvcLocator* svcLoc)
       : Transformer(name, svcLoc, 
        {KeyValues("InputCollection", {"SimCaloHitsCollection"})},
        {KeyValues("OutputCollection", {"DigitsFloat"})}) {}
 
    StatusCode initialize() override{
        // Decide which pulse shape to use and create the pulse shape and derivative vectors on the fly
        m_samplingInterval = (m_pulseEndTime.value() - m_pulseInitTime.value()) / m_lenSample.value();
        
        if (std::strcmp(m_pulseType.value().c_str(), "Gaussian") == 0) {
            debug() << "Using Gaussian pulse shape!" << endmsg;

            for (size_t i = 0; i < m_lenSample.value(); ++i) {
                PulseShape.push_back(Gaussian(i * m_samplingInterval, m_mu.value(), m_sigma.value()));
                PulseShapeDeriv.push_back(GaussianDerivative(i * m_samplingInterval, m_mu.value(), m_sigma.value()));
                info() << "PulseShape[" << i * m_samplingInterval << "]: " << PulseShape[i] << ", PulseShapeDeriv[" << i * m_samplingInterval << "]: " << PulseShapeDeriv[i] << endmsg;
            }
        } else {
            error() << "Unknown pulse type: " << m_pulseType.value() << endmsg;
        }
        return StatusCode::SUCCESS;
    }
   
   // This is the function that will be called to transform the data
   // Note that the function has to be const, as well as all pointers to collections
   // we get from the input
   edm4hep::TimeSeriesCollection operator()(const edm4hep::SimCalorimeterHitCollection& CaloHits) const override {
     info() << "calorimeter hits collection size: " << CaloHits.size() << endmsg;

     edm4hep::TimeSeriesCollection DigitsCollection;

    // Loop over CaloHits to accumulate energy
    for (const auto& hit : CaloHits) {
        uint64_t cellID = hit.getCellID();
        float energy = hit.getEnergy();
        debug() << "Hit energy: " << energy << endmsg;
        auto Contributions = hit.getContributions();   

        std::vector<float> DigitVectorSum(PulseShape.size(), 0.0f);

        // Loop over contributions to get the energy and time
        for (const auto& contribution : Contributions) {
            auto DigitVector = Hit2DigitFloat(contribution.getEnergy(), contribution.getTime(), PulseShape, PulseShapeDeriv);
            
            // Sum the contributions
            std::transform(DigitVectorSum.begin(), DigitVectorSum.end(), DigitVector.begin(), DigitVectorSum.begin(), std::plus<float>());

            // Print the contribution info
            debug() << "Contribution (energy, time): (" << contribution.getEnergy() << ", " << contribution.getTime() << ")" << endmsg;
        }
        auto Digit = DigitsCollection.create();
        Digit.setCellID(cellID);
        Digit.setTime(0.0); // Placeholder for time info
        Digit.setInterval(m_samplingInterval); // Set the interval for the digitized pulse in ns

        for (unsigned int i = 0; i < DigitVectorSum.size(); i++) {
            Digit.addToAmplitude(DigitVectorSum[i]);
        }

    }
     return DigitsCollection;
   }

   std::vector<float> Hit2DigitFloat(float Energy, float time, const std::vector<float>& pulseShape, const std::vector<float>& pulseShapeDeriv) const
   {     
        // Basing this off of LArDigitization pub note: https://inspirehep.net/files/8b92316ea8786d3954cccbfaa0d10ada
        std::vector<float> DigitVector(pulseShape.size(), 0.0f);

        // Loop over the DigitVector
        for (unsigned int i = 0; i < DigitVector.size(); ++i) {
            int j = i - std::rint(time / m_samplingInterval);
            // if j < 0, it means the pulse is not in the time window so set it to 0
            if (j >= 0) {
                float DeltaT = time - m_samplingInterval * std::rint(time / m_samplingInterval);
                                
                DigitVector[i] = Energy * (pulseShape[j] - DeltaT * pulseShapeDeriv[j]);
            }
        }
        return DigitVector;
   }

   StatusCode finalize() override{return StatusCode::SUCCESS;}

   float Gaussian(float x, float mu = 0.0, float sigma = 1.0) const {
        return (1.0 / (sigma * sqrt(2 * M_PI))) * exp(-0.5 * pow((x - mu) / sigma, 2));
    }

    float GaussianDerivative(float x, float mu = 0.0, float sigma = 1.0) const {
        return - (x - mu) / (sigma * sigma) * Gaussian(x, mu, sigma);
    }


 
 private:
  // Type of pulse to create
  Gaudi::Property<std::string> m_pulseType {this, "pulseType", "Gaussian", "Type of pulse to create" };
  // Initial time of the pulse
  Gaudi::Property<float> m_pulseInitTime {this, "pulseInitTime", 0.0, "Initial time of the pulse" };
  // End time of the pulse
  Gaudi::Property<float> m_pulseEndTime {this, "pulseEndTime", 750.0, "End time of the pulse" };
  // Number of samples in pulse
  Gaudi::Property<int> m_lenSample {this, "pulseSamplingLength", 30, "Number of samples in pulse" };
  // Gaussian pulse properties
  Gaudi::Property<float> m_mu {this, "mu", 100.0, "Mean of Gaussian pulse" };
  Gaudi::Property<float> m_sigma {this, "sigma", 20.0, "Sigma of Gaussian pulse" };
  
  std::vector<float> PulseShape;
  std::vector<float> PulseShapeDeriv;

  float m_samplingInterval;

 };
 
 DECLARE_COMPONENT(CaloDigitizerFunc)