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

#include <ostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
 
 struct CaloFilterFunc final
     : k4FWCore::Transformer<edm4hep::TimeSeriesCollection(const edm4hep::TimeSeriesCollection&)> {
        CaloFilterFunc(const std::string& name, ISvcLocator* svcLoc)
       : Transformer(name, svcLoc, 
        {KeyValues("InputCollection", {"TimeSeriesCollection"})},
        {KeyValues("OutputCollection", {"FilteredDigits"})}) {}
 
    StatusCode initialize() override{
           return StatusCode::SUCCESS;
    }
   
    // This is the function that will be called to transform the data
   // Note that the function has to be const, as well as all pointers to collections
   // we get from the input
   edm4hep::TimeSeriesCollection operator()(const edm4hep::TimeSeriesCollection& DigitsPulse) const override {
     info() << "Digitized pulse collection size: " << DigitsPulse.size() << endmsg;

    edm4hep::TimeSeriesCollection FilteredDigitsCollection;
    
    // Loop over DigitsPulse to extract the pulse amplitudes
    for (const auto& Digit : DigitsPulse) {
        std::vector<double> amplitude;
        const auto InputPulse = Digit.getAmplitude();

        auto FilteredDigit = FilteredDigitsCollection.create();
        FilteredDigit.setCellID(Digit.getCellID());

        FilteredDigit.setTime(0.0); // Placeholder for time info
        FilteredDigit.setInterval(Digit.getInterval()); // Set the interval for the digitized pulse in ns

        // Apply the matched filter to the pulse
        std::vector<float> FilterTemplate;
        if (m_filterName.value() == "MatchedDirac") {
            debug() << "Using the matched filter: Matched Dirac" << endmsg;
            FilterTemplate = {0.0, 1.0, 0.0};
        } else {
            error() << "Unknown filter name: " << m_filterName.value() << endmsg;
        }

        auto Out = applyMatchedFilter(InputPulse, FilterTemplate);

        for (unsigned int i = 0; i < Out.size(); i++) {
            FilteredDigit.addToAmplitude(Out[i]);
        }

    }

    return FilteredDigitsCollection;
}

   StatusCode finalize() override{return StatusCode::SUCCESS;}


 
 private:  
  /// Map to be used for the lookup of the pulse shapes
  Gaudi::Property<std::string> m_filterName{this, "filterName", "Matching", "Name of the filter to apply" };

  // Method to apply matched filter
  std::vector<float> applyMatchedFilter(const  podio::RelationRange<float>& pulse, std::vector<float>& filter) const {
      std::vector<float> filteredPulse(pulse.size(), 0.0f);
      // Reverse filter
        std::reverse(filter.begin(), filter.end());
      for (size_t i = 0; i < pulse.size(); ++i) {
          for (size_t j = 0; j < filter.size(); ++j) {
              if (i + j < pulse.size()) {
                  filteredPulse[i] += pulse[i - j] * filter[j];
              }
          }
      }
      return filteredPulse;
  }

 };
 
 DECLARE_COMPONENT(CaloFilterFunc)