/*
* Copyright (c) 2014-2025 Key4hep-Project.
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

#include "podio/UserDataCollection.h"

#include "k4FWCore/Transformer.h"

#include <iostream>
#include <ostream>
#include <string>

#include "TFile.h"
#include "TSystem.h"
#include <TMatrixD.h>
#include <TVectorD.h>
#include "TMatrixDSym.h"

using DigitsColl = edm4hep::TimeSeriesCollection;

using FilterColl = edm4hep::TimeSeriesCollection;
using MatchedSampleIdxColl = podio::UserDataCollection<int>;
using MatchedSampleEnergyColl = podio::UserDataCollection<float>;

struct CaloFilterFunc final
: k4FWCore::MultiTransformer<std::tuple<FilterColl, MatchedSampleIdxColl, MatchedSampleEnergyColl>(const DigitsColl&)>{
        CaloFilterFunc(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                        {KeyValues("InputCollection", {"TimeSeriesCollection"})},

                        {
                            KeyValues("OutputCollectionFilteredPulse", {"FilteredDigitsCollection"}),
                            KeyValues("OutputCollectionMatchedSampleIdx", {"MatchedSampleIdx"}),
                            KeyValues("OutputCollectionMatchedSampleEnergy", {"MatchedSampleEnergy"})
                        }) {}


    StatusCode initialize() override{
        // Get matched filter template
        if (m_filterName.value() == "Matched_Gaussian"){
            info() << "Using the matched filter: Matched Gaussian" << endmsg;          

            // Number of samples you want the known pulse (i.e. filter) to consider
            // Note that an odd number of samples > 1 is the best for a Gaussian. Otherwise, migration effects will happen

            // If even number of samples then up count > down count
            int DownCnt = (m_filterTemplateSize.value() - 1)/2;
            int UpCnt = m_filterTemplateSize.value() - DownCnt - 1;

            // We will create a Gaussian shape w/ mean = m_mu and sigma = m_sigma. The filter will be defined as e^(-0.5 * ((x - m_mu)/m_sigma)^2) for x in [maximum - DownCnt, maximum + UpCnt].

            int samplingInterval = (m_pulseEndTime.value() - m_pulseInitTime.value()) / m_lenSample.value();
            std::vector<float> PulseShape(m_lenSample.value(), 0.0);
            for (int i = 0; i < m_lenSample.value(); ++i) {
                PulseShape[i] = Gaussian(i * samplingInterval, m_mu.value(), m_sigma.value());
                info() << "PulseShape[" << i << "] = " << PulseShape[i] << endmsg;
            }

            // Get the filter template between MaxIdx - DownCnt and MaxIdx + UpCnt
            auto MaxIdx = std::distance(PulseShape.begin(), std::max_element(PulseShape.begin(), PulseShape.end()));

            int DownIdx = MaxIdx - DownCnt;
            int UpIdx = MaxIdx + UpCnt;
            // Validate indices
            if (DownIdx < 0 || UpIdx >= static_cast<int>(PulseShape.size())) {
                error() << "Index out of bounds: DownIdx=" << DownIdx << ", UpIdx=" << UpIdx
                        << ", PulseShape.size()=" << PulseShape.size() << endmsg;
                return StatusCode::FAILURE;
            }
            
            info() << "MaxIdx: " << MaxIdx << ", DownIdx: " << DownIdx << ", UpIdx: " << UpIdx << endmsg;

            // Create padded FilterTemplate with small part that is signal
            TVectorD FilterTemplatePadded(m_lenSample.value());
            for (int i = DownIdx; i <= UpIdx; ++i) {
                FilterTemplatePadded[i] = PulseShape[i];
            }

            for (int i = 0; i < FilterTemplatePadded.GetNrows(); ++i) {
                info() << "FilterTemplatePadded[" << i << "] = " << FilterTemplatePadded[i] << endmsg;
            }

            auto GetInvCorrMatStatus = GetInvCorrMat();
            if (GetInvCorrMatStatus != StatusCode::SUCCESS) {
                return GetInvCorrMatStatus;
            }

            auto Temp = (*invCorrMat) * FilterTemplatePadded;
            float Norm = FilterTemplatePadded * Temp;

            info() << "Norm for matched filter (PaddedSignal^T * InvCorrMat * PaddedSignal): " << Norm << endmsg;

            for (int i = 0; i < Temp.GetNrows(); ++i) {
                info() << "Temp[" << i << "] = " << Temp[i] << endmsg;
            }
            FilterTemplate = new std::vector<float>(m_filterTemplateSize.value(), 0.0f);
            for (int i = DownIdx; i <= UpIdx; ++i) {
                (*FilterTemplate)[i - DownIdx] = Temp[i] / Norm;
            } // Should be subvector of Temp --> Temp[DownIdx:UpIdx], divided by Norm of course

            for (int i = 0; i < FilterTemplate->size(); ++i) {
                info() << "FilterTemplate[" << i << "] = " << (*FilterTemplate)[i] << endmsg;
            }

        } else {
            error() << "Unknown filter name: " << m_filterName.value() << endmsg;
            return StatusCode::FAILURE;
        }
        
        return StatusCode::SUCCESS;
    }
   
    // This is the function that will be called to transform the data
    // Note that the function has to be const, as well as all pointers to collections
    // we get from the input
    std::tuple<FilterColl, MatchedSampleIdxColl, MatchedSampleEnergyColl> operator()(const DigitsColl& DigitsPulse) const override {
        info() << "Digitized pulse collection size: " << DigitsPulse.size() << endmsg;

        DigitsColl FilteredDigitsCollection;
        MatchedSampleIdxColl MaxIdxCollection;
        MatchedSampleEnergyColl EnergyCollection;
    
        // Loop over DigitsPulse to extract the pulse amplitudes
        for (const auto& Digit : DigitsPulse) {
            const auto InputPulse = Digit.getAmplitude();

            auto FilteredDigit = FilteredDigitsCollection.create();
            FilteredDigit.setCellID(Digit.getCellID());

            FilteredDigit.setTime(0.0); // Placeholder for time info
            FilteredDigit.setInterval(Digit.getInterval()); // Set the interval for the digitized pulse in ns

            // Apply matched filter
            auto Out = applyMatchedFilter(InputPulse, *FilterTemplate);

            for (unsigned int i = 0; i < Out.size(); i++) {
                FilteredDigit.addToAmplitude(Out[i]);
            }

            // Calculate the energy of the matched filter and the matched sample index
            auto MaxIdx = std::distance(Out.begin(), std::max_element(Out.begin(), Out.end()));
            auto Energy = *std::max_element(Out.begin(), Out.end());

            debug() << "Cell ID" << Digit.getCellID() << ", MaxIdx: " << MaxIdx << ", Energy - max val: " << Energy << endmsg;

            // Store the matched sample index and energy
            MaxIdxCollection.push_back(MaxIdx);
            EnergyCollection.push_back(Energy);


        }

    return std::make_tuple(std::move(FilteredDigitsCollection), std::move(MaxIdxCollection), std::move(EnergyCollection));
    }

   StatusCode finalize() override{return StatusCode::SUCCESS;}

   StatusCode GetInvCorrMat(){
        // Check if file exists
        if (m_noiseInfoFileName.empty()) {
            error() << "Name of the file with the noise info not provided!" << endmsg;
            return StatusCode::FAILURE;
        }
        if (gSystem->AccessPathName(m_noiseInfoFileName.value().c_str())) {
            error() << "Provided file with the noise info not found!" << endmsg;
            error() << "File path: " << m_noiseInfoFileName.value() << endmsg;
            return StatusCode::FAILURE;
        }

        // Open the file with the noise info
        std::unique_ptr<TFile> NoiseInfoFile(TFile::Open(m_noiseInfoFileName.value().c_str(), "READ"));
        if (NoiseInfoFile->IsZombie()) {
            error() << "Unable to open the file with the noise info!" << endmsg;
            error() << "File path: " << m_noiseInfoFileName.value() << endmsg;
            return StatusCode::FAILURE;
        } else {
            info() << "Using the following file with noise info: "
                << m_noiseInfoFileName.value() << endmsg;
        }

        // Load CorrMat into memory
        NoiseInfoFile->GetObject(m_invCorrMatName.value().c_str(), invCorrMat);
        // invCorrMat->Print();
        if (!invCorrMat) {
            error() << "Unable to load inverse of correlation matrix from file!" << endmsg;
            return StatusCode::FAILURE;
        }

        return StatusCode::SUCCESS;
    }

    // Method to apply matched filter
    std::vector<float> applyMatchedFilter(const  podio::RelationRange<float>& pulse, const std::vector<float>& Cfilter) const {
        // Reverse filter
        std::vector<float> filter = Cfilter;
        std::reverse(filter.begin(), filter.end());
        // Convolve the pulse with the filter
        return ConvolveSame(pulse, filter);
        }

    float Gaussian(float x, float mu = 0.0, float sigma = 1.0) const {
        return (1.0 / (sigma * sqrt(2 * M_PI))) * exp(-0.5 * pow((x - mu) / sigma, 2));
    }

    float GaussianNoNorm(float x, float mu = 0.0, float sigma = 1.0) const {
        return exp(-0.5 * pow((x - mu) / sigma, 2));
    }

    std::vector<float> ConvolveSame(const podio::RelationRange<float>& pulse, const std::vector<float>& filter) const {
        int n = pulse.size();
        int m = filter.size();
        int conv_size = n + m - 1;
        int start = (m - 1) / 2;
        int end = start + n;

        std::vector<float> result(conv_size, 0.0);
        for (size_t i = 0; i < pulse.size(); ++i) {
            for (size_t j = 0; j < filter.size(); ++j) {
                result[i + j] += pulse[i] * filter[j];
            }
        }
        std::vector<float> result_same(result.begin() + start, result.begin() + end);
        return result_same;
    }


 private:  
  /// Map to be used for the lookup of the pulse shapes
  Gaudi::Property<std::string> m_filterName{this, "filterName", "Matched_Gaussian", "Name of the filter to apply" };
    
  Gaudi::Property<float> m_mu{this, "mu", 50, "Mean of Gaussian pulse" };
  Gaudi::Property<float> m_sigma{this, "sigma", 20, "Sigma of Gaussian pulse" };
  
  // Number of samples in the filter template
  Gaudi::Property<int> m_filterTemplateSize{this, "filterTemplateSize", 5, "Size of the filter template" };

  // Initial time of the pulse
  Gaudi::Property<float> m_pulseInitTime {this, "pulseInitTime", 0.0, "Initial time of the pulse" };
  // End time of the pulse
  Gaudi::Property<float> m_pulseEndTime {this, "pulseEndTime", 750.0, "End time of the pulse" };
  // Number of samples in pulse
  Gaudi::Property<int> m_lenSample {this, "pulseSamplingLength", 30, "Number of samples in pulse" };
  Gaudi::Property<std::string> m_noiseInfoFileName{this, "noiseInfoFileName", "NoiseInfo.root", "Name of file to load noise corr matrix from"};
  Gaudi::Property<std::string> m_invCorrMatName{this, "invCorrMatName", "InvertedCorrelationMatrix", "Name of ROOT TMatrixD that represents the inverse of the correlation matrix of the noise"};

  TMatrixDSym* invCorrMat = nullptr;
  std::vector<float>* FilterTemplate = nullptr;

 };
 
 DECLARE_COMPONENT(CaloFilterFunc)