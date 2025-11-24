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
#include <TMatrixD.h>
#include <TVectorD.h>

// For correlation calculation
#include <TMatrixDSym.h>
#include <cmath>
#include <vector>

 struct CaloAddNoise2Digits final
     : k4FWCore::Transformer<edm4hep::TimeSeriesCollection(const edm4hep::TimeSeriesCollection&)> {
        CaloAddNoise2Digits(const std::string& name, ISvcLocator* svcLoc)
       : Transformer(name, svcLoc, 
        {KeyValues("InputCollection", {"DigitsFloat"})},
        {KeyValues("OutputCollection", {"DigitsWithNoiseFloat"})}) {}
 
    StatusCode initialize() override{
        r3->SetSeed(m_noiseSeed.value()); // Set the seed for the random number generator

        TMatrixD NoiseSampleMat(m_noiseSimSamples.value(), m_lenSample.value());
        std::vector<double> means(m_lenSample.value(), 0.0); // Pulse size to calculate means
        createNoiseSampleMatrix(&NoiseSampleMat, &means);
        
        TMatrixDSym covMat = ComputeCovarianceMatrix(NoiseSampleMat, means);
        TMatrixDSym corMat = ComputeCorrelationMatrix(covMat);
        TMatrixDSym invCorrMat(corMat);
        invCorrMat.Invert(); // Invert correlation matrix to then save

        TVectorD meansVec(m_lenSample.value());
        for (int i = 0; i < m_lenSample.value(); ++i) {
            meansVec[i] = means[i];
            info() << "MeansVec[" << i << "]: " << means[i] << endmsg;
        }

        // Check if file exists
        if (m_noiseInfoFileName.empty()) {
            error() << "Name of the file with the pulse shapes not provided!" << endmsg;
            return StatusCode::FAILURE;
        }

        info() << "Saving noise sample matrix, means vector, covariance matrix, correlation matrix, and inverse correlation matrix to: " << m_noiseInfoFileName.value() << endmsg;
        auto f = TFile::Open(m_noiseInfoFileName.value().c_str(), "RECREATE");
        f->cd();
        NoiseSampleMat.Write("NoiseSampleMatrix");

        meansVec.Write("NoiseSampleMat_MeansVector");
        covMat.Write("CovarianceMatrix");
        corMat.Write("CorrelationMatrix");
        invCorrMat.Write("InvertedCorrelationMatrix");
        f->Close();

        info() << "Noise sample matrix created with size: " << NoiseSampleMat.GetNrows() << "x" << NoiseSampleMat.GetNcols() << endmsg;
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
            OutVector[i] = Digits[i] + r3->Gaus(NoiseEnergy, NoiseWidth);
        }
        return OutVector;
   }


    // Computes the correlation coefficient matrix from a TMatrixD of observations
    // data: TMatrixD with rows = observations, cols = variables
    TMatrixDSym ComputeCovarianceMatrix(const TMatrixD& data, const std::vector<double>& means) {
        // Compute covariance matrix
        int nVars = data.GetNcols();
        TMatrixDSym cov(nVars);
        for (int i = 0; i < nVars; ++i) {
            for (int j = i; j < nVars; ++j) {
                double sum = 0.0;
                for (int k = 0; k < m_noiseSimSamples.value(); ++k) {
                    sum += (data(k, i) - means[i]) * (data(k, j) - means[j]);
                }
                cov(i, j) = sum / (m_noiseSimSamples.value() - 1);
                if (i != j) cov(j, i) = cov(i, j);
            }
        }
        return cov;
    }

    TMatrixDSym ComputeCorrelationMatrix(const TMatrixDSym& cov){
        // Compute correlation matrix from cov mat
        int nVars = cov.GetNcols();
        TMatrixDSym corr(nVars);
        for (int i = 0; i < nVars; ++i) {
            for (int j = i; j < nVars; ++j) {
                double denom = std::sqrt(cov(i, i) * cov(j, j));
                debug() << "Covariance: " << cov(i, j) << ", Denominator: " << denom << std::endl;
                double val = denom != 0.0 ? cov(i, j) / denom : 0.0;
                corr(i, j) = val;
                if (i != j) corr(j, i) = val;
            }
        }
        return corr;
    }
    

    void createNoiseSampleMatrix(TMatrixD* NoiseSampleMat, std::vector<double>* means) const {
        for (int j = 0; j < m_lenSample.value(); ++j) {
            for (int i = 0; i < m_noiseSimSamples.value(); ++i) {
                (*NoiseSampleMat)(i, j) = r3->Gaus(m_noiseEnergy.value(), m_noiseWidth.value());
                debug() << "NoiseSampleMat(" << i << ", " << j << ") = " << (*NoiseSampleMat)(i, j) << std::endl;
                (*means)[j] += (*NoiseSampleMat)(i, j);
            }
            (*means)[j] /= m_noiseSimSamples.value(); // Calculate mean for each pulse sample
        }
    }

    StatusCode finalize() override{
        return StatusCode::SUCCESS;
    }


 private:
  Gaudi::Property<std::string> m_noiseInfoFileName{this, "noiseInfoFileName", "NoiseInfo.root", "Name of file to store noise samples, means, cov, and corr matrices"};
  Gaudi::Property<float> m_noiseEnergy {this, "noiseEnergy", 0.1, "Noise energy for Gaussian (mean) - in GeV" };
  Gaudi::Property<float> m_noiseWidth {this, "noiseWidth", 0.05, "Noise width - in GeV as well" };
  Gaudi::Property<int> m_noiseSimSamples {this, "noiseSimSamples", 1000, "Number of samples for noise simulation" };
  Gaudi::Property<int> m_noiseSeed {this, "noiseSeed", 32, "Seed for the random number generator" };
  Gaudi::Property<int> m_lenSample {this, "pulseSamplingLength", 30, "Number of samples in pulse" };

  TRandom *r3 = new TRandom3();

 };
 
 DECLARE_COMPONENT(CaloAddNoise2Digits)
