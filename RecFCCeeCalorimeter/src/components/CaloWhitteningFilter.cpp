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

#include "TFile.h"
#include "TSystem.h"

#include <TMatrixD.h>
#include <TVectorD.h>

// For correlation calculation
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <cmath>

 struct CaloWhitening final
     : k4FWCore::Transformer<edm4hep::TimeSeriesCollection(const edm4hep::TimeSeriesCollection&)> {
        CaloWhitening(const std::string& name, ISvcLocator* svcLoc)
       : Transformer(name, svcLoc, 
        {KeyValues("InputCollection", {"DigitsFloat"})},
        {KeyValues("OutputCollection", {"WhitenedDigitsCollection"})}) {}
 
    StatusCode initialize() override{
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
        TMatrixDSym* invCorrMat = nullptr;
        NoiseInfoFile->GetObject(m_invCorrMatName.value().c_str(), invCorrMat);
        if (!invCorrMat) {
            error() << "Unable to load inverse of correlation matrix from file!" << endmsg;
            return StatusCode::FAILURE;
        }

        // Load muVec into memory
        NoiseInfoFile->GetObject(m_muVecName.value().c_str(), MuVec);
        if (!MuVec) {
            error() << "Unable to load mean vector from file!" << endmsg;
            return StatusCode::FAILURE;
        }

        if (m_filterName.value() == "ZCA"){
            info() << "Using the ZCA filter!" << endmsg;
            // Define ZCA whitening matrix --> corr^(-1/2)
            // Compute eigen decomposition
            TMatrixDSymEigen eig(*invCorrMat);

            info() << "Doing eigendecomposition of inverse correlation matrix to get the sqrt" << endmsg;
            // Get eigenvalues and eigenvectors
            TVectorD eigenVals = eig.GetEigenValues();
            TMatrixD eigenVecs = eig.GetEigenVectors();

            // Build diagonal matrix of sqrt(eigenvalues)
            TMatrixD sqrtDiag(eigenVals.GetNrows(), eigenVals.GetNrows());
            sqrtDiag.Zero();
            for (int i = 0; i < eigenVals.GetNrows(); ++i) {
                sqrtDiag(i, i) = std::sqrt(eigenVals[i]);
            }

            // Compute square root: V * sqrt(D) * V^T
            // TMatrixD WhiteningMat = TMatrixD(eigenVecs, TMatrixD::kMult, sqrtDiag); // Does not overwrite anything. ROOT is stupid......
            // WhiteningMat *= eigenVecs.T(); // Overwrites WhiteningMatrix w/ WhiteningMat *= eigenVecs.T()
            // WhiteningMatrix = &WhiteningMat; // Assign to the member variable

            WhiteningMatrix = new TMatrixD(eigenVecs, TMatrixD::kMult, sqrtDiag);
            (*WhiteningMatrix) *= eigenVecs.T();

            // info() << "Eigenvalues: " << endmsg;
            // eigenVals.Print();

            // info() << "Eigenvecs: " << endmsg;
            // eigenVecs.Print();

            // info() << "sqrt(corr^-1): " << endmsg;
            // WhiteningMatrix->Print();

            info() << "Sqrt of Corr^(-1) done!" << endmsg;
        }

        // else if (m_filterName.value().c_str() == "Cholesky"){

        // }

        else{
            error() << "Unknown filter name: " << m_filterName.value() << endmsg;
        }

        return StatusCode::SUCCESS;
    }
   


   // This is the function that will be called to transform the data
   // Note that the function has to be const, as well as all pointers to collections
   // we get from the input
   edm4hep::TimeSeriesCollection operator()(const edm4hep::TimeSeriesCollection& DigitsPulse) const override {
    info() << "Input digitized pulse collection size: " << DigitsPulse.size() << endmsg;

    edm4hep::TimeSeriesCollection WhitenedDigitsCollection;

    // Loop over DigitsPulse to extract the pulse amplitudes
    for (const auto& Digit : DigitsPulse) {
        const auto InputPulse = Digit.getAmplitude();

        auto WhitenedDigit = WhitenedDigitsCollection.create();
        WhitenedDigit.setCellID(Digit.getCellID());

        WhitenedDigit.setTime(0.0); // Placeholder for time info
        WhitenedDigit.setInterval(Digit.getInterval()); // Set the interval for the digitized pulse in ns

        // Apply whitening to the digitized pulse
        auto WhitenedPulse = ApplyWhiteningMat(InputPulse);



        for (unsigned int i = 0; i < WhitenedPulse.size(); i++) {
            WhitenedDigit.addToAmplitude(WhitenedPulse[i]);
            // std::cout << "WhitenedPulse[" << i << "] = " << WhitenedPulse[i] << std::endl;
        }
    }
     return WhitenedDigitsCollection;
   }

   TVectorD ApplyWhiteningMat(const podio::RelationRange<float> Digits) const
   {
       TVectorD Pulse(Digits.size());
       // Loop over the DigitVector
       for (unsigned int i = 0; i < Digits.size(); ++i) {
           Pulse[i] = Digits[i] - (*MuVec)[i]; // Subtract mean
        //    std::cout << "Digits w/ noise, no mu subtraction[" << i << "] = " << Digits[i] << std::endl;
        //    std::cout << "Pulse after Mu subtraction[" << i << "] = " << Pulse[i] << std::endl;
       }
       auto WhitenedPulse = (*WhiteningMatrix) * Pulse;
       return WhitenedPulse;
   }


    StatusCode finalize() override{
        return StatusCode::SUCCESS;
    }


 private:
  Gaudi::Property<std::string> m_noiseInfoFileName{this, "noiseInfoFileName", "NoiseInfo.root", "Name of file to load noise corr matrix from"};
  Gaudi::Property<std::string> m_invCorrMatName{this, "invCorrMatName", "InvertedCorrelationMatrix", "Name of ROOT TMatrixD that represents the inverse of the correlation matrix of the noise"};
  Gaudi::Property<std::string> m_muVecName{this, "muVecName", "NoiseSampleMat_MeansVector", "Name of ROOT TVectorD that represents the mean vector of the noise"};
  Gaudi::Property<std::string> m_filterName{this, "filterName", "ZCA", "Name of the whitening filter to apply"};
  TMatrixD* WhiteningMatrix = nullptr;
  TVectorD* MuVec = nullptr;

};
 
 DECLARE_COMPONENT(CaloWhitening)