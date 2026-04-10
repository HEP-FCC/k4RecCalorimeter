/** @class CaloWhitening
 * Gaudi Transformer to apply a whitening filter to a TimeSeriesCollection
 *
 * @author Sahibjeet Singh
 * @date   2025-12-01
 *
 * Gaudi Transformer that applies a whitening filter to a TimeSeriesCollection of digitized pulses per cell.
 *
 * The whitening filter is used to decorrelate noise between samples in the digitized pulses. The whitening matrix is computed from the inverse of the noise correlation matrix, which is provided via a ROOT file during initialization. The mean vector of the noise is also provided via the same ROOT file and is subtracted from each digitized pulse before applying the whitening filter. The whitening operation is performed using the ZCA (Zero-phase Component Analysis) method but there are other methods that can be implemented in the future.
 *
 * In mathematical terms, the ZCA method is: 
 * \f$ \vec{W} = \Sigma^{-1/2} (\vec{D} - \vec{\mu})\f$
 * \f$ \vec{W} \f$ is a vector of the whitened pulse.
 * \f$ \Sigma^{-1/2} \f$ is the whitening matrix, defined as the inverse square root of the noise correlation matrix. This should be precomputed and provided via the ROOT file.
 * \f$ \vec{D} \f$ is a vector of the digitized pulse.
 * \f$ \vec{\mu} \f$ is the mean vector of the noise. This should also be provided via the ROOT file.
 * Note that the vectors \f$ \vec{W} \f$, \f$ \vec{D} \f$, and \f$ \vec{\mu} \f$ have a length equal to the number of samples in the digitized pulse and are thus defined per sample.
 * More details for v01 be found at: https://indico.cern.ch/event/1580025/contributions/6686602/attachments/3133485/5559196/FCCDigitization4BNLWorkshopEndOfWeekUpdate.pdf.
 *
 *The unit system used here is GeV and ns throughout.
 * Inputs:
 *     - TimeSeriesCollection: Collection of digitized pulses per cell (see CaloDigitizerFunc for more detail).
 *
 * Properties:
 *     - @param m_noiseInfoFileName: The name of the ROOT file to load the noise correlation matrix from.
 *     - @param m_invCorrMatName: The name of the ROOT TMatrixD that represents the inverse of the correlation matrix of the noise.
 *     - @param m_muVecName: The name of the ROOT TVectorD that represents the mean vector of the noise.
 *    - @param m_filterName: The name of the whitening filter to apply. Currently only "ZCA" is implemented. 
 *
 * Outputs:
 *     - TimeSeriesCollection: Collection of digitized pulses after applying the whitening filter.
 *
 * LIMITATIONS: (status 01/12/2025)
 *     - Only the ZCA whitening method is currently implemented.
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
 
    /**
    * \brief Computes the quantities necessary for the filter.
    *
    * This function computes the stuff needed for the whitening filter (currently only ZCA whitening is implemented) and the provided parameters. The inverse correlation matrix is also loaded from a ROOT file for use in the filter computation.
    *
    * \param None -> See class parameters for inputs.
    * \return StatusCode indicating success or failure.
    */
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
            WhiteningMatrix = new TMatrixD(eigenVecs, TMatrixD::kMult, sqrtDiag);
            (*WhiteningMatrix) *= eigenVecs.T();

            if (msgLevel(MSG::DEBUG)) {
                debug() << "Eigenvalues: " << endmsg;
                eigenVals.Print();

                debug() << "Eigenvecs: " << endmsg;
                eigenVecs.Print();

                debug() << "sqrt(corr^-1): " << endmsg;
                WhiteningMatrix->Print();
            }
            info() << "Sqrt of Corr^(-1) done!" << endmsg;
        }

        // Cholesky decomposition to be added in the future
        // else if (m_filterName.value().c_str() == "Cholesky"){

        // }

        else{
            error() << "Unknown filter name: " << m_filterName.value() << endmsg;
        }

        return StatusCode::SUCCESS;
    }
   


   /**
    * \brief Applies the whitening filter to the digitized pulses.
    *
    * This function takes as input a TimeSeriesCollection of digitized pulses per cell and applies the whitening filter to each pulse. The resulting TimeSeriesCollection with whitened pulses is returned.
    *
    * \param DigitsPulse: The input TimeSeriesCollection of digitized pulses.
    * \return TimeSeriesCollection of whitened pulses.
    */
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
            debug() << "WhitenedPulse[" << i << "] = " << WhitenedPulse[i] << endmsg;
        }
    }
     return WhitenedDigitsCollection;
   }

   /**
    * \brief Applies the whitening filter to the digitized pulses.
    *
    * This function applies the whitening filter to a single digitized pulse by subtracting the mean vector and multiplying by the whitening matrix.
    * \f$ \vec{W} = \Sigma^{-1/2} (\vec{D} - \vec{\mu})\f$
    *
    * \param Digits: The input digitized pulse.
    * \return TVectorD representing the whitened pulse.
    */
   TVectorD ApplyWhiteningMat(const podio::RelationRange<float> Digits) const
   {
       TVectorD Pulse(Digits.size());
       // Loop over the DigitVector
       for (unsigned int i = 0; i < Digits.size(); ++i) {
           Pulse[i] = Digits[i] - (*MuVec)[i]; // Subtract mean
            debug() << "Digits w/ noise, no mu subtraction[" << i << "] = " << Digits[i] << endmsg;
           debug() << "Pulse after Mu subtraction[" << i << "] = " << Pulse[i] << endmsg;
       }
       auto WhitenedPulse = (*WhiteningMatrix) * Pulse;
       return WhitenedPulse;
   }


    /**    
    * \brief Finalizes the transformer but does nothing in this case.
    * \return StatusCode indicating success or failure.
    */
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
