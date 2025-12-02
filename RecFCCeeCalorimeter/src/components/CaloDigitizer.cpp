/** @class CaloDigitizerFunc
 * Gaudi Transformer for ALLEGRO digitization
 *
 * @author Sahibjeet Singh
 * @date   2025-12-01
 *
 * Gaudi Transformer that digitises a SimCalorimeterHitCollection from noble liquid calorimeter to edm4hep::TimeSeriesCollection of digitized hits.
 *
 * The digitization uses a Gaussian pulse shape to simulate the response of the readout electronics. The mean and sigma of the Gaussian pulse can be defined via properties.
 * The digitization in v01 follows the LAr digitization as defined in https://cds.cern.ch/record/1057879/files/larg-pub-2007-011.pdf?version=1. Note that currently, the pedestal is not conisdered, noise is a work in progress, the ADC2MEV conversion is not considered, and the sampling fraction correction is a work in progress.
 * A more detailed explaination of the digitization procedure in v01, with the noise simulation can be found at: https://indico.cern.ch/event/1580025/contributions/6686602/attachments/3133485/5559196/FCCDigitization4BNLWorkshopEndOfWeekUpdate.pdf.
 * Note that fundamentally, each hit belonging to the geometric area of a given cell is converted to a digitized pulse in the corresponding cell, which are then summed to give the final digitized pulse in that cell.
 *
 *The unit system used here is GeV and ns throughout.
 * Inputs:
 *     - SimCalorimeterHitCollection: Collection of simulated calorimeter hits from ddsim/G4.
 *
 * Properties:
 *     - @param m_pulseType: The name of pulse shape to use for digitization. Currently only "Gaussian" is implemented.
 *     - @param m_mu: Mean of the Gaussian pulse shape in ns.
 *     - @param m_sigma: Standard deviation of the Gaussian pulse shape in ns.
 *     - @param m_lenSample: Number of samples in the digitized pulse.
 *     - @param m_pulseInitTime: Start time of the digitized pulse in ns.
 *     - @param m_pulseEndTime: End time of the digitized pulse in ns.
 *
 * Outputs:
 *     - TimeSeriesCollection: Digitised hits collection for each cell
 *
 * LIMITATIONS: (status 01/12/2025)
 *     - The digitization does not yet include the full noise simulation but a simpler one, in conjunction with AddNoise2Digits.cpp.
 *     - The digitization does not yet include the sampling fraction correction, but this will be added in the future via a separate algorithm.
 *     - The digitization does not yet include the ADC2MEV conversion
 *     - A more realistic pulse shape will be added in the future.
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
        
    /**
    * \brief Computes the pulse shape and its derivative.
    *
    * This function computes the pulse shape and its derivative based on the specified pulse type (currently only Gaussian is implemented) and the provided mean and sigma values. The pulse shape and its derivative are stored in member variables for later use during digitization.
    *
    * \param None -> See class parameters for inputs.
    * \return StatusCode indicating success or failure.
    */
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
   
   /**
    * \brief Digitizes the simulated calorimeter hits from Geant4.
    *
    * This function takes as input a SimCalorimeterHitCollection and digitizes each hit using the predefined pulse shape and its derivative. The digitized pulses are stored in a TimeSeriesCollection, which is returned as output.
    *
    * \param CaloHits: The input SimCalorimeterHitCollection.
    * \return TimeSeriesCollection of digitized hits.
    */
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

   /**
    * \brief Function to convert a hit to a digitized pulse.
    *
    * This function takes as input the energy and time of a hit, along with the pulse shape and its derivative, and computes the corresponding digitized pulse samples. The digitized pulse is represented as a vector of floats. See ATLAS LAr digitization note for details.
    *
    * \param Energy: The energy of the hit.
    * \param time: The time of the hit.
    * \param pulseShape: The pulse shape vector.
    * \param pulseShapeDeriv: The derivative of the pulse shape vector.
    * \return Vector<float> representing the digitized pulse samples.
    */
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

    /**
    * \brief Finalizes the transformer but does nothing in this case.
    * \return StatusCode indicating success or failure.
    */
   StatusCode finalize() override{return StatusCode::SUCCESS;}

    /**
    * \brief Computes the value of a Gaussian function analytically.
    *
    * This function computes \f$ G(x) = \frac{1}{\sigma \sqrt{2 \pi}} e^{-\frac{1}{2} \left( \frac{x - \mu}{\sigma} \right)^2} \f$ for given x, mu, and sigma.
    *
    * \param x: The input variable.
    * \param mu: The mean of the Gaussian.
    * \param sigma: The standard deviation of the Gaussian.
    * \return The value of the Gaussian function at x.
    */
   float Gaussian(float x, float mu = 0.0, float sigma = 1.0) const {
        return (1.0 / (sigma * sqrt(2 * M_PI))) * exp(-0.5 * pow((x - mu) / sigma, 2));
    }


    /**
    * \brief Computes the derivative of a Gaussian function analytically.
    * 
    * This function computes the derivative of the Gaussian \f$ G'(x) = -\left(\frac{x - \mu}{\sigma^2}\right) \cdot G(x) \f$
    *
    * \param x: The input variable.
    * \param mu: The mean of the Gaussian.
    * \param sigma: The standard deviation of the Gaussian.
    * \return The value of the derivative of the Gaussian function at x.
    */
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