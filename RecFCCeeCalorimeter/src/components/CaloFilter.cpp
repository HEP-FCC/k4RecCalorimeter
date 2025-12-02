/** @class CaloFilterFunc
 * Gaudi MultiTransformer to apply a filter to a digitized pulse, defined as a TimeSeriesCollection object.
 *
 * @author Sahibjeet Singh
 * @date   2025-12-01
 *
 * Gaudi MultiTransformer that applies a filter to a TimeSeriesCollection of digitized pulses per cell, in order to extract the output of this filter, the most likely energy, and the timing information for each cell as an index of the pulse.
 *
 * This is a MultiTransformer that in v01 applies a matched filter to the digitized pulses. The matched filter is defined using the time-reversed and conjugate of a known pulse shape (currently a Gaussian) of length \f$N\f$. This is then convolved with the digitized pulse in order to extract the best estimate of the energy and timing information for each cell. The matched filter provides the best signal to noise ratio in the presence of additive stochastic noise. More information on the matched filter can be found at: https://en.wikipedia.org/wiki/Matched_filter and for this specific application, see: https://indico.cern.ch/event/1580025/contributions/6686602/attachments/3133485/5559196/FCCDigitization4BNLWorkshopEndOfWeekUpdate.pdf.
 *
 * In mathematical terms, a correctly normalized matched filter is defined as: 
 * \f$ \vec{M} = \frac{\Sigma^{-1} \vec{s}}{s^T \Sigma^{-1} s} \f$
 * \f$ \vec{M} \f$ is the vector defining the matched filter.
 * \f$ \Sigma^{-1} \f$ is the inverse of the noise correlation matrix. Must be provided via a ROOT file during initialization.
 * \f$ \vec{s} \f$ is the time reversed conjugate of the known signal shape. Currently assumed to be a Gausian shape with a size < the number of samples in the digitized pulse.
 *
 * The output of the matched filter is then computed by convolving \f$ \vec{M} \f$ with the digitized pulse \f$ \vec{D} \f$ as: \f$ \vec{O} = \vec{M} * \vec{D} \f$
 *
 * Note that the length of \f$ \vec{O} \f$ is set to be of the same size as \f$ \vec{D} \f$ in the convolution.
 * The maximum value of \f$ \vec{O} \f$ is taken as the best estimate of the energy in the cell, and the index of this maximum value is taken as the best estimate of the timing.
 *
 * Future iterations of this MultiTransformer will include other filter types.
 *
 *The unit system used here is GeV and ns throughout.
 * Inputs:
 *     - TimeSeriesCollection: Collection of digitized pulses per cell (see CaloDigitizerFunc for more detail).
 *
 * Properties:
 *     - @param m_filterName: The name of the filter to apply. Currently only "Matched_Gaussian" is implemented.
 *     - @param m_mu: Mean of the Gaussian pulse shape in ns.
 *     - @param m_sigma: Standard deviation of the Gaussian pulse shape in ns.
 *     - @param m_filterTemplateSize: The size of the known pulse shape to use for the matched filter. This should be a number less than the number of samples in the digitized pulse. Odd number are preferred since even numbers can lead to migration effects. This is typically (~4 or 5 in LAr).
 *     - @param m_pulseInitTime: Start time of the digitized pulse in ns.
 *     - @param m_pulseEndTime: End time of the digitized pulse in ns
 *     - @param m_lenSample: Number of samples in the full digitized pulse (typically ~30 in LAr).
 *     - @param m_noiseInfoFileName: The name of the ROOT file to load the noise correlation matrix from.
 *     - @param m_invCorrMatName: The name of the ROOT TMatrixD that represents the inverse of the correlation matrix of the noise.
 *
 * Outputs:
 *     - TimeSeriesCollection: Collection of digitized pulses after applying the filter.
 *     - int: Index representing the best estimate of the timing information.
 *     - float: The best estimate of the energy in the cell.
 *
 * LIMITATIONS: (status 01/12/2025)
 *     - Only the matched filter, assuming a Gaussian signal pulse shape, is currently implemented. The OFC method might be added in the future.
 *     - A nicer way to provide the signal shape would be good to have in the future, rather than hardcoding a Gaussian shape.
 *     - The matched filter should technically be a conjugate but since we are dealing with real numbers only, this is not done. Beware if complex pulse shapes are to be used in the future.
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

#include <cstddef>
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



    /**
    * \brief Computes the quantities necessary for the filter.
    *
    * This function computes the filter template (signal like part) based on the specified filter type (currently only matched filter with Gaussian pulse shape is implemented) and the provided parameters. The inverse correlation matrix is also loaded from a ROOT file for use in the filter computation.
    *
    * \param None -> See class parameters for inputs.
    * \return StatusCode indicating success or failure.
    */
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

            for (size_t i = 0; i < FilterTemplate->size(); ++i) {
                info() << "FilterTemplate[" << i << "] = " << (*FilterTemplate)[i] << endmsg;
            }

        } else {
            error() << "Unknown filter name: " << m_filterName.value() << endmsg;
            return StatusCode::FAILURE;
        }
        
        return StatusCode::SUCCESS;
    }
   
    /**
    * \brief Applies the matched filter to the digitized pulses.
    *
    * This function takes as input a collection of digitized pulses and applies the matched filter to each pulse. The filtered pulses, along with the best estimate of the energy (max(filtered sample)) and timing (sample index of the energy) for each pulse, are returned.
    *
    * \param DigitsPulse: The input collection of digitized pulses.
    * \return A tuple containing the filtered pulses, sample index representing the timing estimate, and the best estimate for the energy.
    */
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


    /**    
    * \brief Finalizes the transformer but does nothing in this case.
    * \return StatusCode indicating success or failure.
    */
    StatusCode finalize() override{return StatusCode::SUCCESS;}

    /**
    * \brief Read inverse correlation matrix from file.
    *
    * \param None -> uses class parameters for inputs.
    * \return StatusCode indicating success or failure.
    */
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

    /**
    * \brief Apply the matched filter to a digitized pulse.
    *
    * This function takes as input a digitized pulse and the matched filter template. It convolves the digitized pulse with the time reversed version of the matched filter template and returns the filtered pulse. 
    *
    * \param pulse: The input digitized pulse.
    * \param Cfilter: The matched filter template. This should NOT be time reversedyet.
    * \return Vector<float> of the filtered pulse after applying the matched filter.
    */
    std::vector<float> applyMatchedFilter(const  podio::RelationRange<float>& pulse, const std::vector<float>& Cfilter) const {
        // Reverse filter
        std::vector<float> filter = Cfilter;
        std::reverse(filter.begin(), filter.end());
        // Convolve the pulse with the filter
        return ConvolveSame(pulse, filter);
        }

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
    float GaussianNoNorm(float x, float mu = 0.0, float sigma = 1.0) const {
        return exp(-0.5 * pow((x - mu) / sigma, 2));
    }

    /**
    * \brief Computes the convolution of two vectors and returns a result of the same size as the input pulse.
    *
    * This function computes the convolution of the input pulse with the provided filter and returns a vector of the same size as the input pulse. This is the equivalent of np.convolve(..., mode='same') in Python.
    * \param pulse: The input digitized pulse.
    * \param filter: The filter to convolve with the pulse.
    * \return Vector<float> representing the convolved pulse of the same size as the input pulse.
    */
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