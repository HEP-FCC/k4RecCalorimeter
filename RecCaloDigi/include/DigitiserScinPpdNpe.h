#ifndef DIGITISERSCINPPDNPE_H
#define DIGITISERSCINPPDNPE_H 1

#include "BaseDigitiser.h"

/** === DigitiserScinPpdNpe algorithm === <br>
    digitisation of scint+PPD (SiPM, MPPC) calorimeter hits
    input: deposited GeV, output: saturated photo-electrons (NPE)
    D.Jeans 02/2016.
*/

struct DigitiserScinPpdNpe : public BaseDigitiser {

public:
  DigitiserScinPpdNpe(const std::string& name, ISvcLocator* svcLoc);

protected:
  // apply scin+PPD specific effects
  float digitiseDetectorEnergy(float energy) const override;
  // convert energy from the input scale to the output (NPE) scale
  float convertEnergy(float energy, EnergyScale inputUnit) const override;
  bool canConvertFrom(EnergyScale inputUnit) const override;

  Gaudi::Property<float> m_PPD_pe_per_mip{this, "ppd_mipPe", 10.0f,
                                          "# Photo-electrons per MIP (scintillator): used to Poisson smear #PEs if >0"};
  Gaudi::Property<int> m_PPD_n_pixels{this, "ppd_npix", 10000,
                                      "Total number of MPPC/SiPM pixels for implementation of saturation effect"};
  Gaudi::Property<float> m_misCalibNpix{this, "ppd_npix_uncert", 0.05f,
                                        "Fractional uncertainty of effective total number of MPPC/SiPM pixels"};
  Gaudi::Property<float> m_pixSpread{this, "ppd_pix_spread", 0.05f,
                                     "Variation of PPD pixel signal (as a fraction: 0.01=1%)"};
};

#endif
