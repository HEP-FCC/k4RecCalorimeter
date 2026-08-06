#ifndef DIGITIZER_DDCCALODIGISCINT_H
#define DIGITIZER_DDCCALODIGISCINT_H 1

#include "RealisticCaloDigi.h"


/** === RealisticCaloDigiScinPpd Processor === <br>
    realistic digitisation of scint+PPD (SiPM, MPPC) calorimeter hits
    D.Jeans 02/2016.
*/

struct RealisticCaloDigiScinPpd : public RealisticCaloDigi {
  
 public:
  RealisticCaloDigiScinPpd(const std::string& name, ISvcLocator* svcLoc);

 protected:
  int getMyUnit() const {return NPE;}
  float digitiseDetectorEnergy(float energy) const ; // apply scin+PPD specific effects
  float convertEnergy( float energy, int inputUnit ) const; // convert energy from input to output scale

  Gaudi::Property<float> m_PPD_pe_per_mip{this, "ppd_mipPe", 10.0f, "# Photo-electrons per MIP (scintillator): used to Poisson smear #PEs if >0"};
  Gaudi::Property<int> m_PPD_n_pixels{this, "ppd_npix", 10000, "Total number of MPPC/SiPM pixels for implementation of saturation effect"};
  Gaudi::Property<float> m_misCalibNpix{this, "ppd_npix_uncert", 0.05f, "Fractional uncertainty of effective total number of MPPC/SiPM pixels"};
  Gaudi::Property<float> m_pixSpread{this, "ppd_pix_spread", 0.05f, "Variation of PPD pixel signal (as a fraction: 0.01=1%)"};
} ;

#endif



