// Calorimeter digitiser for the IDC ECAL and HCAL
// For other detectors/models SimpleCaloDigi should be used
#include "DigitiserSiliconMip.h"

#include <string>

using namespace std;

DECLARE_COMPONENT(DigitiserSiliconMip)

DigitiserSiliconMip::DigitiserSiliconMip(const std::string& name, ISvcLocator* svcLoc) : BaseDigitiser(name, svcLoc) {}

// a MIP-scale output can be reached from deposited GeV or from MIPs, but there is no way to
// get there from a photo-electron count: that needs the scintillator+PPD response model
bool DigitiserSiliconMip::canConvertFrom(EnergyScale inUnit) const {
  return inUnit == EnergyScale::MIP || inUnit == EnergyScale::GEVDEP;
}

// convert energy from input to output scale (MIP)
float DigitiserSiliconMip::convertEnergy(float energy, EnergyScale inUnit) const {
  // only the scales accepted by canConvertFrom() reach this point
  if (inUnit == EnergyScale::GEVDEP)
    return energy / m_calib_mip;
  return energy; // already on the MIP scale
}

float DigitiserSiliconMip::digitiseDetectorEnergy(float energy) const {
  // applies extra digitisation to silicon hits
  //  input energy in deposited GeV
  //  output is MIP scale
  float smeared_energy(energy);
  if (m_ehEnergy > 0) {
    // calculate #e-h pairs
    float nehpairs = 1e9 * energy / m_ehEnergy; // check units of energy! _ehEnergy is in eV, energy in GeV
    // fluctuate it by Poisson (actually an overestimate: Fano factor actually makes it smaller, however even this
    // overstimated effect is tiny for our purposes)
    smeared_energy *= m_engine.Poisson(nehpairs) / nehpairs;
  }

  return smeared_energy / m_calib_mip; // convert to MIP units
}
