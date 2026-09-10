#ifndef DIGITISERSILICONMIP_H
#define DIGITISERSILICONMIP_H 1

#include "BaseDigitiser.h"

/** === DigitiserSiliconMip algorithm === <br>
    digitisation of silicon calorimeter hits
    input: deposited GeV, output: MIP scale
    D.Jeans 02/2016.
*/

struct DigitiserSiliconMip : public BaseDigitiser {

public:
  DigitiserSiliconMip(const std::string& name, ISvcLocator* svcLoc);

protected:
  // convert energy from the input scale to the output (MIP) scale
  float convertEnergy(float energy, EnergyScale inputUnit) const override;
  bool canConvertFrom(EnergyScale inputUnit) const override;
  // apply silicon-specific digitisation
  float digitiseDetectorEnergy(float energy) const override;

  Gaudi::Property<float> m_ehEnergy{this, "silicon_pairEnergy", 3.6f,
                                    "energy required to create e-h pair in silicon (in eV)"};
};

#endif
