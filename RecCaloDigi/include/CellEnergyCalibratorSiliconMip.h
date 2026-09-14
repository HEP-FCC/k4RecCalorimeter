#ifndef CELLENERGYCALIBRATORSILICONMIP_H
#define CELLENERGYCALIBRATORSILICONMIP_H 1

#include "BaseCellEnergyCalibrator.h"

/** === CellEnergyCalibratorSiliconMip algorithm === <br>
    calibration of digitised silicon calorimeter hits
    input: MIP scale, output: shower GeV
    D.Jeans 02/2016.

    24 March 2016: removed gap corrections - to be put into separate processor

*/

struct CellEnergyCalibratorSiliconMip final : BaseCellEnergyCalibrator {

public:
  CellEnergyCalibratorSiliconMip(const std::string& name, ISvcLocator* svcLoc);

protected:
  float reconstructEnergy(const edm4hep::CalorimeterHit& hit, int layer) const override;
};

#endif
