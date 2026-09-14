#include "CellEnergyCalibratorSiliconMip.h"

DECLARE_COMPONENT(CellEnergyCalibratorSiliconMip)

CellEnergyCalibratorSiliconMip::CellEnergyCalibratorSiliconMip(const std::string& name, ISvcLocator* svcLoc)
    : BaseCellEnergyCalibrator(name, svcLoc) {}

float CellEnergyCalibratorSiliconMip::reconstructEnergy(const edm4hep::CalorimeterHit& hit, int layer) const {
  // here the input energy should be in MIPs
  float energy = hit.getEnergy();
  // now correct for sampling fraction
  energy *= getLayerCalib(layer);
  return energy;
}
