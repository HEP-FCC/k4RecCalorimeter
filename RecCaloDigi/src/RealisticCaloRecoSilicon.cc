#include "RealisticCaloRecoSilicon.h"
#include <algorithm>
#include <cassert>
#include <iostream>


DECLARE_COMPONENT(RealisticCaloRecoSilicon)

RealisticCaloRecoSilicon::RealisticCaloRecoSilicon(const std::string& name, ISvcLocator* svcLoc) : RealisticCaloReco(name, svcLoc) {}

float RealisticCaloRecoSilicon::reconstructEnergy(const edm4hep::CalorimeterHit* hit, int layer) const{
  // here the input energy should be in MIPs
  float energy = hit->getEnergy();
  // now correct for sampling fraction
  energy *= getLayerCalib( layer );
  return energy;
}
