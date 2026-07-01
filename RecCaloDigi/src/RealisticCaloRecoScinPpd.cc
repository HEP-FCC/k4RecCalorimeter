#include "RealisticCaloRecoScinPpd.h"
#include <cmath>
#include <iostream>

DECLARE_COMPONENT(RealisticCaloRecoScinPpd)

RealisticCaloRecoScinPpd::RealisticCaloRecoScinPpd(const std::string& name, ISvcLocator* svcLoc) : RealisticCaloReco(name, svcLoc) {}

float RealisticCaloRecoScinPpd::reconstructEnergy(const edm4hep::CalorimeterHit* hit, int layer) const{
  // here the input energy should be in NPE
  float energy = hit->getEnergy();

  // first de-saturate PPD response
  // this is the fraction of SiPM pixels fired above which a linear continuation of the saturation-reconstruction function is used. 
  // 0.95 of nPixel corresponds to a energy correction of factor ~3.
  const float r = 0.95;
  if (energy < r*m_PPD_n_pixels){ //current hit below linearisation threshold, reconstruct energy normally:
    energy = -m_PPD_n_pixels * std::log ( 1. - ( energy / m_PPD_n_pixels ) );
  } else { //current hit is aove linearisation threshold, reconstruct using linear continuation function:
    energy = 1/(1-r)*(energy-r*m_PPD_n_pixels)-m_PPD_n_pixels*std::log(1-r);
  }
  // then go back to MIP scale
  energy/=m_PPD_pe_per_mip;

  // now correct for sampling fraction (calibration from MIP -> shower GeV)
  energy *= getLayerCalib( layer );

  return energy;
}


