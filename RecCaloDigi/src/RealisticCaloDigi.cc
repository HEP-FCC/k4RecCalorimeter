// Calorimeter digitiser
#include "RealisticCaloDigi.h"

#include <edm4hep/MutableCalorimeterHit.h>
#include <edm4hep/CaloHitSimCaloHitLinkCollection.h>
#include <edm4hep/CaloHitContribution.h>
#include <CalorimeterHitType.h>

#include <DDSegmentation/BitFieldCoder.h>

#include <k4Interface/IUniqueIDGenSvc.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <assert.h>
#include <cmath>
#include <set>

#include "CLHEP/Units/PhysicalConstants.h"


using namespace std;
using namespace std::placeholders;

struct MCC {
  float energy {0.f};
  float time {0.f};
};

RealisticCaloDigi::RealisticCaloDigi(const std::string& name, ISvcLocator* svcLoc) : MultiTransformer(name, svcLoc, 
    { KeyValues("inputHitCollections", {"SimCalorimeterHits"}),
      KeyValues("inputHeaderCollections", {"EventHeader"}) },
    { KeyValues("outputHitCollections", {"CalorimeterHits"}),
     KeyValues("outputRelationCollections", {"CaloHitLinks"}) }) {}

StatusCode RealisticCaloDigi::initialize() {
  m_geoSvc = serviceLocator()->service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to retrieve the GeoSvc" << endmsg;
    return StatusCode::FAILURE;
  }
  m_uidSvc = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
  if (!m_uidSvc) {
    error() << "Unable to get UniqueIDGenSvc" << endmsg;
  }


  // unit in which threshold is specified
  if (m_threshold_unit.value().compare("MIP") == 0){
    m_threshold_iunit=MIP;
  } else if (m_threshold_unit.value().compare("GeV") == 0){
    m_threshold_iunit=GEVDEP;
  } else if (m_threshold_unit.value().compare("px") == 0){
    m_threshold_iunit=NPE;
  } else {
    error() << "could not identify threshold unit. Please use \"GeV\", \"MIP\" or \"px\"! Aborting." << endmsg;
  }
  
  // convert the threshold to the approriate units (i.e. MIP for silicon, NPE for scint)
  m_threshold_value = convertEnergy( m_threshold_value, m_threshold_iunit );

  // deal with timing calculations  
  std::map<std::string, integr_function> integrations = {
    {"Standard", std::bind(&RealisticCaloDigi::StandardIntegration, this, _1)},
    {"ROC", std::bind(&RealisticCaloDigi::ROCIntegration, this, _1)}
  };  
  auto findIter = integrations.find( m_integration_method ) ;
  if(integrations.end() == findIter) {
    error() << "Could not guess timing calculation method!" << endmsg;
    error() << "Available are: Standard, ROC. Provided: " << m_integration_method << endmsg;
    error() << "Aborting..." << endmsg;
  }
  m_integr_function = findIter->second;
  
  // check if parameters are correctly set for the ROC integration
  if("ROC" == m_integration_method) {
    if(m_fast_shaper == 0.0f || m_slow_shaper == 0.0f) {
      error() << "Fast/slow shaper parameter(s) not set. Required for ROC integration!" << endmsg;
      error() << "Aborting..." << endmsg;
    }
  }

  return StatusCode::SUCCESS;
}



std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection> RealisticCaloDigi::operator()(
      const edm4hep::SimCalorimeterHitCollection& inputSim,
      const edm4hep::EventHeaderCollection& headers) const {
  auto seed = m_uidSvc->getUniqueID(headers[0].getEventNumber(), headers[0].getRunNumber(), this->name());
  debug() << "Using seed " << seed << " for event " << headers[0].getEventNumber() << " and run "
          << headers[0].getRunNumber() << endmsg;
  m_engine.SetSeed(seed);

  // decide on this event's correlated miscalibration
  float event_correl_miscalib = ( m_misCalib_correl>0 ) ? m_engine.Gaus(1.0, m_misCalib_correl) : 0;

  edm4hep::CalorimeterHitCollection newcol;
  edm4hep::CaloHitSimCaloHitLinkCollection relcol;

  CHT::CaloType cht_type = caloTypeFromString(m_calo_type);
  CHT::CaloID   cht_id   = caloIDFromString(m_calo_id);
  CHT::Layout   cht_lay  = layoutFromString(m_calo_layout);

  std::string initString;
  initString = m_geoSvc->constantAsString(m_encodingStringVariable.value());
  dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(initString);  // check!

  debug() << "Number of elements = " << inputSim.size() << endmsg;
  // loop over input hits
  for (int j=0; j < inputSim.size(); ++j) {
    edm4hep::SimCalorimeterHit simhit0 = inputSim.at( j );
    edm4hep::SimCalorimeterHit *simhit = &simhit0;

    // deal with energy integration and timing aspects
    auto integrationResult = Integrate(simhit);
    if( ! integrationResult.has_value() ) {
      continue;
    }
    float time      = integrationResult.value().first;
    float energyDep = integrationResult.value().second;
    // apply extra energy digitisation onto the energy
    float energyDig = EnergyDigi(energyDep, event_correl_miscalib);

    if (energyDig > m_threshold_value) { // write out this hit
      edm4hep::MutableCalorimeterHit newhit = newcol.create();
      newhit.setCellID( simhit->getCellID() );
      newhit.setTime( time );
      newhit.setPosition( simhit->getPosition() );
	    newhit.setEnergy( energyDig );
      
	    int layer = bitFieldCoder.get(simhit->getCellID(), "layer");
	    newhit.setType( CHT( cht_type, cht_id, cht_lay, layer ) );

      debug() << "orig/new hit energy: " << simhit->getEnergy() << " " << newhit.getEnergy() << endmsg;

      edm4hep::MutableCaloHitSimCaloHitLink rel = relcol.create();
      rel.setTo(simhit0);
      rel.setFrom(newhit);
      rel.setWeight(1.0);

    } // theshold
  } // input hits 

  return std::make_tuple(std::move(newcol), std::move(relcol));
}

//------------------------------------------------------------------------------

StatusCode RealisticCaloDigi::finalize(){
  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------

RealisticCaloDigi::integr_res_opt RealisticCaloDigi::Integrate( const edm4hep::SimCalorimeterHit * hit ) const {
  return m_integr_function(hit);
}

//------------------------------------------------------------------------------

float RealisticCaloDigi::EnergyDigi(float energy, float event_correl_miscalib) const{
  // some extra digi effects
  // controlled by _applyDigi = 0 (none), 1 (apply)
  // input parameters: hit energy ( in any unit: effects are all relative )
  // returns energy ( in units determined by the overloaded digitiseDetectorEnergy )

  float e_out(energy);
  e_out = digitiseDetectorEnergy(energy); // this is an overloaded method, provides energy in technology-dependent units

  // the following make only relative changes to the energy

  // random miscalib, uncorrelated in cells
  if (m_misCalib_uncorrel>0) {
    float miscal(0);
    miscal = m_engine.Gaus(1.0, m_misCalib_uncorrel);
    e_out*=miscal;
  }

  // random miscalib, correlated across cells in one event
  if (m_misCalib_correl>0) e_out*=event_correl_miscalib;

  float oneMipInMyUnits = convertEnergy( 1.0, MIP );
  // limited electronics dynamic range
  if ( m_elec_rangeMip > 0 ) e_out = std::min ( e_out, m_elec_rangeMip*oneMipInMyUnits );
  // add electronics noise
  if ( m_elec_noiseMip > 0 ) {
    e_out += m_engine.Gaus(0, m_elec_noiseMip*oneMipInMyUnits);
  }

  // random cell kill
  if (m_deadCell_fraction>0) {
    if (m_engine.Uniform(0., 1.) < m_deadCell_fraction ) e_out=0;
  }
  return e_out;
}

//------------------------------------------------------------------------------

RealisticCaloDigi::integr_res_opt RealisticCaloDigi::StandardIntegration( const edm4hep::SimCalorimeterHit * hit ) const {
  // apply timing cuts on simhit contributions
  // outputs a (time,energy) pair
  float timeCorrection(0);
  if ( m_time_correctForPropagation ) { // time of flight from IP to this point
    float r = pow(hit->getPosition().x,2) + pow(hit->getPosition().y,2) + pow(hit->getPosition().z,2); 
    timeCorrection = sqrt(r)/CLHEP::c_light; // [speed of light in mm/ns]
  }
  // this is Oskar's simple (and probably the most correct) method for treatment of timing
  //  - collect energy in some predefined time window around collision time (possibly corrected for TOF)
  //  - assign time of earliest contribution to hit
  float energySum = 0;
  float earliestTime=std::numeric_limits<float>::max();
  for(edm4hep::CaloHitContribution contribution : hit->getContributions()){ // loop over all contributions
    float timei   = contribution.getTime(); //absolute hit timing of current subhit
    float energyi = contribution.getEnergy(); //energy of current subhit
    float relativetime = timei - timeCorrection; // wrt time of flight
    if (relativetime>m_time_windowMin && relativetime<m_time_windowMax){
      energySum += energyi;
      if (relativetime<earliestTime){
	       earliestTime = relativetime; //use earliest hit time for simpletimingcut
      }
    }
  }
  if(earliestTime > m_time_windowMin && earliestTime < m_time_windowMax){ //accept this hit
    return integr_res{SmearTime(earliestTime), energySum};
  }
  return std::nullopt;
}

//------------------------------------------------------------------------------

RealisticCaloDigi::integr_res_opt RealisticCaloDigi::ROCIntegration( const edm4hep::SimCalorimeterHit * hit ) const {
  const unsigned int ncontrib = hit->contributions_size() ;
  // Sort MC contribution by time
  std::vector<MCC> mcconts{ncontrib};
  for(int i=0; i<hit->getContributions().size();i++){
    mcconts[i].energy = hit->getContributions(i).getEnergy();
    mcconts[i].time = hit->getContributions(i).getTime();
  }
  std::sort(mcconts.begin(), mcconts.end(), [](auto lhs, auto rhs){
    return (lhs.time < rhs.time);
  });
  // Accumulate energy until threshold is reached.
  // The first MC contriubtion after the threshold has been reached sets the hit time 
  bool passThreshold = false;
  float epar=0.f, hitTime=0.f;
  unsigned int thresholdIndex=0;
  // First determine the hit time (hitTime) and the initial hit index 
  // at which we need to start the integration (thresholdIndex)
  for(unsigned int i=0; i<ncontrib ; ++i) {
    const auto timei = mcconts[i].time;
    thresholdIndex = i ;
    epar = 0.f;
    for(unsigned int j=i; j<ncontrib ; ++j) {
      const auto timej = mcconts[j].time;
      if( (timej-timei) < m_fast_shaper) {
        epar += mcconts[j].energy;
      }
      else {
        break;
      }
      if( convertEnergy(epar, GEVDEP) > m_threshold_value ) {
        hitTime = timej;
        passThreshold = true ;
        break;
      }
    }
    if(passThreshold) {
      break;
    }
  }
  // check hit time 
  const float thresholdTime = mcconts[thresholdIndex].time;
  if( not (thresholdTime>m_time_windowMin && thresholdTime<m_time_windowMax) ) {
    return std::nullopt;
  }
  // If we've found a hit above the threshold, accumulate the energy until
  // until the maximum time given by the slow shaper
  if(passThreshold) {
    float energySum = 0.f ; 
    for(unsigned int i=thresholdIndex ; i<ncontrib ; ++i) {
      if( mcconts[i].time < thresholdTime + m_slow_shaper) {
        energySum += mcconts[i].energy;
      }
    }
    hitTime = SmearTime(hitTime);
    return integr_res{hitTime, energySum};
  }
  // else no hit dude !
  else {
    return std::nullopt;
  }
}

//------------------------------------------------------------------------------

float RealisticCaloDigi::SmearTime(float time) const{
  return m_time_resol>0.f ? time + m_engine.Gaus(0, m_time_resol) : time;
}

