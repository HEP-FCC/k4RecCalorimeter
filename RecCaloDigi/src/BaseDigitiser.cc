// Calorimeter digitiser
#include "BaseDigitiser.h"

#include "CalorimeterHitType.h"
#include <edm4hep/CaloHitContribution.h>
#include <edm4hep/CaloHitSimCaloHitLinkCollection.h>
#include <edm4hep/MutableCalorimeterHit.h>

#include <k4Interface/IUniqueIDGenSvc.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <string>

#include "CLHEP/Units/PhysicalConstants.h"

using namespace std;
using namespace std::placeholders;

struct MCC {
  float energy{0.f};
  float time{0.f};
};

BaseDigitiser::BaseDigitiser(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(
          name, svcLoc,
          {KeyValue("inputHitCollection", "SimCalorimeterHits"), KeyValue("inputHeaderCollection", "EventHeader")},
          {KeyValue("outputHitCollection", "CalorimeterHits"), KeyValue("outputRelationCollection", "CaloHitLinks")}) {}

StatusCode BaseDigitiser::initialize() {
  m_geoSvc = serviceLocator()->service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to retrieve the GeoSvc" << endmsg;
    return StatusCode::FAILURE;
  }
  m_uidSvc = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
  if (!m_uidSvc) {
    error() << "Unable to get UniqueIDGenSvc" << endmsg;
    return StatusCode::FAILURE;
  }

  // the cellID encoding does not change during the job, so decode it once here
  m_bitFieldCoder = dd4hep::DDSegmentation::BitFieldCoder(m_geoSvc->constantAsString(m_encodingStringVariable.value()));

  // both algorithms decode the layer number out of the cellID, so the encoding must provide it.
  // Catching this here turns a per-event exception into a clear configuration error.
  try {
    m_bitFieldCoder.index("layer");
  } catch (const std::exception&) {
    error() << "The encoding string '" << m_encodingStringVariable.value() << "' ("
            << m_bitFieldCoder.fieldDescription() << ") has no \"layer\" field" << endmsg;
    return StatusCode::FAILURE;
  }

  // unit in which threshold is specified
  if (m_threshold_unit.value().compare("MIP") == 0) {
    m_threshold_iunit = EnergyScale::MIP;
  } else if (m_threshold_unit.value().compare("GeV") == 0) {
    m_threshold_iunit = EnergyScale::GEVDEP;
  } else if (m_threshold_unit.value().compare("px") == 0) {
    m_threshold_iunit = EnergyScale::NPE;
  } else {
    error() << "Could not identify threshold unit '" << m_threshold_unit.value()
            << "'. Please use \"GeV\", \"MIP\" or \"px\"" << endmsg;
    return StatusCode::FAILURE;
  }

  // convert the threshold to the approriate units (i.e. MIP for silicon, NPE for scint)
  if (!canConvertFrom(m_threshold_iunit)) {
    error() << "thresholdUnit '" << m_threshold_unit.value()
            << "' cannot be expressed on the output scale of this digitiser" << endmsg;
    return StatusCode::FAILURE;
  }
  m_threshold_value = convertEnergy(m_threshold_value, m_threshold_iunit);

  // deal with timing calculations
  std::map<std::string, IntegrationFunction> integrations = {
      {"Standard", std::bind(&BaseDigitiser::standardIntegration, this, _1)},
      {"ROC", std::bind(&BaseDigitiser::rocIntegration, this, _1)}};
  auto findIter = integrations.find(m_integration_method);
  if (integrations.end() == findIter) {
    error() << "Could not guess timing calculation method!" << endmsg;
    error() << "Available are: Standard, ROC. Provided: " << m_integration_method << endmsg;
    return StatusCode::FAILURE;
  }
  m_integr_function = findIter->second;

  // Decode the calorimeter type/ID/layout once: these are fixed for the whole job.
  // The *FromString helpers silently fall back to a default when nothing matches, so without
  // these checks a typo in one of the properties would go unnoticed until someone wondered
  // why every hit came out typed "unknown".
  m_chtType = caloTypeFromString(m_calo_type);
  if (m_chtType == CHT::CaloType::unknown) {
    warning() << "CaloType '" << m_calo_type.value() << "' matches none of em, had, muon: hits will be typed as unknown"
              << endmsg;
  }

  m_chtId = caloIDFromString(m_calo_id);
  if (m_chtId == CHT::CaloID::unknown) {
    warning() << "CaloID '" << m_calo_id.value()
              << "' matches none of ecal, hcal, yoke, lcal, lhcal, bcal: hits will be typed as unknown" << endmsg;
  }

  m_chtLayout = layoutFromString(m_calo_layout);
  if (m_chtLayout == CHT::Layout::any) {
    warning() << "CaloLayout '" << m_calo_layout.value()
              << "' matches none of barrel, endcap, plug, ring: hits will be typed as any" << endmsg;
  }

  // check if parameters are correctly set for the ROC integration
  if ("ROC" == m_integration_method) {
    if (m_fast_shaper == 0.0f || m_slow_shaper == 0.0f) {
      error() << "Fast/slow shaper parameter(s) not set. Required for ROC integration!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>
BaseDigitiser::operator()(const edm4hep::SimCalorimeterHitCollection& inputSim,
                          const edm4hep::EventHeaderCollection& headers) const {
  m_engine.SetSeed(m_uidSvc->getUniqueID(headers, this->name()));

  // decide on this event's correlated miscalibration
  float event_correl_miscalib = (m_misCalib_correl > 0) ? m_engine.Gaus(1.0, m_misCalib_correl) : 0;

  edm4hep::CalorimeterHitCollection newcol;
  edm4hep::CaloHitSimCaloHitLinkCollection relcol;

  debug() << "Number of elements = " << inputSim.size() << endmsg;

  std::size_t nOutsideWindow = 0;
  std::size_t nBelowThreshold = 0;

  // loop over input hits
  for (const auto& simhit : inputSim) {
    // deal with energy integration and timing aspects
    const auto integrationResult = integrate(simhit);
    if (!integrationResult.has_value()) {
      // No contribution fell inside the integration window. This is expected for
      // out-of-time hits, so it is not worth more than a debug message.
      ++nOutsideWindow;
      debug() << "Dropping hit with cellID " << simhit.getCellID() << ": no contribution inside the integration window"
              << endmsg;
      continue;
    }
    float time = integrationResult->time;
    float energyDep = integrationResult->energy;
    // apply extra energy digitisation onto the energy
    float energyDig = energyDigi(energyDep, event_correl_miscalib);

    if (energyDig > m_threshold_value) { // write out this hit
      edm4hep::MutableCalorimeterHit newhit = newcol.create();
      newhit.setCellID(simhit.getCellID());
      newhit.setTime(time);
      newhit.setPosition(simhit.getPosition());
      newhit.setEnergy(energyDig);

      int layer = m_bitFieldCoder.get(simhit.getCellID(), "layer");
      newhit.setType(CHT(m_chtType, m_chtId, m_chtLayout, layer));

      debug() << "orig/new hit energy: " << simhit.getEnergy() << " " << newhit.getEnergy() << endmsg;

      auto rel = relcol.create();
      rel.setTo(simhit);
      rel.setFrom(newhit);
      rel.setWeight(1.0);

    } // threshold
    else {
      ++nBelowThreshold;
      debug() << "Dropping hit with cellID " << simhit.getCellID() << ": digitised energy " << energyDig
              << " below threshold " << m_threshold_value.value() << endmsg;
    }
  } // input hits

  debug() << "Digitised " << newcol.size() << " of " << inputSim.size() << " hits (" << nOutsideWindow
          << " outside the integration window, " << nBelowThreshold << " below threshold)" << endmsg;

  return std::make_tuple(std::move(newcol), std::move(relcol));
}

//------------------------------------------------------------------------------

BaseDigitiser::IntegrationResult BaseDigitiser::integrate(const edm4hep::SimCalorimeterHit& hit) const {
  return m_integr_function(hit);
}

//------------------------------------------------------------------------------

float BaseDigitiser::energyDigi(float energy, float event_correl_miscalib) const {
  // some extra digi effects
  // controlled by _applyDigi = 0 (none), 1 (apply)
  // input parameters: hit energy ( in any unit: effects are all relative )
  // returns energy ( in units determined by the overloaded digitiseDetectorEnergy )

  float e_out(energy);
  e_out = digitiseDetectorEnergy(energy); // this is an overloaded method, provides energy in technology-dependent units

  // the following make only relative changes to the energy

  // random miscalib, uncorrelated in cells
  if (m_misCalib_uncorrel > 0) {
    float miscal(0);
    miscal = m_engine.Gaus(1.0, m_misCalib_uncorrel);
    e_out *= miscal;
  }

  // random miscalib, correlated across cells in one event
  if (m_misCalib_correl > 0)
    e_out *= event_correl_miscalib;

  float oneMipInMyUnits = convertEnergy(1.0, EnergyScale::MIP);
  // limited electronics dynamic range
  if (m_elec_rangeMip > 0)
    e_out = std::min(e_out, m_elec_rangeMip * oneMipInMyUnits);
  // add electronics noise
  if (m_elec_noiseMip > 0) {
    e_out += m_engine.Gaus(0, m_elec_noiseMip * oneMipInMyUnits);
  }

  // random cell kill
  if (m_deadCell_fraction > 0) {
    if (m_engine.Uniform(0., 1.) < m_deadCell_fraction)
      e_out = 0;
  }
  return e_out;
}

//------------------------------------------------------------------------------

BaseDigitiser::IntegrationResult BaseDigitiser::standardIntegration(const edm4hep::SimCalorimeterHit& hit) const {
  // apply timing cuts on simhit contributions
  // outputs the accepted time and integrated energy
  float timeCorrection(0);
  if (m_time_correctForPropagation) { // time of flight from IP to this point
    float r = pow(hit.getPosition().x, 2) + pow(hit.getPosition().y, 2) + pow(hit.getPosition().z, 2);
    timeCorrection = sqrt(r) / CLHEP::c_light; // [speed of light in mm/ns]
  }
  // this is Oskar's simple (and probably the most correct) method for treatment of timing
  //  - collect energy in some predefined time window around collision time (possibly corrected for TOF)
  //  - assign time of earliest contribution to hit
  float energySum = 0;
  float earliestTime = std::numeric_limits<float>::max();
  for (const edm4hep::CaloHitContribution& contribution : hit.getContributions()) { // loop over all contributions
    float timei = contribution.getTime();        // absolute hit timing of current subhit
    float energyi = contribution.getEnergy();    // energy of current subhit
    float relativetime = timei - timeCorrection; // wrt time of flight
    if (relativetime > m_time_windowMin && relativetime < m_time_windowMax) {
      energySum += energyi;
      if (relativetime < earliestTime) {
        earliestTime = relativetime; // use earliest hit time for simpletimingcut
      }
    }
  }
  if (earliestTime > m_time_windowMin && earliestTime < m_time_windowMax) { // accept this hit
    return TimeEnergy{smearTime(earliestTime), energySum};
  }
  return std::nullopt;
}

//------------------------------------------------------------------------------

BaseDigitiser::IntegrationResult BaseDigitiser::rocIntegration(const edm4hep::SimCalorimeterHit& hit) const {
  // Collect the MC contributions, then sort them by time
  std::vector<MCC> mcconts;
  mcconts.reserve(hit.contributions_size());
  for (const auto& contribution : hit.getContributions()) {
    mcconts.push_back({contribution.getEnergy(), contribution.getTime()});
  }
  if (mcconts.empty()) {
    return std::nullopt;
  }
  const std::size_t ncontrib = mcconts.size();
  std::sort(mcconts.begin(), mcconts.end(), [](const MCC& lhs, const MCC& rhs) { return (lhs.time < rhs.time); });
  // Accumulate energy until threshold is reached.
  // The first MC contriubtion after the threshold has been reached sets the hit time
  bool passThreshold = false;
  float epar = 0.f, hitTime = 0.f;
  std::size_t thresholdIndex = 0;
  // First determine the hit time (hitTime) and the initial hit index
  // at which we need to start the integration (thresholdIndex)
  for (std::size_t i = 0; i < ncontrib; ++i) {
    const auto timei = mcconts[i].time;
    thresholdIndex = i;
    epar = 0.f;
    for (std::size_t j = i; j < ncontrib; ++j) {
      const auto timej = mcconts[j].time;
      if ((timej - timei) < m_fast_shaper) {
        epar += mcconts[j].energy;
      } else {
        break;
      }
      if (convertEnergy(epar, EnergyScale::GEVDEP) > m_threshold_value) {
        hitTime = timej;
        passThreshold = true;
        break;
      }
    }
    if (passThreshold) {
      break;
    }
  }
  // check hit time
  const float thresholdTime = mcconts[thresholdIndex].time;
  if (not(thresholdTime > m_time_windowMin && thresholdTime < m_time_windowMax)) {
    return std::nullopt;
  }
  // If we've found a hit above the threshold, accumulate the energy until
  // until the maximum time given by the slow shaper
  if (passThreshold) {
    float energySum = 0.f;
    for (std::size_t i = thresholdIndex; i < ncontrib; ++i) {
      if (mcconts[i].time < thresholdTime + m_slow_shaper) {
        energySum += mcconts[i].energy;
      }
    }
    hitTime = smearTime(hitTime);
    return TimeEnergy{hitTime, energySum};
  }
  // else no hit dude !
  else {
    return std::nullopt;
  }
}

//------------------------------------------------------------------------------

float BaseDigitiser::smearTime(float time) const {
  return m_time_resol > 0.f ? time + m_engine.Gaus(0, m_time_resol) : time;
}
