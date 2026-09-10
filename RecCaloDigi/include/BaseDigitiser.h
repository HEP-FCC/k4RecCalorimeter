#ifndef BASEDIGITISER_H
#define BASEDIGITISER_H 1

#include <edm4hep/CaloHitSimCaloHitLinkCollection.h>
#include <edm4hep/CalorimeterHitCollection.h>
#include <edm4hep/EventHeaderCollection.h>
#include <edm4hep/SimCalorimeterHit.h>
#include <edm4hep/SimCalorimeterHitCollection.h>
#include <k4FWCore/Transformer.h>

#include "CalorimeterHitType.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

#include <DDSegmentation/BitFieldCoder.h>

#include "TRandom3.h"

#include <functional>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

/** === BaseDigitiser algorithm === <br>
    Digitisation of calorimeter hits
    e.g. timing, dead cells, miscalibrations
    this is an abstract base class, technology-blind
    technology-specific classes inherit from this one
    D. Jeans 02/2016, rewrite of parts of ILDCaloDigi, DDCaloDigi
    R. Ete 11/2020, rewrite of charge integration and extension of timing treatment
 */

struct BaseDigitiser : k4FWCore::MultiTransformer<
                           std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>(
                               const edm4hep::SimCalorimeterHitCollection&, const edm4hep::EventHeaderCollection&)> {
public:
  BaseDigitiser(const std::string& name, ISvcLocator* svcLoc);
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  StatusCode initialize();

  /** Called for every event.
   */
  std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>
  operator()(const edm4hep::SimCalorimeterHitCollection& inputSim, const edm4hep::EventHeaderCollection& headers) const;

protected:
  // energy scales we know about
  enum class EnergyScale { MIP, GEVDEP, NPE };
  // result of integrating the contributions of one simulated hit
  struct TimeEnergy {
    float time{0.f};
    float energy{0.f};
  };
  using IntegrationResult = std::optional<TimeEnergy>;
  using IntegrationFunction = std::function<IntegrationResult(const edm4hep::SimCalorimeterHit&)>;

  float energyDigi(float energy, float event_correl_miscalib) const;
  IntegrationResult integrate(const edm4hep::SimCalorimeterHit& hit) const;

  IntegrationResult standardIntegration(const edm4hep::SimCalorimeterHit& hit) const;
  IntegrationResult rocIntegration(const edm4hep::SimCalorimeterHit& hit) const;
  float smearTime(float time) const;

  // virtual methods to be be overloaded in tech-specific derived classes
  virtual float digitiseDetectorEnergy(float energy) const = 0;
  // convert energy from the input scale to the scale of the derived class
  virtual float convertEnergy(float energy, EnergyScale inScale) const = 0;
  // whether convertEnergy() accepts that input scale. Checked in initialize() so that a
  // threshold configured in units this technology cannot express fails at startup, rather
  // than at the point of conversion.
  virtual bool canConvertFrom(EnergyScale inScale) const = 0;

  // timing
  Gaudi::Property<int> m_time_apply{this, "timingCut", 0, "Use hit times"};
  Gaudi::Property<int> m_time_correctForPropagation{this, "timingCorrectForPropagation", 0,
                                                    "Correct hit times for propagation: radial distance/c"};
  Gaudi::Property<float> m_time_windowMin{this, "timingWindowMin", -10.0f, "Time Window minimum time in ns"};
  Gaudi::Property<float> m_time_windowMax{this, "timingWindowMax", 100.0f, "Time Window maximum time in ns"};
  Gaudi::Property<std::string> m_integration_method{
      this, "integrationMethod", "Standard", "Energy integration and time calculation method. Options: Standard, ROC"};
  Gaudi::Property<float> m_fast_shaper{this, "fastShaper", 0.f, "Fast shaper value. Unit in ns"};
  Gaudi::Property<float> m_slow_shaper{this, "slowShaper", 0.f, "Slow shaper value. Unit in ns"};
  Gaudi::Property<float> m_time_resol{this, "timingResolution", 0.f,
                                      "Time resolution to apply (gaussian smearing). Unit in ns"};
  // additional digi effects
  Gaudi::Property<float> m_calib_mip{this, "calibration_mip", 1.0e-4f,
                                     "Average G4 deposited energy by MIP for calibration"};
  Gaudi::Property<float> m_misCalib_uncorrel{this, "miscalibration_uncorrel", 0.0f,
                                             "Uncorrelated random Gaussian miscalibration (as a fraction: 1.0 = 100%)"};
  Gaudi::Property<float> m_misCalib_correl{this, "miscalibration_correl", 0.0f,
                                           "Correlated random Gaussian miscalibration (as a fraction: 1.0 = 100%)"};
  Gaudi::Property<float> m_deadCell_fraction{this, "deadCell_fraction", 0.0f,
                                             "Random dead cell fraction (as a fraction: 0->1)"};
  // simple model of electronics properties
  Gaudi::Property<float> m_elec_noiseMip{this, "elec_noise_mip", 0.0f, "Typical electronics noise (in MIP units)"};
  Gaudi::Property<float> m_elec_rangeMip{this, "elec_range_mip", 2500.0f,
                                         "Maximum of dynamic range of electronics (in MIPs)"};
  // name of the DD4hep constant holding the calorimeter cellID encoding
  Gaudi::Property<std::string> m_encodingStringVariable{
      this, "EncodingStringParameterName", "GlobalCalorimeterReadoutID",
      "The name of the DD4hep constant that contains the encoding string for calorimeters"};
  // energy threshold
  Gaudi::Property<float> m_threshold_value{this, "threshold", 0.5f, "Threshold for Hit"};
  Gaudi::Property<std::string> m_threshold_unit{
      this, "thresholdUnit", std::string("MIP"),
      "Unit for threshold. Can be \"GeV\", \"MIP\" or \"px\". MIP and px need properly set calibration constants"};
  // id parameters
  Gaudi::Property<std::string> m_calo_type{this, "CaloType", "em", "Calorimeter Type: em, had, muon"};
  Gaudi::Property<std::string> m_calo_id{this, "CaloID", "ecal", "Calorimeter ID: ecal, hcal, yoke, lcal, lhcal, bcal"};
  Gaudi::Property<std::string> m_calo_layout{this, "CaloLayout", "barrel",
                                             "Calorimeter Layout: barrel, endcap, ring, plug"};

  EnergyScale m_threshold_iunit{EnergyScale::MIP};
  // decoded once in initialize(): the calorimeter type, ID and layout are the same
  // for every hit of every event
  CHT::CaloType m_chtType{CHT::CaloType::unknown};
  CHT::CaloID m_chtId{CHT::CaloID::unknown};
  CHT::Layout m_chtLayout{CHT::Layout::any};
  // Re-seeded from the UniqueIDGenSvc at the start of every event, so results depend only on
  // (run, event, algorithm name) and not on how events are distributed over threads.
  inline static thread_local TRandom3 m_engine;
  SmartIF<IGeoSvc> m_geoSvc;
  SmartIF<IUniqueIDGenSvc> m_uidSvc;
  // built once in initialize() from the geometry encoding string
  dd4hep::DDSegmentation::BitFieldCoder m_bitFieldCoder{};

  IntegrationFunction m_integr_function{};
};

#endif
