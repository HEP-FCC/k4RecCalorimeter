#ifndef REALISTICCALODIGI_H
#define REALISTICCALODIGI_H 1

#include <k4FWCore/Transformer.h>
#include <edm4hep/SimCalorimeterHit.h>
#include <edm4hep/SimCalorimeterHitCollection.h>
#include <edm4hep/CalorimeterHitCollection.h>
#include <edm4hep/EventHeaderCollection.h>
#include <edm4hep/CaloHitSimCaloHitLinkCollection.h>

#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

#include "TRandom2.h"

#include <string>
#include <tuple>
#include <vector>
#include <functional>
#include <optional>


/** === RealisticCaloDigi Processor === <br>
    Digitisation of calorimeter hits
    e.g. timing, dead cells, miscalibrations
    this is virtual class, technology-blind
    technology-specific classes can inherit from this one
    D. Jeans 02/2016, rewrite of parts of ILDCaloDigi, DDCaloDigi
    R. Ete 11/2020, rewrite of charge integration and extension of timing treatment
 */

struct RealisticCaloDigi : k4FWCore::MultiTransformer<
    std::tuple<edm4hep::CalorimeterHitCollection, 
    edm4hep::CaloHitSimCaloHitLinkCollection>(
      const edm4hep::SimCalorimeterHitCollection&,
      const edm4hep::EventHeaderCollection&)> {
  public:
    RealisticCaloDigi(const std::string& name, ISvcLocator* svcLoc);
    /** Called at the begin of the job before anything is read.
    * Use to initialize the processor, e.g. book histograms.
    */
    StatusCode initialize();
  
    /** Called for every run.
     */
    std::tuple<edm4hep::CalorimeterHitCollection, 
               edm4hep::CaloHitSimCaloHitLinkCollection> operator()(
              const edm4hep::SimCalorimeterHitCollection& inputSim,
              const edm4hep::EventHeaderCollection& headers) const; 
  
    /** Called after data processing for clean up.
     */
    StatusCode finalize();


  protected:
  
    // energy scales we know about
    enum { MIP, GEVDEP, NPE };
    // integration result types
    using integr_res = std::pair<float,float>;
    using integr_res_opt = std::optional<integr_res>;
    using integr_function = std::function<integr_res_opt(const edm4hep::SimCalorimeterHit*)>;

   virtual float EnergyDigi(float energy, float event_correl_miscalib) const;
   virtual integr_res_opt Integrate( const edm4hep::SimCalorimeterHit * hit ) const;
   
   integr_res_opt StandardIntegration( const edm4hep::SimCalorimeterHit * hit ) const ;
   integr_res_opt ROCIntegration( const edm4hep::SimCalorimeterHit * hit ) const ;
   float SmearTime(float time) const;

   // virtual methods to be be overloaded in tech-specific derived classes
   virtual int    getMyUnit() const = 0 ;
   virtual float digitiseDetectorEnergy(float energy) const = 0 ;
   virtual float convertEnergy( float energy, int inScale ) const = 0; // convert energy from input to output scale

   // timing
   Gaudi::Property<int> m_time_apply{this, "timingCut", 0, "Use hit times"};
   Gaudi::Property<int> m_time_correctForPropagation{this, "timingCorrectForPropagation", 0, "Correct hit times for propagation: radial distance/c"};
   Gaudi::Property<float> m_time_windowMin{this, "timingWindowMin", -10.0f, "Time Window minimum time in ns"};
   Gaudi::Property<float> m_time_windowMax{this, "timingWindowMax", 100.0f, "Time Window maximum time in ns"};
   Gaudi::Property<std::string> m_integration_method{this, "integrationMethod", "Standard", "Energy integration and time calculation method. Options: Standard, ROC"};
   Gaudi::Property<float> m_fast_shaper{this, "fastShaper", 0.f, "Fast shaper value. Unit in ns"};
   Gaudi::Property<float> m_slow_shaper{this, "slowShaper", 0.f, "Slow shaper value. Unit in ns"};
   Gaudi::Property<float> m_time_resol{this, "timingResolution", 0.f, "Time resolution to apply (gaussian smearing). Unit in ns"};
   // additional digi effects
   Gaudi::Property<float> m_calib_mip{this, "calibration_mip", 1.0e-4f, "Average G4 deposited energy by MIP for calibration"};
   Gaudi::Property<float> m_misCalib_uncorrel{this, "miscalibration_uncorrel", 0.0f, "Uncorrelated random Gaussian miscalibration (as a fraction: 1.0 = 100%)"};
   Gaudi::Property<float> m_misCalib_correl{this, "miscalibration_correl", 0.0f, "Correlated random Gaussian miscalibration (as a fraction: 1.0 = 100%)"};
   Gaudi::Property<float> m_deadCell_fraction{this, "deadCell_fraction", 0.0f, "Random dead cell fraction (as a fraction: 0->1)"};
   // simple model of electronics properties
   Gaudi::Property<float> m_elec_noiseMip{this, "elec_noise_mip", 0.0f, "Typical electronics noise (in MIP units)"};
   Gaudi::Property<float> m_elec_rangeMip{this, "elec_range_mip", 2500.0f, "Maximum of dynamic range of electronics (in MIPs)"};
   // code for layer info for cellID decoder
   Gaudi::Property<std::string> m_encodingStringVariable{this, "EncodingStringParameterName", "GlobalTrackerReadoutID", "The name of the DD4hep constant that contains the Encoding string for tracking detectors"};
   // energy threshold
   Gaudi::Property<float> m_threshold_value{this, "threshold", 0.5f, "Threshold for Hit"};
   Gaudi::Property<std::string> m_threshold_unit{this, "thresholdUnit", std::string("MIP"), "Unit for threshold. Can be \"GeV\", \"MIP\" or \"px\". MIP and px need properly set calibration constants"};
   // id parameters
   Gaudi::Property<std::string> m_calo_type {this, "CaloType", "em", "Calorimeter Type: em, had, mu"};
   Gaudi::Property<std::string> m_calo_id {this, "CaloID", "ecal", "Calorimeter ID: ecal, hcal, yoke, lcal, lhcal, bcal"};
   Gaudi::Property<std::string> m_calo_layout {this, "CaloLayout", "barrel", "Calorimeter Layout: barrel, endcap, ring, plug"};
   

   int m_threshold_iunit{};  
   inline static thread_local TRandom2 m_engine;
   SmartIF<IGeoSvc>                    m_geoSvc;
   SmartIF<IUniqueIDGenSvc>            m_uidSvc;
 

   integr_function m_integr_function{};

} ;

#endif



