//==========================================================================
//  AIDA Detector description implementation 
//--------------------------------------------------------------------------
// Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
//
// Author     : M.Frank
//
//==========================================================================
#ifndef EXAMPLES_DDDualCrys_SRC_DualCrysCalorimeterHIT_H
#define EXAMPLES_DDDualCrys_SRC_DualCrysCalorimeterHIT_H

/// Framework include files
#include "DDG4/Geant4Data.h"
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"


typedef ROOT::Math::XYZVector Position;
typedef ROOT::Math::XYZVector Direction;



namespace CalVision {


  const int finenbin=40;
  const int coarsenbin=4;


    //const int truthnbin=1000;

  /// This is the hit definition.
  /** I took here the same definition of the default Geant4Tracker class,
   *  (see DDG4/Geant4Data.h)  but it could be anything else as well.
   *
   *  Please note:
   *  ============
   *  The MC truth handling as implemented in the Geant4ParticleHandler
   *  will not work with this class if the object(s) are saved with 
   *  the standard Geant4Output2ROOT event action. If the hit is 
   *  specialized, the output writing also must be specialized if
   *  MC truth handling should be supported.
   *  Otherwise it is sufficient to provide a ROOT dictionary as long as the
   *  base class dd4hep::sim::Geant4HitData is kept.
   *
   *  \author  M.Frank
   *  \version 1.0
   *  \ingroup DD4HEP_SIMULATION
   */
  class DualCrysCalorimeterHit : public dd4hep::sim::Geant4Calorimeter::Hit   {

  public:
    int n_inelastic;
    int ncerenkov,nscintillator;
    float edeprelativistic;
    float edepepgam;
    float wavelenmin=300;
    float wavelenmax=1000;


    int nfinebin=finenbin;
    float timemin=0;
    float timemax=400;
    float timemaxz=40;
    //    std::array<int,finenbin>  ncerwave;
    // std::array<int,finenbin> nscintwave;
    //std::array<int,finenbin>  ncertime;
    //std::array<int,finenbin> nscinttime;
    //std::array<int,finenbin>  ncertimez;
    //std::array<int,finenbin> nscinttimez;
    //std::array<float,finenbin> edeptime;
    //std::array<float,finenbin> ereldeptime;
    //float xmax=10;
    //float ymax=10;
    //float xmin=-10;
    //float ymin=-10;
    //int ncoarsebin=coarsenbin;
    // std::array<std::array<int,coarsenbin>,coarsenbin> cerhitpos;
    //std::array<std::array<int,coarsenbin>,coarsenbin> scinthitpos;

    std::vector<std::pair<float, float>> HitCer;
    std::vector<std::pair<float, float>> HitScin;
    //    std::vector<float> CerTime;
    //std::vector<float> ScinTime;



    //int ntruthbin=truthnbin;
    //std::array<float,truthnbin> contribBeta;
    //std::array<float,truthnbin> contribCharge;
    std::vector<float> contribBeta;
    std::vector<float> contribCharge;


  public:
    /// Default constructor
    DualCrysCalorimeterHit() = default;
    /// Initializing constructor
    DualCrysCalorimeterHit(const Position& cell_pos)
  	: dd4hep::sim::Geant4Calorimeter::Hit(cell_pos)
  	, n_inelastic(0)          // first int member declared
  	, ncerenkov(0)
  	, nscintillator(0)
  	, edeprelativistic(0.)    // then float members in declared order
  	, edepepgam(0.)
  	, wavelenmin(300)
  	, wavelenmax(1000)
    {
      for (int i = 0; i < finenbin; i++) {
        //ncerwave[i]=0;
        //nscintwave[i]=0;
        //ncertime[i]=0;
        //nscinttime[i]=0;
        //ncertimez[i]=0;
        //nscinttimez[i]=0;
	//edeptime[i]=0.;
	//ereldeptime[i]=0.;
      }
      //      for( int i=0;i<coarsenbin;i++ ) {
      //  for( int j=0;j<coarsenbin;j++ ) {
      //    cerhitpos[i][j]=0;
      //    scinthitpos[i][j]=0;
      //  }
      // }

}

    /// Default destructor
    virtual ~DualCrysCalorimeterHit() = default;
    /// Assignment operator
    //DualCrysCalorimeterHit& operator=(const DualCrysCalorimeterHit& c);
  };

  /// Helper to dump data file
  /**
   *  Usage:  
   *  $> root.exe
   *  ....
   *  root [0] gSystem->Load("libDDG4Plugins.so");
   *  root [1] gSystem->Load("libDDG4_MySensDet.so");
   *  root [2] CalVision::Dump::dumpData(<num-ebents>,<file-name>);
   *
   */
  class Dump   {
  public:
    /// Standalone function to dump data from a root file
    static int DualCrysCalorimeterdumpData(int num_evts, const char* file_name);
  };
}

// CINT configuration
#if defined(__CINT__) || defined(__MAKECINT__) || defined(__CLING__) || defined(__ROOTCLING__)
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

/// Define namespaces
#pragma link C++ namespace dd4hep;
#pragma link C++ namespace dd4hep::sim;
#pragma link C++ namespace CalVision;
#pragma link C++ class     CalVision::DualCrysCalorimeterHit+;
#pragma link C++ class     CalVision::DualCrysCalorimeterDump;
#endif

#endif // EXAMPLES_DDDualCrys_SRC_DualCrysCalorimeterHIT_H

