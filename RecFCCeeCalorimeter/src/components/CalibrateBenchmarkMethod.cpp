#include "CalibrateBenchmarkMethod.h"

// Key4HEP
#include "k4Interface/IGeoSvc.h"

// datamodel
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DD4hep/Readout.h"

#include "GaudiKernel/ITHistSvc.h"

#include "TH1F.h"
#include "TMath.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

/// vectors containing the energy deposits to be used for minimization
std::vector<double> vecEgenerated;
std::vector<double> vecEinECaltotal;
std::vector<double> vecEinHCaltotal;
std::vector<double> vecEinHCalfirstLayer;
std::vector<double> vecEinECallastLayer;

DECLARE_COMPONENT(CalibrateBenchmarkMethod)


CalibrateBenchmarkMethod::CalibrateBenchmarkMethod(const std::string& aName, ISvcLocator* aSvcLoc)
    : GaudiAlgorithm(aName, aSvcLoc), 
    m_geoSvc("GeoSvc", aName), 
    m_histSvc("THistSvc", aName),
    m_totalEnergyECal(nullptr),
    m_totalEnergyHCal(nullptr),
    m_totalEnergyBoth(nullptr),
    m_parameters(nullptr)
    {
      declareProperty("ecalBarrelCells", m_ecalBarrelCells, "");
      declareProperty("hcalBarrelCells", m_hcalBarrelCells, "");
    }


CalibrateBenchmarkMethod::~CalibrateBenchmarkMethod() {}


StatusCode CalibrateBenchmarkMethod::initialize() {
  if (GaudiAlgorithm::initialize().isFailure()) return StatusCode::FAILURE;

  // Check geometry service
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service! "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }

  // Check if readouts exist
  {
    bool readoutMissing = false;
    for (size_t i = 0; i < m_readoutNames.size(); ++i) {
      auto readouts = m_geoSvc->lcdd()->readouts();
      if (readouts.find(m_readoutNames.value().at(i)) == readouts.end()) {
        readoutMissing = true;
      }
    }
    if (readoutMissing) {
      error() << "Missing readout, exiting!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // book histograms for checks
  m_totalEnergyECal = new TH1F("totalEnergyECal", "Total deposited energy in ECal", 1000, 0, 1.5 * m_energy);
  if (m_histSvc->regHist("/rec/ecal_total", m_totalEnergyECal).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }

  m_totalEnergyHCal = new TH1F("totalEnergyHCal", "Total deposited energy in HCal", 1000, 0, 1.5 * m_energy);
  if (m_histSvc->regHist("/rec/hcal_total", m_totalEnergyHCal).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }

  m_totalEnergyBoth = new TH1F("totalEnergyBoth", "Total deposited energy in ECal+HCal", 1000, 0, 1.5 * m_energy);
  if (m_histSvc->regHist("/rec/both_total", m_totalEnergyBoth).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }

  m_parameters = new TH1F("parameters","", 4, 0,4);
  if (m_histSvc->regHist("/rec/parameters", m_parameters).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }

  // clear vectors that will be later used for fitting
  vecEgenerated.clear();
  vecEinECaltotal.clear();
  vecEinHCaltotal.clear();
  vecEinHCalfirstLayer.clear();
  vecEinECallastLayer.clear();


  m_energyInECalLayer.resize(m_numLayersECal);
  m_energyInHCalLayer.resize(m_numLayersHCal);

  return StatusCode::SUCCESS;
}


StatusCode CalibrateBenchmarkMethod::execute() {
  double energyInFirstHCalLayer = 0; 
  double energyInLastECalLayer = 0; 
  double energyInBoth = 0.;

  int ecal_index;
  int hcal_index;

  for (double& eneEcal : m_energyInECalLayer){
    eneEcal = 0.;
  }

  for (double& eneHcal : m_energyInHCalLayer){
    eneHcal = 0.;
  }

  // identify ECal and HCal readout position in the input parameters when running the calibration
  for (uint j=0; j<m_readoutNames.size(); j++){
    if (m_SystemID[j]==m_ECalSystemID){
      ecal_index = j; 
    }
    else if (m_SystemID[j]==m_HCalSystemID){
      hcal_index = j; 
    }
  }

  // decoders for ECal and HCal
  auto decoder_ECal = m_geoSvc->lcdd()->readout(m_readoutNames[ecal_index]).idSpec().decoder();
  auto decoder_HCal = m_geoSvc->lcdd()->readout(m_readoutNames[hcal_index]).idSpec().decoder();

  std::vector<size_t> cellIDs; 

  // Get the calorimeter hit collection 
  const edm4hep::CalorimeterHitCollection* ecalBarrelCells = m_ecalBarrelCells.get();
  verbose() << "Input Ecal barrel cell collection size: " << ecalBarrelCells->size() << endmsg;
  const edm4hep::CalorimeterHitCollection* hcalBarrelCells = m_hcalBarrelCells.get();
  verbose() << "Input hadronic barrel cell collection size: " << hcalBarrelCells->size() << endmsg;


  // Loop over a collection of ECal cells and get energy in each ECal layer
  for (const auto& iCell : *ecalBarrelCells) {
    dd4hep::DDSegmentation::CellID cellID = iCell.getCellID(); 
    size_t layerIDecal = decoder_ECal->get(cellID, "layer");
    m_energyInECalLayer.at(layerIDecal) += iCell.getEnergy();
  } 

  // Loop over a collection of HCal cells and get energy in each HCal layer
  for (const auto& iCell : *hcalBarrelCells) {
    dd4hep::DDSegmentation::CellID cellID = iCell.getCellID(); 
    size_t layerIDhcal = decoder_HCal->get(cellID, "layer");
    m_energyInHCalLayer.at(layerIDhcal) += iCell.getEnergy();
  }

  // Energy deposited in the whole ECal
  const double energyInECal = std::accumulate(m_energyInECalLayer.begin(), m_energyInECalLayer.end(), 0);

  // Energy deposited in the last ECal layer
  energyInLastECalLayer = m_energyInECalLayer.at(m_numLayersECal-1);

  // Energy deposited in the whole HCal
  const double energyInHCal = std::accumulate(m_energyInHCalLayer.begin(), m_energyInHCalLayer.end(), 0);

  // Energy deposited in the first HCal layer
  energyInFirstHCalLayer = m_energyInHCalLayer.at(m_firstLayerHCal);

  // Total energy deposited in ECal and HCal
  energyInBoth = energyInECal+energyInHCal;

  // Fill histograms with ECal and HCal energy
  m_totalEnergyECal->Fill(energyInECal);
  m_totalEnergyHCal->Fill(energyInHCal);
  m_totalEnergyBoth->Fill(energyInBoth);
  

  // Fill vectors that will be later used for the energy fitting - length of the vector = number of events
  vecEgenerated.push_back(m_energy);
  vecEinECaltotal.push_back(energyInECal); 
  vecEinHCaltotal.push_back(energyInHCal); 
  vecEinHCalfirstLayer.push_back(energyInFirstHCalLayer);
  vecEinECallastLayer.push_back(energyInLastECalLayer);

  // Printouts for checks
  verbose() << "********************************************************************" << endmsg;
  verbose() << "Energy in ECAL and HCAL: " << energyInECal << " GeV" << endmsg;
  verbose() << "Energy in ECAL: " << energyInECal << " GeV" << endmsg;
  verbose() << "Energy in HCAL: " << energyInHCal << " GeV" << endmsg;
  verbose() << "Energy in ECAL last layer: " << energyInLastECalLayer << " GeV" << endmsg;
  verbose() << "Energy in HCAL first layer: " << energyInFirstHCalLayer << " GeV" << endmsg;
 
  return StatusCode::SUCCESS;
}


// minimisation function for the benchmark method
Double_t chiSquareFitBarrel(const Double_t *par) {
  Double_t fitvalue = 0.;
  // loop over all events, vector of a size of #evts is filled with the energies 
  // ECal calibrated to EM scale, HCal calibrated to HAD scale
  for(uint i=0; i<vecEinECaltotal.size(); i++){
    double E_generated = vecEgenerated.at(i);
    double E_ECal_total = vecEinECaltotal.at(i);
    double E_ECal_lastLayer = vecEinECallastLayer.at(i);
    double E_HCal_total = vecEinHCaltotal.at(i);
    double E_HCal_firstLayer = vecEinHCalfirstLayer.at(i);
    
    Double_t E0 = E_ECal_total*par[0] + E_HCal_total*par[1] + par[2]*sqrt(abs(E_ECal_lastLayer*par[0]*E_HCal_firstLayer*par[1])) + par[3]*pow(E_ECal_total*par[0],2);
    
    fitvalue += pow((E_generated-E0),2)/E_generated; 
  }
  return fitvalue;    
}


StatusCode CalibrateBenchmarkMethod::finalize() {
  // the actual minimisation is running here 
  std::cout << "Running minimisation for " << vecEgenerated.size() << " #events! \n";
  if (vecEgenerated.size() != vecEinECaltotal.size()){
    std::cout << "Something's wrong!!" << std::endl; 
  } 
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  
  // set tolerance 
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
  min->SetMaxIterations(10000);  // for GSL 
  min->SetTolerance(0.001);
  min->SetPrintLevel(1);
    
  // create funciton wrapper for minimizer
  int n_param = 4;
  ROOT::Math::Functor f(&chiSquareFitBarrel,n_param); 
  static const std::vector<double> steps = {0.001, 0.001, 0.001, 0.001};
  static const std::vector<double>  variable = {1., 1., .5, .01};
  min->SetFunction(f);
    
  // Variables to be minimized 
  min->SetVariable(0,"p0", variable[0], steps[0] ); 
  min->SetVariable(1,"p1", variable[1], steps[1] ); 
  min->SetVariable(2,"p2", variable[2], steps[2] ); 
  min->SetVariable(3,"p3", variable[3], steps[3] );

  // fix par1 because HCal is already calibrated to HAD scale
  min->FixVariable(1); 
  //  min->SetVariableLimits(2, 0, 1e6);
   
  // do the minimization
  min->Minimize(); 
    
  const Double_t *xs = min->X();
  const Double_t *ys = min->Errors();
  std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "," << xs[2] <<  "," << xs[3] << "): "  << min->MinValue()  << std::endl;
  std::cout << "Minimum: #delta f(" << ys[0] << "," << ys[1] << "," << ys[2] << "," << ys[3] <<std::endl;

  m_parameters->SetBinContent(1, xs[0]);
  m_parameters->SetBinContent(2, xs[1]);
  m_parameters->SetBinContent(3, xs[2]);
  m_parameters->SetBinContent(4, xs[3]);

  m_parameters->SetBinError(1, ys[0]);
  m_parameters->SetBinError(2, ys[1]);
  m_parameters->SetBinError(3, ys[2]);
  m_parameters->SetBinError(4, ys[3]);

  return GaudiAlgorithm::finalize();
}
