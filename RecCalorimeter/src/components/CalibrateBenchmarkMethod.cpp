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

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "TH1F.h"
#include "TMath.h"

// Include the <cmath> header for std::fabs
#include <cmath>

DECLARE_COMPONENT(CalibrateBenchmarkMethod)

CalibrateBenchmarkMethod::CalibrateBenchmarkMethod(const std::string& aName, ISvcLocator* aSvcLoc)
    : Gaudi::Algorithm(aName, aSvcLoc), m_geoSvc("GeoSvc", aName), m_histSvc("THistSvc", aName),
      m_totalEnergyECal(nullptr), m_totalEnergyHCal(nullptr), m_totalEnergyBoth(nullptr), m_parameters(nullptr) {
  declareProperty("ecalBarrelCells", m_ecalBarrelCells, "");
  declareProperty("hcalBarrelCells", m_hcalBarrelCells, "");
}

CalibrateBenchmarkMethod::~CalibrateBenchmarkMethod() {}

StatusCode CalibrateBenchmarkMethod::initialize() {
  if (Gaudi::Algorithm::initialize().isFailure())
    return StatusCode::FAILURE;

  // Check geometry service
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service! "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }
  // Check histogram service
  if (!m_histSvc) {
    error() << "Unable to locate Histogram Service" << endmsg;
    return StatusCode::FAILURE;
  }

  // Check if readouts exist
  {
    bool readoutMissing = false;
    for (size_t i = 0; i < m_readoutNames.size(); ++i) {
      auto readouts = m_geoSvc->getDetector()->readouts();
      if (readouts.find(m_readoutNames.value().at(i)) == readouts.end()) {
        readoutMissing = true;
      }
    }
    if (readoutMissing) {
      error() << "Missing readout, exiting!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // number of benchmark parameters
  int n_param = 6;

  // book histograms
  m_totalEnergyECal = new TH1F("totalEnergyECal", "Total deposited energy in ECal", 1000, 0, 1.5 * m_energy);
  m_totalEnergyHCal = new TH1F("totalEnergyHCal", "Total deposited energy in HCal", 1000, 0, 1.5 * m_energy);
  m_totalEnergyBoth = new TH1F("totalEnergyBoth", "Total deposited energy in ECal+HCal", 1000, 0, 1.5 * m_energy);
  m_parameters = new TH1F("parameters", "Benchmark parameters", n_param, 0, n_param);

  registerHistogram("/rec/ecal_total", m_totalEnergyECal);
  registerHistogram("/rec/hcal_total", m_totalEnergyHCal);
  registerHistogram("/rec/both_total", m_totalEnergyBoth);
  registerHistogram("/rec/parameters", m_parameters);

  // clear vectors that will be later used for fitting
  m_vecGeneratedEnergy.clear();
  m_vecTotalEnergyinECal.clear();
  m_vecTotalEnergyinHCal.clear();
  m_vecEnergyInFirstLayerHCal.clear();
  m_vecEnergyInLastLayerECal.clear();
  m_vecEnergyInFirstLayerECal.clear();

  m_energyInLayerECal.resize(m_numLayersECal);
  m_energyInLayerHCal.resize(m_numLayersHCal);

  return StatusCode::SUCCESS;
}

StatusCode CalibrateBenchmarkMethod::execute(const EventContext&) const {
  double totalEnergyInECal = 0.;
  double totalEnergyInHCal = 0.;
  double energyInFirstLayerECal = 0;
  double energyInLastLayerECal = 0;
  double energyInFirstLayerHCal = 0;
  double energyInBoth = 0.;

  int indexECal = -1;
  int indexHCal = -1;

  for (double& eneECal : m_energyInLayerECal) {
    eneECal = 0.;
  }

  for (double& eneHCal : m_energyInLayerHCal) {
    eneHCal = 0.;
  }

  // identify ECal and HCal readout position in the input parameters when running the calibration
  for (uint j = 0; j < m_readoutNames.size(); j++) {
    if (m_systemID[j] == m_systemIDECal) {
      indexECal = j;
    } else if (m_systemID[j] == m_systemIDHCal) {
      indexHCal = j;
    }
  }

  // decoders for ECal and HCal
  auto decoderECal = m_geoSvc->getDetector()->readout(m_readoutNames[indexECal]).idSpec().decoder();
  auto decoderHCal = m_geoSvc->getDetector()->readout(m_readoutNames[indexHCal]).idSpec().decoder();

  std::vector<size_t> cellIDs;

  // Get the calorimeter hit collection
  const edm4hep::CalorimeterHitCollection* ecalBarrelCells = m_ecalBarrelCells.get();
  verbose() << "Input Ecal barrel cell collection size: " << ecalBarrelCells->size() << endmsg;
  const edm4hep::CalorimeterHitCollection* hcalBarrelCells = m_hcalBarrelCells.get();
  verbose() << "Input hadronic barrel cell collection size: " << hcalBarrelCells->size() << endmsg;

  // Loop over a collection of ECal cells and get energy in each ECal layer
  for (const auto& iCell : *ecalBarrelCells) {
    dd4hep::DDSegmentation::CellID cellID = iCell.getCellID();
    size_t layerIDecal = decoderECal->get(cellID, "layer");
    m_energyInLayerECal.at(layerIDecal) += iCell.getEnergy();
  }

  // Loop over a collection of HCal cells and get energy in each HCal layer
  for (const auto& iCell : *hcalBarrelCells) {
    dd4hep::DDSegmentation::CellID cellID = iCell.getCellID();
    size_t layerIDhcal = decoderHCal->get(cellID, "layer");
    m_energyInLayerHCal.at(layerIDhcal) += iCell.getEnergy();
  }

  // Energy deposited in the whole ECal
  for (size_t i = 0; i < m_energyInLayerECal.size(); ++i) {
    totalEnergyInECal += m_energyInLayerECal.at(i);
  }

  // Energy deposited in the first/last ECal layer
  energyInFirstLayerECal = m_energyInLayerECal.at(m_firstLayerECal);
  energyInLastLayerECal = m_energyInLayerECal.at(m_numLayersECal - 1);

  // Energy deposited in the whole HCal
  for (size_t i = 0; i < m_energyInLayerHCal.size(); ++i) {
    totalEnergyInHCal += m_energyInLayerHCal.at(i);
  }

  // Energy deposited in the first HCal layer
  energyInFirstLayerHCal = m_energyInLayerHCal.at(m_firstLayerHCal);

  // Total energy deposited in ECal and HCal
  energyInBoth = totalEnergyInECal + totalEnergyInHCal;

  // Fill histograms with ECal and HCal energy
  m_totalEnergyECal->Fill(totalEnergyInECal);
  m_totalEnergyHCal->Fill(totalEnergyInHCal);
  m_totalEnergyBoth->Fill(energyInBoth);

  // Fill vectors that will be later used for the energy fitting - length of the vector = number of events
  m_vecGeneratedEnergy.push_back(m_energy);
  m_vecTotalEnergyinECal.push_back(totalEnergyInECal);
  m_vecTotalEnergyinHCal.push_back(totalEnergyInHCal);
  m_vecEnergyInFirstLayerECal.push_back(energyInFirstLayerECal);
  m_vecEnergyInLastLayerECal.push_back(energyInLastLayerECal);
  m_vecEnergyInFirstLayerHCal.push_back(energyInFirstLayerHCal);

  // Printouts for checks
  verbose() << "********************************************************************" << endmsg;
  verbose() << "Energy in ECAL and HCAL: " << energyInBoth << " GeV" << endmsg;
  verbose() << "Energy in ECAL: " << totalEnergyInECal << " GeV" << endmsg;
  verbose() << "Energy in HCAL: " << totalEnergyInHCal << " GeV" << endmsg;
  verbose() << "Energy in ECAL last layer: " << energyInLastLayerECal << " GeV" << endmsg;
  verbose() << "Energy in HCAL first layer: " << energyInFirstLayerHCal << " GeV" << endmsg;

  return StatusCode::SUCCESS;
}

// minimisation function for the benchmark method
double CalibrateBenchmarkMethod::chiSquareFitBarrel(const Double_t* parameter) const {
  double fitvalue = 0.;
  // loop over all events, vector of a size of #evts is filled with the energies
  // ECal calibrated to EM scale, HCal calibrated to HAD scale
  for (uint i = 0; i < m_vecTotalEnergyinECal.size(); i++) {
    double generatedEnergy = m_vecGeneratedEnergy.at(i);
    double totalEnergyInECal = m_vecTotalEnergyinECal.at(i);
    double totalEnergyInHCal = m_vecTotalEnergyinHCal.at(i);
    double energyInFirstLayerECal = m_vecEnergyInFirstLayerECal.at(i);
    double energyInLastLayerECal = m_vecEnergyInLastLayerECal.at(i);
    double energyInFirstLayerHCal = m_vecEnergyInFirstLayerHCal.at(i);

    Double_t benchmarkEnergy =
        parameter[0] * totalEnergyInECal + parameter[1] * totalEnergyInHCal +
        parameter[2] *
            std::sqrt(std::fabs(energyInLastLayerECal * parameter[0] * energyInFirstLayerHCal * parameter[1])) +
        parameter[3] * std::pow(totalEnergyInECal * parameter[0], 2) + parameter[4] * energyInFirstLayerECal +
        parameter[5];

    // general formula below gives possibility to fit also parameter[1] to scale E_HCal_total and a contant term
    // parameter[5] to include residuals parameter[0] scales ECal to HAD scale parameter[1] scales HCal to HAD scale
    // (fixed to 1 if HCal is already calibrated at HAD scale at the input level) parameter[2] energy losses between
    // ECal and HCal parameter[3] corrects for non-compensation of ECal parameter[4] upstream correction parameter[5]
    // residuals

    fitvalue += std::pow((generatedEnergy - benchmarkEnergy), 2) / generatedEnergy;
  }
  return fitvalue;
}

StatusCode CalibrateBenchmarkMethod::finalize() {
  // the actual minimisation is running here
  std::cout << "Running minimisation for " << m_vecGeneratedEnergy.size() << " #events! \n";
  if (m_vecGeneratedEnergy.size() != m_vecTotalEnergyinECal.size()) {
    std::cout << "Something's wrong! Vector size does not match for generated energy and total energy in ECal."
              << std::endl;
  }

  // number of benchmark parameters
  int n_param = 6;
  // initial values for the minimization and steps while running the minimization
  static const std::vector<double> variable = {1., 1., .5, .01, 0.5, 0.};
  static const std::vector<double> steps = {0.001, 0.001, 0.001, 0.001, 0.001, 0.001};

  runMinimization(n_param, variable, steps, m_fixedParameters);

  return Gaudi::Algorithm::finalize();
}

void CalibrateBenchmarkMethod::registerHistogram(const std::string& path, TH1F*& histogramName) {
  if (m_histSvc->regHist(path, histogramName).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    throw std::runtime_error("Histogram registration failure");
  }
}

void CalibrateBenchmarkMethod::runMinimization(int n_param, const std::vector<double>& variable,
                                               const std::vector<double>& steps,
                                               const std::vector<int>& fixedParameters) const {
  // Set up the functor (the function to be minimized)
  ROOT::Math::Functor f(this, &CalibrateBenchmarkMethod::chiSquareFitBarrel, n_param);

  // Initialize the minimizer
  ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  minimizer->SetMaxFunctionCalls(1000000);
  minimizer->SetMaxIterations(10000);
  minimizer->SetTolerance(0.001);
  minimizer->SetPrintLevel(1);

  // Set the function to minimize
  minimizer->SetFunction(f);

  // Set parameters for the minimization
  for (int i = 0; i < n_param; ++i) {
    minimizer->SetVariable(i, "p" + std::to_string(i), variable[i], steps[i]);

    // Fix parameters listed in fixedParameters that will not be included in the fit minimization
    if (std::find(fixedParameters.begin(), fixedParameters.end(), i) != fixedParameters.end()) {
      minimizer->FixVariable(i);
    }
  }

  // Perform the minimization
  minimizer->Minimize();

  // Access the results if needed
  const Double_t* xs = minimizer->X();
  const Double_t* ys = minimizer->Errors();

  std::cout << "Minimum: f(";
  for (int i = 0; i < n_param; ++i) {
    std::cout << xs[i];
    if (i < n_param - 1)
      std::cout << ",";
  }
  std::cout << "): " << minimizer->MinValue() << std::endl;

  std::cout << "Minimum: #delta f(";
  for (int i = 0; i < n_param; ++i) {
    std::cout << ys[i];
    if (i < n_param - 1)
      std::cout << ",";
  }
  std::cout << std::endl;

  // histogram->SetBinContent(1, xs[0]);
  // histogram->SetBinError(1, ys[0]);

  // each fitted parameter is stored in one bin of the output histogram
  // bin 1 contains the value for par[0] etc
  for (int i = 0; i < n_param; ++i) {
    verbose() << "inside filling histos: " << xs[i] << " GeV" << endmsg;
    m_parameters->SetBinContent(i + 1, xs[i]);
    m_parameters->SetBinError(i + 1, ys[i]);
  }
}
