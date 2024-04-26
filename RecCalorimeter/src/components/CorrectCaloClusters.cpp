#include "CorrectCaloClusters.h"

// Gaudi
#include "GaudiKernel/ITHistSvc.h"

// Key4HEP
#include "k4Interface/IGeoSvc.h"

// FCC Detectors
#include "detectorCommon/DetUtils_k4geo.h"
#include "detectorSegmentations/FCCSWGridPhiEta_k4geo.h"
#include "detectorSegmentations/FCCSWGridPhiTheta_k4geo.h"
#include "detectorSegmentations/FCCSWGridModuleThetaMerged_k4geo.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DDSegmentation/MultiSegmentation.h"

// our EDM
#include "edm4hep/Cluster.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/VertexCollection.h"

// ROOT
#include "TF2.h"

// Include the <cmath> header for std::fabs
#include <cmath>  

DECLARE_COMPONENT(CorrectCaloClusters)

CorrectCaloClusters::CorrectCaloClusters(const std::string& name,
                                         ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc),
      m_geoSvc("GeoSvc", "CorrectCaloClusters") {
  declareProperty("inClusters", m_inClusters,
                  "Input cluster collection");
  declareProperty("outClusters", m_outClusters,
                  "Corrected (output) cluster collection");
}

StatusCode CorrectCaloClusters::initialize() {
  {
    StatusCode sc = GaudiAlgorithm::initialize();
    if (sc.isFailure()) {
      return sc;
    }
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

  // Check if readout related variables have the same size
  if (m_systemIDs.size() != m_readoutNames.size()) {
    error() << "Sizes of the systemIDs vector and readoutNames vector does not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_numLayers.size()) {
    error() << "Sizes of systemIDs vector and numLayers vector does not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_firstLayerIDs.size()) {
    error() << "Sizes of systemIDs vector and firstLayerIDs vector does not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_upstreamCorr) {
    if (m_systemIDs.size() != m_upstreamFormulas.size()) {
      error() <<  "Sizes of systemIDs vector and upstreamFormulas vector does not match, exiting!" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_systemIDs.size() != m_upstreamParams.size()) {
      error() <<  "Sizes of systemIDs vector and upstreamParams vector does not match, exiting!" << endmsg;
      return StatusCode::FAILURE;
    }
  }
  if (m_downstreamCorr) { 
    if (m_systemIDs.size() != m_downstreamFormulas.size()) {
      error() <<  "Sizes of systemIDs vector and downstreamFormulas vector does not match, exiting!" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_systemIDs.size() != m_downstreamParams.size()) {
      error() <<  "Sizes of systemIDs vector and downstreamParams vector does not match, exiting!" << endmsg;
      return StatusCode::FAILURE;
    }
  }
  if (m_benchmarkCorr) {  
    if (m_systemIDs.size() != m_benchmarkFormulas.size()) {
      error() <<  "Sizes of systemIDs vector and benchmarkFormulas vector does not match, exiting!" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_systemIDs.size() != m_benchmarkParametrization.size()) {
      error() <<  "Sizes of systemIDs vector and benchmarkParams vector does not match, exiting!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // Prepare upstream and downstream correction functions
  {
    StatusCode sc = initializeCorrFunctions(m_upstreamFunctions, m_upstreamFormulas, m_upstreamParams, "upstream");
    if (sc.isFailure()) {
      error() << "Initialization of upstream correction functions not successful!" << endmsg;
      return sc;
    }
  }
  {
    StatusCode sc = initializeCorrFunctions(m_downstreamFunctions, m_downstreamFormulas, m_downstreamParams, "downstream");
    if (sc.isFailure()) {
      error() << "Initialization of downstream correction functions not successful!" << endmsg;
      return sc;
    }
  }
  {
    StatusCode sc = initializeCorrFunctions(m_benchmarkFunctions, m_benchmarkFormulas, m_benchmarkParametrization, "benchmark"); 
    if (sc.isFailure()) {
      error() << "Initialization of benchmark correction functions not successful!" << endmsg;
      return sc;
    }
  }

  info() << "Initialized following upstream correction functions:" << endmsg;
  for (size_t i = 0; i < m_upstreamFunctions.size(); ++i) {
    for (size_t j = 0; j < m_upstreamFunctions[i].size(); ++j) {
      auto func = m_upstreamFunctions.at(i).at(j);
      info() << "  " << func->GetName() << ": " << func->GetExpFormula() << endmsg;
      for (int k = 0; k < func->GetNpar(); ++k) {
        info() << "    " << func->GetParName(k) << ": " << func->GetParameter(k) << endmsg;
      }
    }
  }

  info() << "Initialized following downstream correction functions:" << endmsg;
  for (size_t i = 0; i < m_downstreamFunctions.size(); ++i) {
    for (size_t j = 0; j < m_downstreamFunctions[i].size(); ++j) {
      auto func = m_downstreamFunctions.at(i).at(j);
      info() << "  " << func->GetName() << ": " << func->GetExpFormula() << endmsg;
      for (int k = 0; k < func->GetNpar(); ++k) {
        info() << "    " << func->GetParName(k) << ": " << func->GetParameter(k) << endmsg;
      }
    }
  }

  info() << "Initialized following benchmark correction functions:" << endmsg;
  for (size_t i = 0; i < m_benchmarkFunctions.size(); ++i) {
    for (size_t j = 0; j < m_benchmarkFunctions[i].size(); ++j) {
      auto func = m_benchmarkFunctions.at(i).at(j);
      info() << "  " << func->GetName() << ": " << func->GetExpFormula() << endmsg;
      for (int k = 0; k < func->GetNpar(); ++k) {
        info() << "    " << func->GetParName(k) << ": " << func->GetParameter(k) << endmsg;
      }
    }
  }

  return StatusCode::SUCCESS;
}


StatusCode CorrectCaloClusters::execute() {
  verbose() << "-------------------------------------------" << endmsg;

  // Get the input collection with clusters
  const edm4hep::ClusterCollection* inClusters = m_inClusters.get();

  // Initialize output clusters
  edm4hep::ClusterCollection* outClusters = initializeOutputClusters(inClusters);
  if (!outClusters) {
    error() << "Something went wrong in initialization of the output cluster collection, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (inClusters->size() != outClusters->size()) {
    error() << "Sizes of input and output cluster collections does not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }

  // Apply upstream correction
  if (m_upstreamCorr) {
    StatusCode sc = applyUpstreamCorr(inClusters, outClusters);
    if (sc.isFailure()) {
      return sc;
    }
    else{
      verbose() << "Running the upstream correction." << endmsg;

    }
  }

  // Apply downstream correction
  if (m_downstreamCorr) {
    StatusCode sc = applyDownstreamCorr(inClusters, outClusters);
    if (sc.isFailure()) {
      return sc;
    }
    else{ 
      verbose() << "Running the downstream correction." << endmsg;
    }
  }

  // Apply benchmark correction
  if (m_benchmarkCorr) { 
    StatusCode sc = applyBenchmarkCorr(inClusters, outClusters);
    if (sc.isFailure()) {
      return sc;
    }
    else{
      verbose() << "Running the benchmark correction." << endmsg;
    }
  }

  if ((m_upstreamCorr && m_downstreamCorr && m_benchmarkCorr) || (m_upstreamCorr && m_benchmarkCorr) || (m_downstreamCorr && m_benchmarkCorr)){
    warning() << "Too many corrections in the house, apply upstream and downstream on ECal standalone and benchmark on combined ECal and HCal simulation." << endmsg;

  }

  return StatusCode::SUCCESS;
}


StatusCode CorrectCaloClusters::finalize() {

  return GaudiAlgorithm::finalize();
}


edm4hep::ClusterCollection* CorrectCaloClusters::initializeOutputClusters(
    const edm4hep::ClusterCollection* inClusters) {

  edm4hep::ClusterCollection* outClusters = m_outClusters.createAndPut();

  for (auto const& inCluster: *inClusters) {
    auto outCluster = inCluster.clone();
    verbose() << "Cluster position:" << endmsg;
    verbose() << "    x: " << outCluster.getPosition().x << endmsg;
    verbose() << "    y: " << outCluster.getPosition().y << endmsg;
    verbose() << "    z: " << outCluster.getPosition().z << endmsg;
    outClusters->push_back(outCluster);
  }

  return outClusters;
}


StatusCode CorrectCaloClusters::initializeCorrFunctions(std::vector<std::vector<TF1*>>& functions,
                                                        std::vector<std::vector<std::string>> formulas,
                                                        std::vector<std::vector<double>> parameters,
                                                        const std::string& funcNameStem) {

  for (size_t i = 0; i < formulas.size(); ++i) {
    auto& formulaVec = formulas[i];
    auto& paramVec = parameters[i];

    bool hasY = false;
    for (auto& formula: formulaVec) {
      if (formula.find("y") != std::string::npos) {
        hasY = true;
      }
    }

    std::vector<TF1*> funcVec;
    for (size_t j = 0; j < formulaVec.size(); ++j) {
      std::string funcName = "func_" + m_readoutNames[i] + "_";
      funcName += funcNameStem + "_";
      funcName += std::to_string(j);
      if (hasY) {
        TF2* func = new TF2(funcName.c_str(), formulaVec.at(j).c_str(), 0., 500., 0., 180.);
        funcVec.emplace_back(func);
      } else {
        TF1* func = new TF1(funcName.c_str(), formulaVec.at(j).c_str(), 0., 500.);
        funcVec.emplace_back(func);
      }
    }

    {
      size_t j = 0;
      for (auto& func: funcVec) {
        for (int k = 0; k < func->GetNpar(); ++k) {
          if (j >= paramVec.size()) {
            error() << "Correction parameter vector is not long enough!" << endmsg;
            return StatusCode::FAILURE;
          }
          func->SetParameter(k, paramVec.at(j));
          std::string parName(1, 'a' + (char) j);
          func->SetParName(k, parName.c_str());
          j += 1;
        }
      }
    }

    functions.emplace_back(funcVec);
  }
  return StatusCode::SUCCESS;
}


StatusCode CorrectCaloClusters::applyUpstreamCorr(const edm4hep::ClusterCollection* inClusters,
                                                  edm4hep::ClusterCollection* outClusters) {
  for (size_t i = 0; i < m_readoutNames.size(); ++i) {
    for (size_t j = 0; j < inClusters->size(); ++j) {
      double energyInFirstLayer = getEnergyInLayer(inClusters->at(j),
                                                   m_readoutNames[i],
                                                   m_systemIDs[i],
                                                   m_firstLayerIDs[i]);
      if (energyInFirstLayer < 0) {
        warning() << "Energy in first calorimeter layer negative, ignoring upstream energy correction!" << endmsg;
        continue;
      }

      const double clusterTheta = getClusterTheta(inClusters->at(j));
      verbose() << "Energy in first layer: " << energyInFirstLayer << endmsg;
      verbose() << "Cluster energy: " << inClusters->at(j).getEnergy() << endmsg;
      verbose() << "Cluster theta: " << clusterTheta << endmsg;

      verbose() << "Upstream correction:" << endmsg;
      double upstreamCorr = 0.;
      for (size_t k = 0; k < m_upstreamFunctions.at(i).size(); ++k) {
        auto func = m_upstreamFunctions.at(i).at(k);
        double corr = func->Eval(inClusters->at(j).getEnergy(), clusterTheta) * std::pow(energyInFirstLayer, k);
        verbose() << "    upsilon_" << k << " * E_firstLayer^" << k << ": " << corr << endmsg;
        upstreamCorr += corr;
      }

      verbose() << "    total: " << upstreamCorr << endmsg;
      outClusters->at(j).setEnergy(outClusters->at(j).getEnergy() + upstreamCorr);
      verbose() << "Corrected cluster energy: " << outClusters->at(j).getEnergy() << endmsg;
    }
  }

  return StatusCode::SUCCESS;
}


StatusCode CorrectCaloClusters::applyDownstreamCorr(const edm4hep::ClusterCollection* inClusters,
                                                    edm4hep::ClusterCollection* outClusters) {
  for (size_t i = 0; i < m_readoutNames.size(); ++i) {
    for (size_t j = 0; j < inClusters->size(); ++j) {
      double energyInLastLayer = getEnergyInLayer(inClusters->at(j),
                                                  m_readoutNames[i],
                                                  m_systemIDs[i],
                                                  m_lastLayerIDs[i]);
      if (energyInLastLayer < 0) {
        warning() << "Energy in last calorimeter layer negative, ignoring downstream energy correction!" << endmsg;
        continue;
      }

      const double clusterTheta = getClusterTheta(inClusters->at(j));
      verbose() << "Energy in last layer: " << energyInLastLayer << endmsg;
      verbose() << "Cluster energy: " << inClusters->at(j).getEnergy() << endmsg;
      verbose() << "Cluster theta: " << clusterTheta << endmsg;

      verbose() << "Downstream correction:" << endmsg;
      double downstreamCorr = 0.;
      for (size_t k = 0; k < m_downstreamFunctions.at(i).size(); ++k) {
        auto func = m_downstreamFunctions.at(i).at(k);
        double corr = func->Eval(inClusters->at(j).getEnergy(), clusterTheta) * std::pow(energyInLastLayer, k);
        verbose() << "    delta_" << k << " * E_lastLayer^" << k << ": " << corr << endmsg;
        downstreamCorr += corr;
      }

      verbose() << "    total: " << downstreamCorr << endmsg;
      outClusters->at(j).setEnergy(outClusters->at(j).getEnergy() + downstreamCorr);
      verbose() << "Corrected cluster energy: " << outClusters->at(j).getEnergy() << endmsg;
    }
  }

  return StatusCode::SUCCESS;
}

StatusCode CorrectCaloClusters::applyBenchmarkCorr(const edm4hep::ClusterCollection* inClusters,
                                                    edm4hep::ClusterCollection* outClusters) {

  const size_t numReadoutNames = m_readoutNames.size();

  int ecal_index = -1;
  int hcal_index = -1;

  // Identify ECal and HCal readout positions in the input parameters when running the calibration
  for (size_t i=0; i<numReadoutNames; ++i){
    if (m_systemIDs[i] == m_systemIDECal){
      ecal_index = i; 
    }
    else if (m_systemIDs[i] == m_systemIDHCal){
      hcal_index = i; 
    }
  }

  for (size_t j = 0; j < inClusters->size(); ++j) {
    double energyInLastLayerECal = getEnergyInLayer(inClusters->at(j),
                                                    m_readoutNames[ecal_index],
                                                    m_systemIDs[ecal_index],
                                                    m_lastLayerIDs[ecal_index]);
  
    double energyInFirstLayerECal = getEnergyInLayer(inClusters->at(j),
                                                     m_readoutNames[ecal_index],
                                                     m_systemIDs[ecal_index],
                                                     m_firstLayerIDs[ecal_index]);

    double energyInFirstLayerHCal = getEnergyInLayer(inClusters->at(j),
                                                     m_readoutNames[hcal_index],
                                                     m_systemIDs[hcal_index],
                                                     m_firstLayerIDs[hcal_index]);

    double totalEnergyInECal = getTotalEnergy(inClusters->at(j),
                                         m_readoutNames[ecal_index],
                                         m_systemIDs[ecal_index]);

    double totalEnergyInHCal = getTotalEnergy(inClusters->at(j),
                                         m_readoutNames[hcal_index],
                                         m_systemIDs[hcal_index]);
                                         

    // calculate approximate benchmark energy using non energy dependent benchmark parameters
    double approximateBenchmarkEnergy = 0.0; 
    approximateBenchmarkEnergy = m_benchmarkParamsApprox[0] * totalEnergyInECal +
                                 m_benchmarkParamsApprox[1] * totalEnergyInHCal +
                                 m_benchmarkParamsApprox[2] * sqrt(abs(energyInLastLayerECal * m_benchmarkParamsApprox[0] * energyInFirstLayerHCal * m_benchmarkParamsApprox[1])) +
                                 m_benchmarkParamsApprox[3] * pow(totalEnergyInECal * m_benchmarkParamsApprox[0], 2) +
                                 m_benchmarkParamsApprox[4] * energyInFirstLayerECal +
                                 m_benchmarkParamsApprox[5];

    // Calculate energy-dependent benchmark parameters p[0]-p[5]
    int benchmarkFormulaIndexHighEne = 0;
    int nParam = m_benchmarkFunctions.at(benchmarkFormulaIndexHighEne).size(); 
    std::vector<double> benchmarkParameters(nParam,-1.); 

    if (m_benchmarkEneSwitch>0.)
    {
      int benchmarkFormulaIndexLowEne = 1;
      // ensure smooth transition between low- and high-energy formulas (parameter l controls how fast the transition is)
      int l=2;
      double transition = 0.5 * (1 + tanh(l * (approximateBenchmarkEnergy - m_benchmarkEneSwitch) / 2));
      verbose() << "Using two formulas for benchmark calibration, the second formula provided will be used to correct energies below benchmarkEneSwitch threshold." << endmsg;
      // number of benchmark parameters must be the same for the two formulas 
      if (m_benchmarkFunctions.at(benchmarkFormulaIndexHighEne).size() != m_benchmarkFunctions.at(benchmarkFormulaIndexLowEne).size()){
        info() << "Size of the benchmarkFormulaIndexHighEne vector and benchmarkFormulaIndexLowEne vector does not match, exiting!" << endmsg;
        return StatusCode::FAILURE;
      }
      for (size_t k = 0; k < m_benchmarkFunctions.at(benchmarkFormulaIndexHighEne).size(); ++k) {
        auto func_low_ene = m_benchmarkFunctions.at(benchmarkFormulaIndexLowEne).at(k);
        auto func_high_ene = m_benchmarkFunctions.at(benchmarkFormulaIndexHighEne).at(k);
        benchmarkParameters[k] = (1 - transition) * func_low_ene->Eval(approximateBenchmarkEnergy) +
                                transition * func_high_ene->Eval(approximateBenchmarkEnergy);
      }
    }
    else{
      verbose() << "Using one formula for benchmark calibration." << endmsg;
      for (size_t k = 0; k < m_benchmarkFunctions.at(benchmarkFormulaIndexHighEne).size(); ++k) {
        auto func = m_benchmarkFunctions.at(benchmarkFormulaIndexHighEne).at(k);
        benchmarkParameters[k] = func->Eval(approximateBenchmarkEnergy);
      }
    }

    // Get final benchmark energy using the energy dependent benchmark parameters 
    double benchmarkEnergy = 0.0;
    
    benchmarkEnergy = benchmarkParameters[0] * totalEnergyInECal +
                      benchmarkParameters[1] * totalEnergyInHCal +
                      benchmarkParameters[2] * std::sqrt(std::fabs(energyInLastLayerECal * benchmarkParameters[0] * energyInFirstLayerHCal * benchmarkParameters[1])) +
                      benchmarkParameters[3] * std::pow(totalEnergyInECal * benchmarkParameters[0], 2) +
                      benchmarkParameters[4] * energyInFirstLayerECal +
                      benchmarkParameters[5];

    // Protection against negative energy (might be improved)
    if (benchmarkEnergy < 0.0) {
      outClusters->at(j).setEnergy(totalEnergyInECal * benchmarkParameters[0] + totalEnergyInHCal * benchmarkParameters[1]);
    } else {
      outClusters->at(j).setEnergy(benchmarkEnergy);
    }

    if (inClusters->at(j).getEnergy()>1.){ 
    debug() << "********************************************************************" << endmsg;
    debug() << "Cluster energy: " << inClusters->at(j).getEnergy() << endmsg;
    debug() << "totalEnergyInECal+HCal from hits: " << totalEnergyInECal+totalEnergyInHCal << endmsg;
    debug() << "********************************************************************" << endmsg;
    debug() << "Corrected cluster energy benchmark: " << outClusters->at(j).getEnergy() << endmsg;
    debug() << "********************************************************************" << endmsg;    
    } 
  }            
  return StatusCode::SUCCESS;
}

double CorrectCaloClusters::getEnergyInLayer(edm4hep::Cluster cluster,
                                             const std::string& readoutName,
                                             int systemID,
                                             int layerID) {
  dd4hep::DDSegmentation::BitFieldCoder* decoder = m_geoSvc->getDetector()->readout(readoutName).idSpec().decoder();

  double energy = 0;
  for (auto cell = cluster.hits_begin(); cell != cluster.hits_end(); ++cell) {
    dd4hep::DDSegmentation::CellID cellID = cell->getCellID();
    if (decoder->get(cellID, "system") != systemID) {
      continue;
    }
    if (decoder->get(cellID, "layer") != layerID) {
      continue;
    }
    energy += cell->getEnergy();
  }

  return energy;
}


double CorrectCaloClusters::getClusterTheta(edm4hep::Cluster cluster) {
  double rxy = std::sqrt(std::pow(cluster.getPosition().x, 2) + std::pow(cluster.getPosition().y, 2));
  double theta = ::fabs(std::atan2(rxy, cluster.getPosition().z));
  theta = 180 * theta / M_PI;

  return theta;
}


double CorrectCaloClusters::getTotalEnergy(edm4hep::Cluster cluster,
                                           const std::string& readoutName,
                                           int systemID) {
  dd4hep::DDSegmentation::BitFieldCoder* decoder = m_geoSvc->getDetector()->readout(readoutName).idSpec().decoder();

  double energy = 0;
  for (auto cell = cluster.hits_begin(); cell != cluster.hits_end(); ++cell) {
    dd4hep::DDSegmentation::CellID cellID = cell->getCellID();
    if (decoder->get(cellID, "system") != systemID) {
      continue;
    }
    energy += cell->getEnergy();
  }
  return energy;
}
