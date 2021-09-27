#include "CorrectCaloClusters.h"

// Gaudi
#include "GaudiKernel/ITHistSvc.h"

// Key4HEP
#include "k4Interface/IGeoSvc.h"

// FCC Detectors
#include "DetCommon/DetUtils.h"
#include "DetSegmentation/FCCSWGridPhiEta.h"

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
  StatusCode sc = GaudiAlgorithm::initialize();
  if (sc.isFailure()) {
    return sc;
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
  if (m_systemIDs.size() != m_upstreamFormulas.size()) {
    error() <<  "Sizes of systemIDs vector and upstreamFormulas vector does not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_upstreamParams.size()) {
    error() <<  "Sizes of systemIDs vector and upstreamParams vector does not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_downstreamFormulas.size()) {
    error() <<  "Sizes of systemIDs vector and downstreamFormulas vector does not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_downstreamParams.size()) {
    error() <<  "Sizes of systemIDs vector and downstreamParams vector does not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }

  // Prepare upstream and downstream correction functions
  initializeCorrFunctions(m_upstreamFunctions, m_upstreamFormulas, m_upstreamParams, "upstream");
  initializeCorrFunctions(m_downstreamFunctions, m_downstreamFormulas, m_downstreamParams, "downstream");

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
  {
    StatusCode sc = applyUpstreamCorr(inClusters, outClusters);
    if (sc.isFailure()) {
      return sc;
    }
  }

  // Apply downstream correction
  {
    StatusCode sc = applyDownstreamCorr(inClusters, outClusters);
    if (sc.isFailure()) {
      return sc;
    }
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
    edm4hep::Cluster outCluster = outClusters->create();

    outCluster.setPosition(inCluster.getPosition());
    verbose() << "Cluster position:" << endmsg;
    verbose() << "    x: " << outCluster.position().x << endmsg;
    verbose() << "    y: " << outCluster.position().y << endmsg;
    verbose() << "    z: " << outCluster.position().z << endmsg;
    outCluster.setEnergy(inCluster.getEnergy());
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


double CorrectCaloClusters::getEnergyInLayer(const edm4hep::Cluster& cluster,
                                             const std::string& readoutName,
                                             size_t systemID,
                                             size_t layerID) {
  dd4hep::DDSegmentation::BitFieldCoder* decoder = m_geoSvc->lcdd()->readout(readoutName).idSpec().decoder();

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


double CorrectCaloClusters::getClusterTheta(const edm4hep::Cluster& cluster) {
  double rxy = std::sqrt(std::pow(cluster.getPosition().x, 2) + std::pow(cluster.getPosition().y, 2));
  double theta = ::fabs(std::atan2(rxy, cluster.getPosition().z));
  theta = 180 * theta / M_PI;

  return theta;
}
