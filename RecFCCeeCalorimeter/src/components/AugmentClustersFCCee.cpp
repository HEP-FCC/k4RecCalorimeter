#include "AugmentClustersFCCee.h"

// k4FWCore
#include "k4Interface/IGeoSvc.h"

// edm4hep
#include "edm4hep/ClusterCollection.h"

// DD4hep
#include "DD4hep/Detector.h"

// ROOT
#include "TLorentzVector.h"
#include "TString.h"
#include "TVector3.h"
#include "TMath.h"

DECLARE_COMPONENT(AugmentClustersFCCee)

AugmentClustersFCCee::AugmentClustersFCCee(const std::string &name, ISvcLocator *svcLoc)
    : Gaudi::Algorithm(name, svcLoc),
      m_geoSvc("GeoSvc", "AugmentClustersFCCee")
{
  declareProperty("inClusters", m_inClusters, "Input clusters");
  declareProperty("outClusters", m_outClusters, "Output clusters");
}

StatusCode AugmentClustersFCCee::initialize()
{
  {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure())
      return sc;
  }

  // check if readouts exist
  for (size_t k = 0; k < m_readoutNames.size(); k++)
  {
    std::string readoutName = m_readoutNames[k];
    if (m_geoSvc->getDetector()->readouts().find(readoutName) == m_geoSvc->getDetector()->readouts().end())
    {
      error() << "Readout <<" << readoutName << ">> does not exist." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // Check if readout related variables have the same size
  if (m_systemIDs.size() != m_readoutNames.size())
  {
    error() << "Sizes of the systemIDs vector and readoutNames vector do not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_numLayers.size())
  {
    error() << "Sizes of systemIDs vector and numLayers vector do not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_layerFieldNames.size())
  {
    error() << "Sizes of systemIDs vector and layerFieldNames vector do not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_thetaFieldNames.size())
  {
    error() << "Sizes of systemIDs vector and thetaFieldNames vector do not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_moduleFieldNames.size())
  {
    error() << "Sizes of systemIDs vector and moduleFieldNames vector do not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_detectorNames.size())
  {
    error() << "Sizes of systemIDs vector and detectorNames vector do not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }

  // initialise the list of metadata for the clusters
  // append to the metadata of the input clusters
  std::vector<std::string> showerShapeDecorations = m_inShapeParameterHandle.get();
  for (size_t k = 0; k < m_detectorNames.size(); k++)
  {
    const char *detector = m_detectorNames[k].c_str();
    for (unsigned layer = 0; layer < m_numLayers[k]; layer++)
    {
      showerShapeDecorations.push_back(Form("energy_fraction_%s_layer_%d", detector, layer));
      showerShapeDecorations.push_back(Form("theta_%s_layer_%d", detector, layer));
      showerShapeDecorations.push_back(Form("phi_%s_layer_%d", detector, layer));
      showerShapeDecorations.push_back(Form("width_theta_%s_layer_%d", detector, layer));
      showerShapeDecorations.push_back(Form("width_module_%s_layer_%d", detector, layer));
      showerShapeDecorations.push_back(Form("Ratio_E_max_2ndmax_%s_layer_%d", detector, layer));
      showerShapeDecorations.push_back(Form("Delta_E_2ndmax_min_%s_layer_%d", detector, layer));
      showerShapeDecorations.push_back(Form("width_theta_3Bin_%s_layer_%d", detector, layer));
    }
  }
  m_showerShapeHandle.put(showerShapeDecorations);

  return StatusCode::SUCCESS;
}

std::pair<std::vector<int>, std::vector<double>> MergeSumAndSort(std::vector<int>& A, std::vector<double>& B) {
  std::unordered_map<int, double> elementSum;
  // traverse vec A, merge the same elements (theta ID) in vec A and sum up the corresponding elements (E_cell) in vec B
  for (size_t i = 0; i < A.size(); i ++) {
    elementSum[A[i]] += B[i];
  }
  std::vector<int> A_new;
  std::vector<double> B_new;
  // re-build vec A and vec B from elementSum
  for (const auto& entry : elementSum) {
    A_new.push_back(entry.first);
    B_new.push_back(entry.second);
  }
  std::vector<int> vec_1(A_new.size());
  std::vector<double> vec_2(B_new.size());
  std::vector<size_t> indices(A_new.size());
  for (size_t i = 0; i < indices.size(); ++ i) {
    indices[i] = i;
  }
  // sort the theta vec, update the E vec simultaneously
  std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
    return A_new[a] < A_new[b];
  });
  for (size_t i = 0; i < indices.size(); ++ i) {
    vec_1[i] = A_new[indices[i]];
    vec_2[i] = B_new[indices[i]];
  }
  return std::make_pair(vec_1, vec_2);
}

StatusCode AugmentClustersFCCee::finalize()
{
  return Gaudi::Algorithm::finalize();
}

StatusCode AugmentClustersFCCee::execute([[maybe_unused]] const EventContext &evtCtx) const
{
  // get the input collection with clusters
  const edm4hep::ClusterCollection *inClusters = m_inClusters.get();

  // create the new output collection
  edm4hep::ClusterCollection *outClusters = m_outClusters.createAndPut();

  // total number of layers
  size_t numLayersTotal = 0;
  for (size_t i = 0; i < m_numLayers.size(); ++i)
  {
    numLayersTotal += m_numLayers[i];
  }
  
  // loop over the clusters, clone them, and calculate the shape parameters to store with them
  for (const auto &cluster : *inClusters)
  {
    // clone original cluster
    auto newCluster = cluster.clone();
    outClusters->push_back(newCluster);

    // calculate the energy deposited in each layer
    // also, find out if cluster is around -pi..pi transition
    double E(0.0);
    std::vector<double> sumEnLayer;
    sumEnLayer.assign(numLayersTotal, 0.);
    double phiMin = 9999.;
    double phiMax = -9999.;
    int module_id_Min = 1536;
    int module_id_Max = -1;

    // loop over each system/readout
    int startPositionToFill = 0;
    for (size_t k = 0; k < m_readoutNames.size(); k++)
    {
      if (k > 0)
      {
        startPositionToFill += m_numLayers[k - 1];  // only 2 sub detectors ?
      }
      int systemID = m_systemIDs[k];
      std::string layerField = m_layerFieldNames[k];
      std::string moduleField = m_moduleFieldNames[k];
      std::string readoutName = m_readoutNames[k];
      dd4hep::DDSegmentation::BitFieldCoder *decoder = m_geoSvc->getDetector()->readout(readoutName).idSpec().decoder();

      // 1st loop over the cells
      for (auto cell = newCluster.hits_begin(); cell != newCluster.hits_end(); cell++)
      {
        dd4hep::DDSegmentation::CellID cID = cell->getCellID();
        int sysId = decoder->get(cID, "system");
        if (sysId != systemID)
          continue;
        uint layer = decoder->get(cID, layerField);
        int module_id = decoder->get(cID, moduleField);
        double eCell = cell->getEnergy();  // to MeV in calculation ??
        sumEnLayer[layer+startPositionToFill] += eCell;
        E += eCell;
        TVector3 v = TVector3(cell->getPosition().x, cell->getPosition().y, cell->getPosition().z);
        double phi = v.Phi();
        if (phi < phiMin)
          phiMin = phi;
        if (phi > phiMax)
          phiMax = phi;

        if (module_id > module_id_Max)    module_id_Max = module_id;
        if (module_id < module_id_Min)    module_id_Min = module_id;

      }  // end of loop over cells
    }  // end of loop over system / readout
    // in the loop above we get E_layer, E_cluster, max and min of cluster phi
  
    // any number close to two pi should do, because if a cluster contains
    // the -pi<->pi transition, phiMin should be close to -pi and phiMax close to pi
    bool isClusterPhiNearPi = false;
    if (phiMax - phiMin > 6.)
      isClusterPhiNearPi = true;
    debug() << "phiMin, phiMax : " << phiMin << " " << phiMax << endmsg;
    debug() << "Cluster is near phi=pi : " << isClusterPhiNearPi << endmsg;

    bool isResetModuleID = false;
    if (module_id_Max - module_id_Min > 1500)    isResetModuleID = true;
    //debug() << "phiMin, phiMax : " << phiMin << " " << phiMax << endmsg;
    //debug() << "Cluster is near phi=pi : " << isClusterPhiNearPi << endmsg;

    // calculate the theta positions with log(E) weighting in each layer
    // for phi use standard E weighting
    std::vector<double> sumThetaLayer;
    std::vector<double> sumPhiLayer;
    std::vector<double> sumWeightLayer;
    sumThetaLayer.assign(numLayersTotal, 0);
    sumPhiLayer.assign(numLayersTotal, 0);
    sumWeightLayer.assign(numLayersTotal, 0);

    // for theta/module width calculation
    std::vector<double> theta2_E_layer(numLayersTotal, 0.);
    std::vector<double> theta_E_layer(numLayersTotal, 0.);
    std::vector<double> module2_E_layer(numLayersTotal, 0.);
    std::vector<double> module_E_layer(numLayersTotal, 0.);

    // E, theta and module ID of cells in each layer
    std::vector<std::vector<double>> vec_E_cell_layer(numLayersTotal, std::vector<double>());
    std::vector<std::vector<int>> vec_theta_cell_layer(numLayersTotal, std::vector<int>());
    std::vector<std::vector<int>> vec_module_cell_layer(numLayersTotal, std::vector<int>());

    //std::vector<double> E_frac_Layer(numLayersTotal, 0.);
    std::vector<double> width_theta(numLayersTotal, 0.);
    std::vector<double> width_module(numLayersTotal, 0.);
    std::vector<double> width_theta_3Bin(numLayersTotal, 0.);
    std::vector<double> Ratio_E_max_2ndmax(numLayersTotal, 0.);
    std::vector<double> Delta_E_2ndmax_min(numLayersTotal, 0.);

    // loop over each system/readout
    startPositionToFill = 0;
    // rather than a loop over the systems, could first determine systemID
    // from cellID, match it against the systemIDs and skip cell if match not found
    // can do it in a separate PR
    for (size_t k = 0; k < m_readoutNames.size(); k++)
    {
      if (k > 0)
      {
        startPositionToFill += m_numLayers[k - 1];
      }
      int systemID = m_systemIDs[k];
      std::string layerField = m_layerFieldNames[k];
      std::string thetaField = m_thetaFieldNames[k];
      std::string moduleField = m_moduleFieldNames[k];
      std::string readoutName = m_readoutNames[k];
      dd4hep::DDSegmentation::BitFieldCoder *decoder = m_geoSvc->getDetector()->readout(readoutName).idSpec().decoder();

      // 2nd loop over the cells
      for (auto cell = newCluster.hits_begin(); cell != newCluster.hits_end(); cell++)
      {
        dd4hep::DDSegmentation::CellID cID = cell->getCellID();
        int sysId = decoder->get(cID, "system");
        if (sysId != systemID)
          continue;
        uint layer = decoder->get(cID, layerField);
        int theta_id = decoder->get(cID, thetaField);
        int module_id = decoder->get(cID, moduleField);

        double eCell = cell->getEnergy();
        double weightLog = std::max(0., m_thetaRecalcLayerWeights[k][layer] + log(eCell / sumEnLayer[layer]));  // TODO layer+startPositionToFill ?
        TVector3 v = TVector3(cell->getPosition().x, cell->getPosition().y, cell->getPosition().z);
        double theta = v.Theta();
        double phi = v.Phi();

        // for clusters that are around the -pi<->pi transition, we want to avoid averaging
        // over phi values that might differ by 2pi. in that case, for cells with negative
        // phi we add two pi, so that we average phi values all close to pi
        if (isClusterPhiNearPi && phi < 0.) {
          phi += TMath::TwoPi();
        }
        if (isResetModuleID && module_id > 1536/2) {
          module_id -= 1536;  // transition of module ID
        }

        if (m_thetaRecalcLayerWeights[k][layer]<0)
          sumThetaLayer[layer+startPositionToFill] += (eCell * theta);
        else
          sumThetaLayer[layer+startPositionToFill] += (weightLog * theta);
        sumWeightLayer[layer+startPositionToFill] += weightLog;
        sumPhiLayer[layer+startPositionToFill] += (eCell * phi);

        vec_E_cell_layer[layer+startPositionToFill].push_back(eCell);
        vec_theta_cell_layer[layer+startPositionToFill].push_back(theta_id);
        vec_module_cell_layer[layer+startPositionToFill].push_back(module_id);

        // sum them for width calculation
        theta2_E_layer[layer+startPositionToFill] += theta_id * theta_id * eCell;
        theta_E_layer[layer+startPositionToFill] += theta_id * eCell;
        module2_E_layer[layer+startPositionToFill] += module_id * module_id * eCell;
        module_E_layer[layer+startPositionToFill] += module_id * eCell;
      }  // end of loop over cells
    }  // end of loop over each system / readout

    std::vector<std::pair<std::vector<int>, std::vector<double>>> theta_E_pair;
    // local maxima (could be more than one) and the corresponding theta
    std::vector<std::vector<double>> local_E_Max(numLayersTotal, std::vector<double>());
    std::vector<std::vector<int>> local_E_Max_theta(numLayersTotal, std::vector<int>());

    std::vector<double> E_cell_Max(numLayersTotal, 0.);
    std::vector<double> E_cell_secMax(numLayersTotal, 0.);
    std::vector<double> E_cell_Max_theta(numLayersTotal, 0.);
    std::vector<double> E_cell_secMax_theta(numLayersTotal, 0.);
    std::vector<double> E_cell_Min(numLayersTotal, std::numeric_limits<double>::max());

    startPositionToFill = 0;
    for (size_t k = 0; k < m_readoutNames.size(); k++) {
      if (k > 0)    startPositionToFill += m_numLayers[k - 1];
      // loop over layers
      for (unsigned layer = 0; layer < m_numLayers[k]; layer++) {
        // in case there's no cell in this layer (sometimes in layer 0)
        if (vec_E_cell_layer[layer+startPositionToFill].empty()) {
          vec_E_cell_layer[layer+startPositionToFill].push_back(0);
          vec_theta_cell_layer[layer+startPositionToFill].push_back(0);
          vec_module_cell_layer[layer+startPositionToFill].push_back(0);
        }
        auto result = MergeSumAndSort(vec_theta_cell_layer[layer+startPositionToFill], vec_E_cell_layer[layer+startPositionToFill]);
        theta_E_pair.push_back(result);

        // loop over theta IDs
        for (size_t i = 0; i < theta_E_pair[layer+startPositionToFill].second.size(); i ++) {
          // find the local E maxima
          // discard the first and last cell in theta
          // local max: E_cell is greater than the 2 neighbors
          if (i != 0 && i != (theta_E_pair[layer+startPositionToFill].second.size()-1) && 
              theta_E_pair[layer+startPositionToFill].second[i] > theta_E_pair[layer+startPositionToFill].second[i-1] && 
              theta_E_pair[layer+startPositionToFill].second[i] > theta_E_pair[layer+startPositionToFill].second[i+1]) {
            local_E_Max[layer+startPositionToFill].push_back(theta_E_pair[layer+startPositionToFill].second[i]);
            local_E_Max_theta[layer+startPositionToFill].push_back(theta_E_pair[layer+startPositionToFill].first[i]);
          }
        }  // end of loop over theta IDs

        if (local_E_Max[layer+startPositionToFill].empty()) {
          E_cell_Max[layer+startPositionToFill] = 0.;
          E_cell_secMax[layer+startPositionToFill] = 0.;
          E_cell_Max_theta[layer+startPositionToFill] = 0.;
          E_cell_secMax_theta[layer+startPositionToFill] = 0.;
          E_cell_Min[layer+startPositionToFill] = 0.;
        } else if (local_E_Max[layer+startPositionToFill].size() < 2) {
          E_cell_Max[layer+startPositionToFill] = local_E_Max[layer+startPositionToFill][0];
          E_cell_secMax[layer+startPositionToFill] = 0.;
          E_cell_Max_theta[layer+startPositionToFill] = local_E_Max_theta[layer+startPositionToFill][0];
          E_cell_secMax_theta[layer+startPositionToFill] = local_E_Max_theta[layer+startPositionToFill][0];
          E_cell_Min[layer+startPositionToFill] = 0.;
        } else {
          std::vector<double> sortedVec = local_E_Max[layer+startPositionToFill];
          // move the top 2 max to the beginning
          std::partial_sort(sortedVec.begin(), sortedVec.begin()+2, sortedVec.end(), std::greater<double>());
          E_cell_Max[layer+startPositionToFill] = sortedVec[0];
          E_cell_secMax[layer+startPositionToFill] = sortedVec[1];
          // get the corresponding theta IDs
          auto it_Max   = std::find(local_E_Max[layer+startPositionToFill].begin(), local_E_Max[layer+startPositionToFill].end(), sortedVec[0]);
          int index_Max = std::distance(local_E_Max[layer+startPositionToFill].begin(), it_Max);
          auto it_secMax   = std::find(local_E_Max[layer+startPositionToFill].begin(), local_E_Max[layer+startPositionToFill].end(), sortedVec[1]);
          int index_secMax = std::distance(local_E_Max[layer+startPositionToFill].begin(), it_secMax);
          E_cell_Max_theta[layer+startPositionToFill]    = local_E_Max_theta[layer+startPositionToFill][index_Max];
          E_cell_secMax_theta[layer+startPositionToFill] = local_E_Max_theta[layer+startPositionToFill][index_secMax];
          // find the E_min inside the theta range of E_cell_Max and E_cell_secMax
          for (size_t i = 0; i < theta_E_pair[layer+startPositionToFill].first.size(); i ++ ) {
            if (theta_E_pair[layer+startPositionToFill].first[i] > std::min(E_cell_Max_theta[layer+startPositionToFill], E_cell_secMax_theta[layer+startPositionToFill])
                 && theta_E_pair[layer+startPositionToFill].first[i] < std::max(E_cell_Max_theta[layer+startPositionToFill], E_cell_secMax_theta[layer+startPositionToFill])) {
              if (theta_E_pair[layer+startPositionToFill].second[i] < E_cell_Min[layer+startPositionToFill]) {
                E_cell_Min[layer+startPositionToFill] = theta_E_pair[layer+startPositionToFill].second[i];
              }
            }
          }
        }
        if (E_cell_Min[layer+startPositionToFill] > 1e12)  E_cell_Min[layer+startPositionToFill] = 0.;  // check E_cell_Min
      }  // end of loop over layers
    }

    // save energy and theta/phi positions per layer in shape parameters
    startPositionToFill = 0;
    for (size_t k = 0; k < m_readoutNames.size(); k++)
    {
      if (k > 0)
      {
        startPositionToFill += m_numLayers[k - 1];
      }
      // loop over layers
      for (unsigned layer = 0; layer < m_numLayers[k]; layer++)
      {
        // theta
        if (m_thetaRecalcLayerWeights[k][layer]<0)
        {
          if (sumEnLayer[layer+startPositionToFill] != 0.0)
          {
            sumThetaLayer[layer+startPositionToFill] /= sumEnLayer[layer+startPositionToFill];
          }
        }
        else
        {
          if (sumWeightLayer[layer+startPositionToFill] != 0.0)
          {
            sumThetaLayer[layer+startPositionToFill] /= sumWeightLayer[layer+startPositionToFill];
          }
        }

        // phi
        if (sumEnLayer[layer+startPositionToFill] != 0.0)
        {
          sumPhiLayer[layer+startPositionToFill] /= sumEnLayer[layer+startPositionToFill];
        }
        // make sure phi is in range -pi..pi
        if (sumPhiLayer[layer+startPositionToFill] > TMath::Pi())
          sumPhiLayer[layer+startPositionToFill] -= TMath::TwoPi();

        newCluster.addToShapeParameters(sumEnLayer[layer+startPositionToFill] / E);  // E fraction of layer
        newCluster.addToShapeParameters(sumThetaLayer[layer+startPositionToFill]);
        newCluster.addToShapeParameters(sumPhiLayer[layer+startPositionToFill]);

        //E_frac_Layer[layer+startPositionToFill] = sumEnLayer[layer+startPositionToFill] / E;
        //if (E == 0)    E_frac_Layer[layer+startPositionToFill] = 0.;
        //newCluster.addToShapeParameters(E_frac_Layer[layer+startPositionToFill]);

        double w_theta = sqrt(fabs(theta2_E_layer[layer+startPositionToFill] / sumEnLayer[layer+startPositionToFill] - std::pow(theta_E_layer[+startPositionToFill] / sumEnLayer[layer+startPositionToFill], 2)));
        if (std::isnan(w_theta))    w_theta = 0.;
        if (w_theta > 40)    w_theta = 40;
        width_theta[layer+startPositionToFill] = w_theta;

        double w_module = sqrt(fabs(module2_E_layer[layer+startPositionToFill] / sumEnLayer[layer+startPositionToFill] - std::pow(module_E_layer[+startPositionToFill] / sumEnLayer[layer+startPositionToFill], 2)));
        if (std::isnan(w_module))    w_module = 0.;
        if (w_module > 40)    w_module = 40;
        width_module[layer+startPositionToFill] = w_module;

        double Ratio_E = (E_cell_Max[layer+startPositionToFill] - E_cell_secMax[layer+startPositionToFill]) /
                         (E_cell_Max[layer+startPositionToFill] + E_cell_secMax[layer+startPositionToFill]);
        if (E_cell_Max[layer+startPositionToFill] + E_cell_secMax[layer+startPositionToFill] == 0)    Ratio_E = 1.;
        Ratio_E_max_2ndmax[layer+startPositionToFill] = Ratio_E;
        Delta_E_2ndmax_min[layer+startPositionToFill] = E_cell_secMax[layer+startPositionToFill] - E_cell_Min[layer+startPositionToFill];

        if (local_E_Max[layer+startPositionToFill].size() > 0) {
          double E_m1 = 0.;
          double E_p1 = 0.;
          int theta_m1 = 0;
          int theta_p1 = 0;
          double theta2_E_3Bin = 0.;
          double theta_E_3Bin = 0.;
          double sum_E_3Bin = 0.;

          auto it_1 = std::find(theta_E_pair[layer+startPositionToFill].second.begin(), theta_E_pair[layer+startPositionToFill].second.end(), E_cell_Max[layer+startPositionToFill]);
          int ind_1 = std::distance(theta_E_pair[layer+startPositionToFill].second.begin(), it_1);
          E_m1 = theta_E_pair[layer+startPositionToFill].second[ind_1-1];
          E_p1 = theta_E_pair[layer+startPositionToFill].second[ind_1+1];
          theta_m1 = theta_E_pair[layer+startPositionToFill].first[ind_1-1];
          theta_p1 = theta_E_pair[layer+startPositionToFill].first[ind_1+1];
          sum_E_3Bin = E_m1 + E_cell_Max[layer+startPositionToFill] + E_p1;
          theta2_E_3Bin = theta_m1 * theta_m1 * E_m1 + E_cell_Max_theta[layer+startPositionToFill] * E_cell_Max_theta[layer+startPositionToFill] * E_cell_Max[layer+startPositionToFill] + theta_p1 * theta_p1 * E_p1;
          theta_E_3Bin = theta_m1 * E_m1 + E_cell_Max_theta[layer+startPositionToFill] * E_cell_Max[layer+startPositionToFill] + theta_p1 * E_p1;

          double _w_theta_3Bin = sqrt(fabs(theta2_E_3Bin / sum_E_3Bin - std::pow(theta_E_3Bin / sum_E_3Bin, 2)));
          if (std::isnan(_w_theta_3Bin))    _w_theta_3Bin = 0.;
          width_theta_3Bin[layer+startPositionToFill] = _w_theta_3Bin;
        } else {
          width_theta_3Bin[layer+startPositionToFill] = 0.;
        }

        newCluster.addToShapeParameters(width_theta[layer+startPositionToFill]);
        newCluster.addToShapeParameters(width_module[layer+startPositionToFill]);
        newCluster.addToShapeParameters(Ratio_E_max_2ndmax[layer+startPositionToFill]);
        newCluster.addToShapeParameters(Delta_E_2ndmax_min[layer+startPositionToFill]);
        newCluster.addToShapeParameters(width_theta_3Bin[layer+startPositionToFill]);

      }  // end of loop over layers
    }  // end of loop over system/readout
  }  // end of loop over clusters

  return StatusCode::SUCCESS;
}
