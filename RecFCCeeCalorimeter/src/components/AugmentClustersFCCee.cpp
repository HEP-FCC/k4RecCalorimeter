#include "AugmentClustersFCCee.h"

// k4FWCore
#include "k4Interface/IGeoSvc.h"

// edm4hep
#include "edm4hep/ClusterCollection.h"

// DD4hep
#include "DD4hep/Detector.h"

// k4geo
#include "detectorSegmentations/FCCSWGridPhiTheta_k4geo.h"
#include "detectorSegmentations/FCCSWGridModuleThetaMerged_k4geo.h"

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

  // retrieve systemID of EMB from "DetID_ECAL_Barrel"
  systemID_EMB = m_geoSvc->getDetector()->constantAsDouble("DetID_ECAL_Barrel");

  // retrieve some information from the segmentation for later use
  // - number of modules/phi bins => needed to account for module/phi periodicity
  // - number of merged modules and theta cells per layer => needed to check if two cells are neighbours
  for (size_t k = 0; k < m_readoutNames.size(); k++)
  {
    std::string readoutName = m_readoutNames[k];
    dd4hep::DDSegmentation::Segmentation *aSegmentation = m_geoSvc->getDetector()->readout(readoutName).segmentation().segmentation();
    std::string segmentationType = aSegmentation->type();
    if (segmentationType == "FCCSWGridModuleThetaMerged_k4geo")
    {
      dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo *moduleThetaSegmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo *>(aSegmentation);
      nModules.push_back(moduleThetaSegmentation->nModules());
      for (size_t iLayer = 0; iLayer < m_numLayers[k]; ++iLayer)
      {
        nMergedThetaCells.push_back(moduleThetaSegmentation->mergedThetaCells(iLayer));
        nMergedModules.push_back(moduleThetaSegmentation->mergedModules(iLayer));
      }
    }
    else if (segmentationType == "FCCSWGridPhiTheta_k4geo")
    {
      dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo *phiThetaSegmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo *>(aSegmentation);
      nModules.push_back(phiThetaSegmentation->phiBins());
      for (size_t iLayer = 0; iLayer < m_numLayers[k]; ++iLayer)
      {
        nMergedThetaCells.push_back(1);
        nMergedModules.push_back(1);
      }
    }
    else
    {
      error() << "Segmentation type not handled, aborting..." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // initialise the list of metadata for the clusters
  // append to the metadata of the input clusters (if any)
  std::vector<std::string> showerShapeDecorations = m_inShapeParameterHandle.get({});
  for (size_t k = 0; k < m_detectorNames.size(); k++)
  {
    const char *detector = m_detectorNames[k].c_str();
    for (unsigned layer = 0; layer < m_numLayers[k]; layer++)
    {
      showerShapeDecorations.push_back(Form("energy_fraction_%s_layer_%d", detector, layer));
      showerShapeDecorations.push_back(Form("theta_%s_layer_%d", detector, layer));
      showerShapeDecorations.push_back(Form("phi_%s_layer_%d", detector, layer));
      showerShapeDecorations.push_back(Form("maxcell_E_%s_layer_%d", detector, layer));
      // pi0/photon shape var only calculated in EMB
      if (m_do_photon_shapeVar && m_systemIDs[k] == systemID_EMB) {
        showerShapeDecorations.push_back(Form("width_theta_%s_layer_%d", detector, layer));
        showerShapeDecorations.push_back(Form("width_module_%s_layer_%d", detector, layer));
        showerShapeDecorations.push_back(Form("Ratio_E_max_2ndmax_%s_layer_%d", detector, layer));
        showerShapeDecorations.push_back(Form("Delta_E_2ndmax_min_%s_layer_%d", detector, layer));
        showerShapeDecorations.push_back(Form("Ratio_E_max_2ndmax_vs_phi_%s_layer_%d", detector, layer));
        showerShapeDecorations.push_back(Form("Delta_E_2ndmax_min_vs_phi_%s_layer_%d", detector, layer));
        showerShapeDecorations.push_back(Form("width_theta_3Bin_%s_layer_%d", detector, layer));
        showerShapeDecorations.push_back(Form("width_theta_5Bin_%s_layer_%d", detector, layer));
        showerShapeDecorations.push_back(Form("width_theta_7Bin_%s_layer_%d", detector, layer));
        showerShapeDecorations.push_back(Form("width_theta_9Bin_%s_layer_%d", detector, layer));
        showerShapeDecorations.push_back(Form("E_fr_side_pm2_%s_layer_%d", detector, layer));
        showerShapeDecorations.push_back(Form("E_fr_side_pm3_%s_layer_%d", detector, layer));
        showerShapeDecorations.push_back(Form("E_fr_side_pm4_%s_layer_%d", detector, layer));
      }
    }
  }
  showerShapeDecorations.push_back("mass"); // cluster invariant mass assuming massless constituents
  showerShapeDecorations.push_back("ncells"); // number of cells in cluster with E>0

  m_showerShapeHandle.put(showerShapeDecorations);

  return StatusCode::SUCCESS;
}

// for all cells in a certain layer: vec A is theta_id, vec B is E_cell
// make a 1D projection to theta: sum up the E_cell over modules with the same theta_id
// then sort the theta vec (for finding theta neighbors), update the E vec simultaneously
// return a pair of vectors: theta_id (in ascending order) and E (summed over modules)
std::pair<std::vector<int>, std::vector<double>> MergeSumAndSort(const std::vector<int>& A, const std::vector<double>& B) {
  std::unordered_map<int, double> elementSum;
  // merge the same elements in vec A and sum up the corresponding elements in vec B
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

    // loop over all cells to:
    // - calculate the cluster invariant mass
    // - calculate the energy deposited in each layer
    // - find out the energy of the cells with largest energy in each layer
    // - find out if cluster is around -pi..pi transition and/or max module .. 0 transition
    double E(0.0);
    TLorentzVector p4cl(0.0, 0.0, 0.0, 0.0);
    unsigned int nCells(0);
    std::vector<double> sumEnLayer(numLayersTotal, 0.);
    std::vector<double> maxCellEnergyInLayer(numLayersTotal, 0.);
    double phiMin = 9999.;
    double phiMax = -9999.;
    int module_id_Min = 9999;
    int module_id_Max = -9999;

    // loop over each system/readout
    int startPositionToFill = 0;
    for (size_t k = 0; k < m_readoutNames.size(); k++)
    {
      if (k > 0)
        startPositionToFill += m_numLayers[k - 1];
      int systemID = m_systemIDs[k];
      std::string layerField = m_layerFieldNames[k];
      std::string moduleField = m_moduleFieldNames[k];
      std::string readoutName = m_readoutNames[k];
      dd4hep::DDSegmentation::BitFieldCoder *decoder = m_geoSvc->getDetector()->readout(readoutName).idSpec().decoder();

      // loop over the cells to get E_layer, E_cluster, max and min phi of cluster cells
      for (auto cell = newCluster.hits_begin(); cell != newCluster.hits_end(); cell++)
      {
        dd4hep::DDSegmentation::CellID cID = cell->getCellID();

        // skip cells with wrong system ID
        int sysId = decoder->get(cID, "system");
        if (sysId != systemID)
          continue;

        // retrieve layer and module ID
        uint layer = decoder->get(cID, layerField);
        int module_id = decoder->get(cID, moduleField);

        // retrieve cell energy, update sum of cell energies per layer, total energy, and max cell energy
        double eCell = cell->getEnergy();
        sumEnLayer[layer+startPositionToFill] += eCell;
        E += eCell;
        if (maxCellEnergyInLayer[layer+startPositionToFill]<eCell)
          maxCellEnergyInLayer[layer+startPositionToFill] = eCell;

        // compute phi and module of cell to determine min/max phi and module ID of cluster
        TVector3 v = TVector3(cell->getPosition().x, cell->getPosition().y, cell->getPosition().z);
        double phi = v.Phi();
        if (phi < phiMin) phiMin = phi;
        if (phi > phiMax) phiMax = phi;
        if (module_id > module_id_Max) module_id_Max = module_id;
        if (module_id < module_id_Min) module_id_Min = module_id;

        // add cell 4-momentum to cluster 4-momentum
        TVector3 pCell = v * (eCell / v.Mag());
        TLorentzVector p4cell(pCell.X(), pCell.Y(), pCell.Z(), eCell);
        p4cl += p4cell;
        nCells++;
      }  // end of loop over cells
    }  // end of loop over system / readout
    //std::cout << "Ecl = " << E << std::endl;

    // any number close to two pi should do, because if a cluster contains
    // the -pi<->pi transition, phiMin should be close to -pi and phiMax close to pi
    bool isClusterPhiNearPi = false;
    if (phiMax - phiMin > 6.)    isClusterPhiNearPi = true;

    debug() << "phiMin, phiMax : " << phiMin << " " << phiMax << endmsg;
    debug() << "Cluster is near phi=pi : " << isClusterPhiNearPi << endmsg;

    bool isResetModuleID = false;
    // near the 1535..0 transition, reset module ID
    if (module_id_Max - module_id_Min > nModules[0] * .9)    isResetModuleID = true;

    // calculate the theta positions with log(E) weighting in each layer
    // for phi use standard E weighting
    std::vector<double> sumThetaLayer;
    std::vector<double> sumPhiLayer;
    std::vector<double> sumWeightLayer;
    sumThetaLayer.assign(numLayersTotal, 0);
    sumPhiLayer.assign(numLayersTotal, 0);
    sumWeightLayer.assign(numLayersTotal, 0);

    // vectors for photon/pi0 discrimination
    // vectors to calculate theta/module width vs layer - using all cells
    std::vector<double> theta2_E_layer(numLayersTotal, 0.);
    std::vector<double> theta_E_layer(numLayersTotal, 0.);
    std::vector<double> module2_E_layer(numLayersTotal, 0.);
    std::vector<double> module_E_layer(numLayersTotal, 0.);

    // E, theta and module ID of cells in each layer
    std::vector<std::vector<double>> vec_E_cell_layer(numLayersTotal, std::vector<double>());
    std::vector<std::vector<int>> vec_theta_cell_layer(numLayersTotal, std::vector<int>());
    std::vector<std::vector<int>> vec_module_cell_layer(numLayersTotal, std::vector<int>());

    // E and theta of 1st and 2nd local max in 1D theta profile, and E of minimum between them
    std::vector<double> E_cell_Max(numLayersTotal, 0.);
    std::vector<double> E_cell_secMax(numLayersTotal, 0.);
    std::vector<int> E_cell_Max_theta(numLayersTotal, 0.);
    std::vector<int> E_cell_secMax_theta(numLayersTotal, 0.);
    std::vector<double> E_cell_Min(numLayersTotal, std::numeric_limits<double>::max());

    // E and theta of 1st and 2nd local max in 1D module profile, and E of minimum between them
    std::vector<double> E_cell_vs_phi_Max(numLayersTotal, 0.);
    std::vector<double> E_cell_vs_phi_secMax(numLayersTotal, 0.);
    std::vector<int> E_cell_vs_phi_Max_module(numLayersTotal, 0.);
    std::vector<int> E_cell_vs_phi_secMax_module(numLayersTotal, 0.);
    std::vector<double> E_cell_vs_phi_Min(numLayersTotal, std::numeric_limits<double>::max());

    // theta/module width using all cells
    std::vector<double> width_theta(numLayersTotal, 0.);
    std::vector<double> width_module(numLayersTotal, 0.);
    // theta width using only cells within deltaThetaBin = +-1, +-2, +-3, +-4
    std::vector<double> width_theta_3Bin(numLayersTotal, 0.);
    std::vector<double> width_theta_5Bin(numLayersTotal, 0.);
    std::vector<double> width_theta_7Bin(numLayersTotal, 0.);
    std::vector<double> width_theta_9Bin(numLayersTotal, 0.);

    // E_fr_side_N, (E around maxE +- N bins) / (E around maxE +- 1 bins) - 1.
    std::vector<double> E_fr_side_pm2(numLayersTotal, 0.);
    std::vector<double> E_fr_side_pm3(numLayersTotal, 0.);
    std::vector<double> E_fr_side_pm4(numLayersTotal, 0.);

    // (Emax - E2ndmax)/(Emax + E2ndmax) where 2nd max must be a local maximum
    std::vector<double> Ratio_E_max_2ndmax(numLayersTotal, 0.);
    // (E2ndmax - Emin) where Emin is the energy with minimum energy in the theta range defined by 1st and 2nd (local) max
    std::vector<double> Delta_E_2ndmax_min(numLayersTotal, 0.);
    // (Emax - E2ndmax)/(Emax + E2ndmax) where 2nd max must be a local maximum - in phi profile
    std::vector<double> Ratio_E_max_2ndmax_vs_phi(numLayersTotal, 0.);
    // (E2ndmax - Emin) where Emin is the energy with minimum energy in the module range defined by 1st and 2nd (local) max
    std::vector<double> Delta_E_2ndmax_min_vs_phi(numLayersTotal, 0.);

    // loop over each system/readout
    startPositionToFill = 0;
    // rather than a loop over the systems, could first determine systemID
    // from cellID, match it against the systemIDs and skip cell if match not found
    // can do it in a separate PR
    for (size_t k = 0; k < m_readoutNames.size(); k++)
    {
      if (k > 0)
        startPositionToFill += m_numLayers[k - 1];
      int systemID = m_systemIDs[k];
      std::string layerField = m_layerFieldNames[k];
      std::string thetaField = m_thetaFieldNames[k];
      std::string moduleField = m_moduleFieldNames[k];
      std::string readoutName = m_readoutNames[k];
      dd4hep::DDSegmentation::BitFieldCoder *decoder = m_geoSvc->getDetector()->readout(readoutName).idSpec().decoder();

      // loop over the cells
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
        double weightLog = std::max(0., m_thetaRecalcLayerWeights[k][layer] + log(eCell / sumEnLayer[layer+startPositionToFill]));
        TVector3 v = TVector3(cell->getPosition().x, cell->getPosition().y, cell->getPosition().z);
        double theta = v.Theta();
        double phi = v.Phi();

        // for clusters that are around the -pi<->pi transition, we want to avoid averaging
        // over phi values that might differ by 2pi. in that case, for cells with negative
        // phi we add two pi, so that we average phi values all close to pi
        if (isClusterPhiNearPi && phi < 0.) {
          phi += TMath::TwoPi();
        }
        if (systemID == systemID_EMB && isResetModuleID && module_id > nModules[k]/2) {
          module_id -= nModules[k];  // transition near 1535..0, reset the module ID
        }

        if (m_thetaRecalcLayerWeights[k][layer]<0)
          sumThetaLayer[layer+startPositionToFill] += (eCell * theta);
        else
          sumThetaLayer[layer+startPositionToFill] += (weightLog * theta);
        sumWeightLayer[layer+startPositionToFill] += weightLog;
        sumPhiLayer[layer+startPositionToFill] += (eCell * phi);

        // do pi0/photon shape var only for EMB
        if (m_do_photon_shapeVar && systemID == systemID_EMB) {
          // E, theta_id, and module_id of cells in layer
          vec_E_cell_layer[layer+startPositionToFill].push_back(eCell);
          vec_theta_cell_layer[layer+startPositionToFill].push_back(theta_id);
          vec_module_cell_layer[layer+startPositionToFill].push_back(module_id);
          // sum them for width in theta/module calculation
	  if (m_do_widthTheta_logE_weights) {
            theta2_E_layer[layer+startPositionToFill] += theta_id * theta_id * weightLog;
            theta_E_layer[layer+startPositionToFill] += theta_id * weightLog;
          } else {
            theta2_E_layer[layer+startPositionToFill] += theta_id * theta_id * eCell;
            theta_E_layer[layer+startPositionToFill] += theta_id * eCell;
          }
          module2_E_layer[layer+startPositionToFill] += module_id * module_id * eCell;
          module_E_layer[layer+startPositionToFill] += module_id * eCell;
        }
      }  // end of loop over cells
    }  // end of loop over each system / readout

    // local maxima (could be more than one) and the corresponding theta
    std::vector<std::pair<std::vector<int>, std::vector<double>>> theta_E_pair;
    std::vector<std::vector<double>> local_E_Max(numLayersTotal, std::vector<double>());
    std::vector<std::vector<int>> local_E_Max_theta(numLayersTotal, std::vector<int>());

    startPositionToFill = 0;
    for (size_t k = 0; k < m_readoutNames.size(); k++) {
      if (!m_do_photon_shapeVar)    break;
      if (m_systemIDs[k] != systemID_EMB)    continue;    // do pi0/photon shape var only for EMB
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

        // fill the zero energy cells in 1D theta-E profile
        for (int i = result.first.front(); i <= result.first.back(); i += nMergedThetaCells[layer+startPositionToFill]) {
          if (std::find(result.first.begin(), result.first.end(), i) == result.first.end()) {
            auto it = std::lower_bound(result.first.begin(), result.first.end(), i);
            int idx = it - result.first.begin();
            result.first.insert(it, i);
            result.second.insert(result.second.begin() + idx, 0);
          }
        }
        theta_E_pair.push_back(result);

        // loop over theta IDs to find the local E maxima
        for (size_t i = 0; i < theta_E_pair[layer+startPositionToFill].second.size(); i ++) {
          //std::cout << i << " " << theta_E_pair[layer+startPositionToFill].first[i] << " " << theta_E_pair[layer+startPositionToFill].second[i] << std::endl;
          if (
            (i == 0 && theta_E_pair[layer+startPositionToFill].second[i] > theta_E_pair[layer+startPositionToFill].second[i+1]) ||
            (i == theta_E_pair[layer+startPositionToFill].second.size()-1 && theta_E_pair[layer+startPositionToFill].second[i] > theta_E_pair[layer+startPositionToFill].second[i-1]) ||
            (i != 0 && i != (theta_E_pair[layer+startPositionToFill].second.size()-1) &&
              theta_E_pair[layer+startPositionToFill].second[i] > theta_E_pair[layer+startPositionToFill].second[i-1] &&
              theta_E_pair[layer+startPositionToFill].second[i] > theta_E_pair[layer+startPositionToFill].second[i+1])
           ) {
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
          // the theta_E_pair are sorted in ascending order of thetaID
          for (size_t i = 0; i < theta_E_pair[layer+startPositionToFill].second.size(); i ++) {
            if (
              (theta_E_pair[layer+startPositionToFill].first[i] > std::min(E_cell_Max_theta[layer+startPositionToFill], E_cell_secMax_theta[layer+startPositionToFill])) &&
              (theta_E_pair[layer+startPositionToFill].first[i] < std::max(E_cell_Max_theta[layer+startPositionToFill], E_cell_secMax_theta[layer+startPositionToFill])) &&
              (theta_E_pair[layer+startPositionToFill].second[i] < E_cell_Min[layer+startPositionToFill])
            )
            {
              E_cell_Min[layer+startPositionToFill] = theta_E_pair[layer+startPositionToFill].second[i];
            }
          }
        }
        if (E_cell_Min[layer+startPositionToFill] > 1e12)    E_cell_Min[layer+startPositionToFill] = 0.;  // check E_cell_Min
      }  // end of loop over layers
    }

    // local maxima of E vs phi (could be more than one) and the corresponding module
    std::vector<std::pair<std::vector<int>, std::vector<double>>> module_E_pair;
    std::vector<std::vector<double>> local_E_Max_vs_phi(numLayersTotal, std::vector<double>());
    std::vector<std::vector<int>> local_E_Max_vs_phi_module(numLayersTotal, std::vector<int>());

    startPositionToFill = 0;
    for (size_t k = 0; k < m_readoutNames.size(); k++) {
      if (!m_do_photon_shapeVar)    break;
      if (m_systemIDs[k] != systemID_EMB)    continue;    // do pi0/photon shape var only for EMB
      if (k > 0)    startPositionToFill += m_numLayers[k - 1];
      // loop over layers
      for (unsigned layer = 0; layer < m_numLayers[k]; layer++) {
        auto result_2 = MergeSumAndSort(vec_module_cell_layer[layer+startPositionToFill], vec_E_cell_layer[layer+startPositionToFill]);
        // fill the zero energy cells in 1D module-E profile
        for (int i = result_2.first.front(); i <= result_2.first.back(); i += nMergedModules[layer+startPositionToFill]) {
          if (std::find(result_2.first.begin(), result_2.first.end(), i) == result_2.first.end()) {
            auto it = std::lower_bound(result_2.first.begin(), result_2.first.end(), i);
            int idx = it - result_2.first.begin();
            result_2.first.insert(it, i);
            result_2.second.insert(result_2.second.begin() + idx, 0);
          }
        }
        module_E_pair.push_back(result_2);

        // loop over module IDs to find the local E maxima
        for (size_t i = 0; i < module_E_pair[layer+startPositionToFill].second.size(); i ++) {
          //std::cout << i << " " << module_E_pair[layer+startPositionToFill].first[i] << " " << module_E_pair[layer+startPositionToFill].second[i] << std::endl;
          if (
            (i == 0 && module_E_pair[layer+startPositionToFill].second[i] > module_E_pair[layer+startPositionToFill].second[i+1]) ||
            (i == module_E_pair[layer+startPositionToFill].second.size()-1 && module_E_pair[layer+startPositionToFill].second[i] > module_E_pair[layer+startPositionToFill].second[i-1]) ||
            (i != 0 && i != (module_E_pair[layer+startPositionToFill].second.size()-1) &&
              module_E_pair[layer+startPositionToFill].second[i] > module_E_pair[layer+startPositionToFill].second[i-1] &&
              module_E_pair[layer+startPositionToFill].second[i] > module_E_pair[layer+startPositionToFill].second[i+1])
           ) {
            local_E_Max_vs_phi[layer+startPositionToFill].push_back(module_E_pair[layer+startPositionToFill].second[i]);
            local_E_Max_vs_phi_module[layer+startPositionToFill].push_back(module_E_pair[layer+startPositionToFill].first[i]);
          }
        }  // end of loop over module IDs

        if (local_E_Max_vs_phi[layer+startPositionToFill].empty()) {
          E_cell_vs_phi_Max[layer+startPositionToFill] = 0.;
          E_cell_vs_phi_secMax[layer+startPositionToFill] = 0.;
          E_cell_vs_phi_Max_module[layer+startPositionToFill] = 0;
          E_cell_vs_phi_secMax_module[layer+startPositionToFill] = 0;
          E_cell_vs_phi_Min[layer+startPositionToFill] = 0.;
        }
        else if (local_E_Max_vs_phi[layer+startPositionToFill].size() < 2) {
          E_cell_vs_phi_Max[layer+startPositionToFill] = local_E_Max_vs_phi[layer+startPositionToFill][0];
          E_cell_vs_phi_secMax[layer+startPositionToFill] = 0.;
          E_cell_vs_phi_Max_module[layer+startPositionToFill] = local_E_Max_vs_phi_module[layer+startPositionToFill][0];
          E_cell_vs_phi_secMax_module[layer+startPositionToFill] = 0;
          E_cell_vs_phi_Min[layer+startPositionToFill] = 0.;
        }
        else {
          std::vector<double> sortedVec = local_E_Max_vs_phi[layer+startPositionToFill];
          // move the top 2 max to the beginning
          std::partial_sort(sortedVec.begin(), sortedVec.begin()+2, sortedVec.end(), std::greater<double>());
          E_cell_vs_phi_Max[layer+startPositionToFill] = sortedVec[0];
          E_cell_vs_phi_secMax[layer+startPositionToFill] = sortedVec[1];
          // get the corresponding module IDs
          auto it_Max   = std::find(local_E_Max_vs_phi[layer+startPositionToFill].begin(), local_E_Max_vs_phi[layer+startPositionToFill].end(), sortedVec[0]);
          int index_Max = std::distance(local_E_Max_vs_phi[layer+startPositionToFill].begin(), it_Max);
          auto it_secMax   = std::find(local_E_Max_vs_phi[layer+startPositionToFill].begin(), local_E_Max_vs_phi[layer+startPositionToFill].end(), sortedVec[1]);
          int index_secMax = std::distance(local_E_Max_vs_phi[layer+startPositionToFill].begin(), it_secMax);
          E_cell_vs_phi_Max_module[layer+startPositionToFill]    = local_E_Max_vs_phi_module[layer+startPositionToFill][index_Max];
          E_cell_vs_phi_secMax_module[layer+startPositionToFill] = local_E_Max_vs_phi_module[layer+startPositionToFill][index_secMax];
          // find the E_min inside the module range of E_cell_Max and E_cell_secMax
          for (size_t i = 0; i < module_E_pair[layer+startPositionToFill].second.size(); i ++) {
            if (
              (module_E_pair[layer+startPositionToFill].first[i] > std::min(E_cell_vs_phi_Max_module[layer+startPositionToFill], E_cell_vs_phi_secMax_module[layer+startPositionToFill])) &&
              (module_E_pair[layer+startPositionToFill].first[i] < std::max(E_cell_vs_phi_Max_module[layer+startPositionToFill], E_cell_vs_phi_secMax_module[layer+startPositionToFill])) &&
              (module_E_pair[layer+startPositionToFill].second[i] < E_cell_vs_phi_Min[layer+startPositionToFill])
            )
            {
              E_cell_vs_phi_Min[layer+startPositionToFill] = module_E_pair[layer+startPositionToFill].second[i];
            }
          }
        }
        if (E_cell_vs_phi_Min[layer+startPositionToFill] > 1e12)  E_cell_vs_phi_Min[layer+startPositionToFill] = 0.;  // check E_cell_Min
      }  // end of loop over layers
    }

    // save energy and theta/phi positions per layer in shape parameters
    startPositionToFill = 0;
    for (size_t k = 0; k < m_readoutNames.size(); k++)
    {
      if (k > 0)
        startPositionToFill += m_numLayers[k - 1];

      int systemID = m_systemIDs[k];
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

        // do pi0/photon shape var only for EMB
        if (m_do_photon_shapeVar && systemID == systemID_EMB) {
          if (m_do_widthTheta_logE_weights) {
            double w_theta2 = theta2_E_layer[layer+startPositionToFill] / sumWeightLayer[layer+startPositionToFill] - std::pow(theta_E_layer[layer+startPositionToFill] / sumWeightLayer[layer+startPositionToFill], 2);
            // Correct very small but negative values caused by computational precision
            if (w_theta2 < 0. && w_theta2 > -1e-9)    w_theta2 = 0. ;
            width_theta[layer+startPositionToFill] = (sumWeightLayer[layer+startPositionToFill] > 0.) ? sqrt(w_theta2) : 0. ;
          } else {
            double w_theta2 = theta2_E_layer[layer+startPositionToFill] / sumEnLayer[layer+startPositionToFill] - std::pow(theta_E_layer[layer+startPositionToFill] / sumEnLayer[layer+startPositionToFill], 2);
            if (w_theta2 < 0. && w_theta2 > -1e-9)    w_theta2 = 0. ;
            width_theta[layer+startPositionToFill] = (sumEnLayer[layer+startPositionToFill] > 0.) ? sqrt(w_theta2) : 0. ;
	  }

          double w_module2 = module2_E_layer[layer+startPositionToFill] / sumEnLayer[layer+startPositionToFill] - std::pow(module_E_layer[layer+startPositionToFill] / sumEnLayer[layer+startPositionToFill], 2);
          // Correct very small but negative values caused by computational precision
          if (w_module2 < 0. && w_module2 > -1e-9)    w_module2 = 0. ;
          width_module[layer+startPositionToFill] = (sumEnLayer[layer+startPositionToFill] > 0.) ? sqrt(w_module2) : 0. ;

          double Ratio_E = (E_cell_Max[layer+startPositionToFill] - E_cell_secMax[layer+startPositionToFill]) /
                           (E_cell_Max[layer+startPositionToFill] + E_cell_secMax[layer+startPositionToFill]);
          if (E_cell_Max[layer+startPositionToFill] + E_cell_secMax[layer+startPositionToFill] == 0)    Ratio_E = 1.;
          Ratio_E_max_2ndmax[layer+startPositionToFill] = Ratio_E;
          Delta_E_2ndmax_min[layer+startPositionToFill] = E_cell_secMax[layer+startPositionToFill] - E_cell_Min[layer+startPositionToFill];

          double Ratio_E_vs_phi = (E_cell_vs_phi_Max[layer+startPositionToFill] - E_cell_vs_phi_secMax[layer+startPositionToFill]) /
                                  (E_cell_vs_phi_Max[layer+startPositionToFill] + E_cell_vs_phi_secMax[layer+startPositionToFill]);
          if (E_cell_vs_phi_Max[layer+startPositionToFill] + E_cell_vs_phi_secMax[layer+startPositionToFill] == 0.)    Ratio_E_vs_phi = 1.;
          Ratio_E_max_2ndmax_vs_phi[layer+startPositionToFill] = Ratio_E_vs_phi;
          Delta_E_2ndmax_min_vs_phi[layer+startPositionToFill] = E_cell_vs_phi_secMax[layer+startPositionToFill] - E_cell_vs_phi_Min[layer+startPositionToFill];

          if (local_E_Max[layer+startPositionToFill].size() > 0) {
            double E_m1 = 0.;
            double E_p1 = 0.;
            int theta_m1 = 0;
            int theta_p1 = 0;
            double E_m2 = 0.;
            double E_p2 = 0.;
            int theta_m2 = 0;
            int theta_p2 = 0;
            double E_m3 = 0.;
            double E_p3 = 0.;
            int theta_m3 = 0;
            int theta_p3 = 0;
            double E_m4 = 0.;
            double E_p4 = 0.;
            int theta_m4 = 0;
            int theta_p4 = 0;
            auto it_1 = std::find(theta_E_pair[layer+startPositionToFill].second.begin(), theta_E_pair[layer+startPositionToFill].second.end(), E_cell_Max[layer+startPositionToFill]);
            int ind_1 = std::distance(theta_E_pair[layer+startPositionToFill].second.begin(), it_1);

            if (ind_1-1 >= 0) {
              E_m1     = theta_E_pair[layer].second[ind_1-1];
              theta_m1 = theta_E_pair[layer].first[ind_1-1];
            } else {
              E_m1     = 0;
              theta_m1 = 0;
            }
            if (static_cast<size_t>(ind_1+1) < theta_E_pair[layer].second.size()) {
              E_p1     = theta_E_pair[layer].second[ind_1+1];
              theta_p1 = theta_E_pair[layer].first[ind_1+1];
            } else {
              E_p1     = 0;
              theta_p1 = 0;
            }

            if (ind_1-2 >= 0) {
              E_m2     = theta_E_pair[layer].second[ind_1-2];
              theta_m2 = theta_E_pair[layer].first[ind_1-2];
            } else {
              E_m2     = 0;
              theta_m2 = 0;
            }
            if (static_cast<size_t>(ind_1+2) < theta_E_pair[layer].second.size()) {
              E_p2     = theta_E_pair[layer].second[ind_1+2];
              theta_p2 = theta_E_pair[layer].first[ind_1+2];
            } else {
              E_p2     = 0;
              theta_p2 = 0;
            }

            if (ind_1-3 >= 0) {
              E_m3     = theta_E_pair[layer].second[ind_1-3];
              theta_m3 = theta_E_pair[layer].first[ind_1-3];
            } else {
              E_m3     = 0;
              theta_m3 = 0;
            }
            if (static_cast<size_t>(ind_1+3) < theta_E_pair[layer].second.size()) {
              E_p3     = theta_E_pair[layer].second[ind_1+3];
              theta_p3 = theta_E_pair[layer].first[ind_1+3];
            } else {
              E_p3     = 0;
              theta_p3 = 0;
            }
 
            if (ind_1-4 >= 0) {
              E_m4     = theta_E_pair[layer].second[ind_1-4];
              theta_m4 = theta_E_pair[layer].first[ind_1-4];
            } else {
              E_m4     = 0;
              theta_m4 = 0;
            }
            if (static_cast<size_t>(ind_1+4) < theta_E_pair[layer].second.size()) {
              E_p4     = theta_E_pair[layer].second[ind_1+4];
              theta_p4 = theta_E_pair[layer].first[ind_1+4];
            } else {
              E_p4     = 0;
              theta_p4 = 0;
            }
            // sum_E_nBin is also used in E fraction side calculation so put that outside the if block
            double sum_E_3Bin = E_m1 + E_cell_Max[layer+startPositionToFill] + E_p1;
            double sum_E_5Bin = sum_E_3Bin + E_m2 + E_p2;
            double sum_E_7Bin = sum_E_5Bin + E_m3 + E_p3;
            double sum_E_9Bin = sum_E_7Bin + E_m4 + E_p4;

            if (m_do_widthTheta_logE_weights) {
              double weightLog_E_max = std::max(0., m_thetaRecalcLayerWeights[k][layer] + log(E_cell_Max[layer+startPositionToFill] / sumEnLayer[layer+startPositionToFill]));
              double weightLog_m1 = std::max(0., m_thetaRecalcLayerWeights[k][layer] + log(E_m1 / sumEnLayer[layer+startPositionToFill]));
              double weightLog_m2 = std::max(0., m_thetaRecalcLayerWeights[k][layer] + log(E_m2 / sumEnLayer[layer+startPositionToFill]));
              double weightLog_m3 = std::max(0., m_thetaRecalcLayerWeights[k][layer] + log(E_m3 / sumEnLayer[layer+startPositionToFill]));
              double weightLog_m4 = std::max(0., m_thetaRecalcLayerWeights[k][layer] + log(E_m4 / sumEnLayer[layer+startPositionToFill]));
              double weightLog_p1 = std::max(0., m_thetaRecalcLayerWeights[k][layer] + log(E_p1 / sumEnLayer[layer+startPositionToFill]));
              double weightLog_p2 = std::max(0., m_thetaRecalcLayerWeights[k][layer] + log(E_p2 / sumEnLayer[layer+startPositionToFill]));
              double weightLog_p3 = std::max(0., m_thetaRecalcLayerWeights[k][layer] + log(E_p3 / sumEnLayer[layer+startPositionToFill]));
              double weightLog_p4 = std::max(0., m_thetaRecalcLayerWeights[k][layer] + log(E_p4 / sumEnLayer[layer+startPositionToFill]));

              double sum_weightLog_3Bin = weightLog_E_max + weightLog_m1 + weightLog_p1;
              double sum_weightLog_5Bin = sum_weightLog_3Bin + weightLog_m2 + weightLog_p2;
              double sum_weightLog_7Bin = sum_weightLog_5Bin + weightLog_m3 + weightLog_p3;
              double sum_weightLog_9Bin = sum_weightLog_7Bin + weightLog_m4 + weightLog_p4;

              double theta2_E_3Bin = theta_m1 * theta_m1 * weightLog_m1 + E_cell_Max_theta[layer+startPositionToFill] * E_cell_Max_theta[layer+startPositionToFill] * weightLog_E_max + theta_p1 * theta_p1 * weightLog_p1;
              double theta_E_3Bin = theta_m1 * weightLog_m1 + E_cell_Max_theta[layer+startPositionToFill] * weightLog_E_max + theta_p1 * weightLog_p1;
              double theta2_E_5Bin = theta2_E_3Bin + theta_m2 * theta_m2 * weightLog_m2 + theta_p2 * theta_p2 * weightLog_p2;
              double theta_E_5Bin = theta_E_3Bin + theta_m2 * weightLog_m2 + theta_p2 * weightLog_p2;
              double theta2_E_7Bin = theta2_E_5Bin + theta_m3 * theta_m3 * weightLog_m3 + theta_p3 * theta_p3 * weightLog_p3;
              double theta_E_7Bin = theta_E_5Bin + theta_m3 * weightLog_m3 + theta_p3 * weightLog_p3;
              double theta2_E_9Bin = theta2_E_7Bin + theta_m4 * theta_m4 * weightLog_m4 + theta_p4 * theta_p4 * weightLog_p4;
              double theta_E_9Bin = theta_E_7Bin + theta_m4 * weightLog_m4 + theta_p4 * weightLog_p4;

              double _w_theta_3Bin2 = theta2_E_3Bin / sum_weightLog_3Bin - std::pow(theta_E_3Bin / sum_weightLog_3Bin, 2);
              double _w_theta_5Bin2 = theta2_E_5Bin / sum_weightLog_5Bin - std::pow(theta_E_5Bin / sum_weightLog_5Bin, 2);
              double _w_theta_7Bin2 = theta2_E_7Bin / sum_weightLog_7Bin - std::pow(theta_E_7Bin / sum_weightLog_7Bin, 2);
              double _w_theta_9Bin2 = theta2_E_9Bin / sum_weightLog_9Bin - std::pow(theta_E_9Bin / sum_weightLog_9Bin, 2);
              // Correct very small but negative values caused by computational precision
              if (_w_theta_3Bin2 < 0. && _w_theta_3Bin2 > -1e-9)    _w_theta_3Bin2 = 0. ;
              if (_w_theta_5Bin2 < 0. && _w_theta_5Bin2 > -1e-9)    _w_theta_5Bin2 = 0. ;
              if (_w_theta_7Bin2 < 0. && _w_theta_7Bin2 > -1e-9)    _w_theta_7Bin2 = 0. ;
              if (_w_theta_9Bin2 < 0. && _w_theta_9Bin2 > -1e-9)    _w_theta_9Bin2 = 0. ;

              width_theta_3Bin[layer+startPositionToFill] = (sum_weightLog_3Bin > 0.) ? sqrt(_w_theta_3Bin2) : 0. ;
              width_theta_5Bin[layer+startPositionToFill] = (sum_weightLog_5Bin > 0.) ? sqrt(_w_theta_5Bin2) : 0. ;
              width_theta_7Bin[layer+startPositionToFill] = (sum_weightLog_7Bin > 0.) ? sqrt(_w_theta_7Bin2) : 0. ;
              width_theta_9Bin[layer+startPositionToFill] = (sum_weightLog_9Bin > 0.) ? sqrt(_w_theta_9Bin2) : 0. ;
            } else {
              double theta2_E_3Bin = theta_m1 * theta_m1 * E_m1 + E_cell_Max_theta[layer+startPositionToFill] * E_cell_Max_theta[layer+startPositionToFill] * E_cell_Max[layer+startPositionToFill] + theta_p1 * theta_p1 * E_p1;
              double theta_E_3Bin = theta_m1 * E_m1 + E_cell_Max_theta[layer+startPositionToFill] * E_cell_Max[layer+startPositionToFill] + theta_p1 * E_p1;
              double theta2_E_5Bin = theta2_E_3Bin + theta_m2 * theta_m2 * E_m2 + theta_p2 * theta_p2 * E_p2;
              double theta_E_5Bin = theta_E_3Bin + theta_m2 * E_m2 + theta_p2 * E_p2;
              double theta2_E_7Bin = theta2_E_5Bin + theta_m3 * theta_m3 * E_m3 + theta_p3 * theta_p3 * E_p3;
              double theta_E_7Bin = theta_E_5Bin + theta_m3 * E_m3 + theta_p3 * E_p3;
              double theta2_E_9Bin = theta2_E_7Bin + theta_m4 * theta_m4 * E_m4 + theta_p4 * theta_p4 * E_p4;
              double theta_E_9Bin = theta_E_7Bin + theta_m4 * E_m4 + theta_p4 * E_p4;

              double _w_theta_3Bin2 = theta2_E_3Bin / sum_E_3Bin - std::pow(theta_E_3Bin / sum_E_3Bin, 2);
              double _w_theta_5Bin2 = theta2_E_5Bin / sum_E_5Bin - std::pow(theta_E_5Bin / sum_E_5Bin, 2);
              double _w_theta_7Bin2 = theta2_E_7Bin / sum_E_7Bin - std::pow(theta_E_7Bin / sum_E_7Bin, 2);
              double _w_theta_9Bin2 = theta2_E_9Bin / sum_E_9Bin - std::pow(theta_E_9Bin / sum_E_9Bin, 2);
              // Correct very small but negative values caused by computational precision
              if (_w_theta_3Bin2 < 0. && _w_theta_3Bin2 > -1e-9)    _w_theta_3Bin2 = 0. ;
              if (_w_theta_5Bin2 < 0. && _w_theta_5Bin2 > -1e-9)    _w_theta_5Bin2 = 0. ;
              if (_w_theta_7Bin2 < 0. && _w_theta_7Bin2 > -1e-9)    _w_theta_7Bin2 = 0. ;
              if (_w_theta_9Bin2 < 0. && _w_theta_9Bin2 > -1e-9)    _w_theta_9Bin2 = 0. ;

              width_theta_3Bin[layer+startPositionToFill] = (sum_E_3Bin > 0.) ? sqrt(_w_theta_3Bin2) : 0. ;
              width_theta_5Bin[layer+startPositionToFill] = (sum_E_5Bin > 0.) ? sqrt(_w_theta_5Bin2) : 0. ;
              width_theta_7Bin[layer+startPositionToFill] = (sum_E_7Bin > 0.) ? sqrt(_w_theta_7Bin2) : 0. ;
              width_theta_9Bin[layer+startPositionToFill] = (sum_E_9Bin > 0.) ? sqrt(_w_theta_9Bin2) : 0. ;
            }
            E_fr_side_pm2[layer+startPositionToFill] = (sum_E_3Bin > 0.) ? (sum_E_5Bin / sum_E_3Bin - 1.) : 0. ;
            E_fr_side_pm3[layer+startPositionToFill] = (sum_E_3Bin > 0.) ? (sum_E_7Bin / sum_E_3Bin - 1.) : 0. ;
            E_fr_side_pm4[layer+startPositionToFill] = (sum_E_3Bin > 0.) ? (sum_E_9Bin / sum_E_3Bin - 1.) : 0. ;
          } else {
            width_theta_3Bin[layer+startPositionToFill] = 0.;
            width_theta_5Bin[layer+startPositionToFill] = 0.;
            width_theta_7Bin[layer+startPositionToFill] = 0.;
            width_theta_9Bin[layer+startPositionToFill] = 0.;
            E_fr_side_pm2[layer+startPositionToFill] = 0.;
            E_fr_side_pm3[layer+startPositionToFill] = 0.;
            E_fr_side_pm4[layer+startPositionToFill] = 0.;
          }
          newCluster.addToShapeParameters(maxCellEnergyInLayer[layer+startPositionToFill]);
          newCluster.addToShapeParameters(width_theta[layer+startPositionToFill]);
          newCluster.addToShapeParameters(width_module[layer+startPositionToFill]);
          newCluster.addToShapeParameters(Ratio_E_max_2ndmax[layer+startPositionToFill]);
          newCluster.addToShapeParameters(Delta_E_2ndmax_min[layer+startPositionToFill]);
          newCluster.addToShapeParameters(Ratio_E_max_2ndmax_vs_phi[layer+startPositionToFill]);
          newCluster.addToShapeParameters(Delta_E_2ndmax_min_vs_phi[layer+startPositionToFill]);
          newCluster.addToShapeParameters(width_theta_3Bin[layer+startPositionToFill]);
          newCluster.addToShapeParameters(width_theta_5Bin[layer+startPositionToFill]);
          newCluster.addToShapeParameters(width_theta_7Bin[layer+startPositionToFill]);
          newCluster.addToShapeParameters(width_theta_9Bin[layer+startPositionToFill]);
          newCluster.addToShapeParameters(E_fr_side_pm2[layer+startPositionToFill]);
          newCluster.addToShapeParameters(E_fr_side_pm3[layer+startPositionToFill]);
          newCluster.addToShapeParameters(E_fr_side_pm4[layer+startPositionToFill]);
        }
      }  // end of loop over layers
    }  // end of loop over system/readout
    newCluster.addToShapeParameters(p4cl.M());
    newCluster.addToShapeParameters(nCells);
  }  // end of loop over clusters

  return StatusCode::SUCCESS;
}
