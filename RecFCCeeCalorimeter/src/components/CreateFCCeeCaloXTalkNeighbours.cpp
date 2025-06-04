#include "CreateFCCeeCaloXTalkNeighbours.h"

// DD4hep
#include "DD4hep/Detector.h"

// k4FWCore
#include "k4Interface/IGeoSvc.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"
#include "detectorSegmentations/FCCSWGridModuleThetaMerged_k4geo.h"
#include "detectorSegmentations/FCCSWGridPhiTheta_k4geo.h"
#include "detectorSegmentations/GridTheta_k4geo.h"

// ROOT
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

DECLARE_COMPONENT(CreateFCCeeCaloXTalkNeighbours)

CreateFCCeeCaloXTalkNeighbours::CreateFCCeeCaloXTalkNeighbours(const std::string& aName, ISvcLocator* aSL)
    : base_class(aName, aSL) {
  declareProperty("outputFileName", m_outputFileName, "Name of the output file");
}

CreateFCCeeCaloXTalkNeighbours::~CreateFCCeeCaloXTalkNeighbours() {}

StatusCode CreateFCCeeCaloXTalkNeighbours::initialize() {
  // Initialize necessary Gaudi components
  if (Service::initialize().isFailure()) {
    error() << "Unable to initialize Service()" << endmsg;
    return StatusCode::FAILURE;
  }
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service. "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }
  if(m_connectBarrels){
    warning() << "connectBarrels feature is not yet implemented!" <<endmsg;
  }
  std::unordered_map<uint64_t, std::vector<std::pair<uint64_t, double>>> map;

  for (uint iSys = 0; iSys < m_readoutNamesSegmented.size(); iSys++) {
    // Check if readout exists
    info() << "Readout: " << m_readoutNamesSegmented[iSys] << endmsg;
    if (m_geoSvc->getDetector()->readouts().find(m_readoutNamesSegmented[iSys]) ==
        m_geoSvc->getDetector()->readouts().end()) {
      error() << "Readout <<" << m_readoutNamesSegmented[iSys] << ">> does not exist." << endmsg;
      return StatusCode::FAILURE;
    }

    // get theta-based segmentation
    dd4hep::DDSegmentation::Segmentation* aSegmentation =
        m_geoSvc->getDetector()->readout(m_readoutNamesSegmented[iSys]).segmentation().segmentation();
    if (aSegmentation == nullptr) {
      error() << "Segmentation does not exist." << endmsg;
      return StatusCode::FAILURE;
    }

    std::string segmentationType = aSegmentation->type();
    info() << "Segmentation type : " << segmentationType << endmsg;

    dd4hep::DDSegmentation::GridTheta_k4geo* segmentation = nullptr;
    dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo* moduleThetaSegmentation = nullptr;
    if (segmentationType == "FCCSWGridModuleThetaMerged_k4geo") {
      segmentation = dynamic_cast<dd4hep::DDSegmentation::GridTheta_k4geo*>(aSegmentation);
      moduleThetaSegmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo*>(aSegmentation);
    } else {
      error() << "Segmentation type not handled." << endmsg;
      return StatusCode::FAILURE;
    }

    if (segmentation == nullptr || moduleThetaSegmentation == nullptr) {
      error() << "Unable to cast segmentation pointer!!!!" << endmsg;
      return StatusCode::FAILURE;
    }

    info() << "Segmentation: size in Theta " << segmentation->gridSizeTheta() << endmsg;
    info() << "Segmentation: offset in Theta " << segmentation->offsetTheta() << endmsg;
    if (segmentationType == "FCCSWGridModuleThetaMerged_k4geo") {
      info() << "Segmentation: bins in Module " << moduleThetaSegmentation->nModules() << endmsg;
    }

    // retrieve decoders for the crosstalk neighbour finding
    auto decoder = m_geoSvc->getDetector()->readout(m_readoutNamesSegmented[iSys]).idSpec().decoder();

    // Loop over all cells in the calorimeter and retrieve existing cellIDs and find neihbours
    if (segmentationType == "FCCSWGridModuleThetaMerged_k4geo") {
      std::vector<std::vector<std::pair<int, int>>> extrema_layer;
      for (unsigned int ilayer = 0; ilayer < m_activeVolumesNumbersSegmented[iSys]; ilayer++) {
        dd4hep::DDSegmentation::CellID volumeId = 0;
        static_cast<dd4hep::DDSegmentation::BitFieldCoder>(*decoder)[m_fieldNamesSegmented[iSys]].set(
            volumeId, m_fieldValuesSegmented[iSys]);
        static_cast<dd4hep::DDSegmentation::BitFieldCoder>(*decoder)[m_activeFieldNamesSegmented[iSys]].set(volumeId,
                                                                                                            ilayer);
        static_cast<dd4hep::DDSegmentation::BitFieldCoder>(*decoder)["theta"].set(volumeId, 0);
        static_cast<dd4hep::DDSegmentation::BitFieldCoder>(*decoder)["module"].set(volumeId, 0);
        auto numCells = det::utils::numberOfCells(volumeId, *moduleThetaSegmentation);

        std::vector<std::pair<int, int>> extrema;
        // extrema[0]: min layer, n layers
        extrema.emplace_back(std::make_pair(0, m_activeVolumesNumbersSegmented[iSys] - 1));
        // extrema[1]: min module number (0), max module number
        extrema.emplace_back(std::make_pair(0, (numCells[0] - 1) * moduleThetaSegmentation->mergedModules(ilayer)));
        // extrema[2]: min theta ID, n (merged) theta cells
        extrema.emplace_back(std::make_pair(
            numCells[2], numCells[2] + (numCells[1] - 1) * moduleThetaSegmentation->mergedThetaCells(ilayer)));
        debug() << "Layer: " << ilayer << endmsg;
        debug() << "Extrema[0]: " << extrema[0].first << " , " << extrema[0].second << endmsg;
        debug() << "Extrema[1]: " << extrema[1].first << " , " << extrema[1].second << endmsg;
        debug() << "Extrema[2]: " << extrema[2].first << " , " << extrema[2].second << endmsg;
        debug() << "Number of segmentation cells in (module,theta): " << numCells << endmsg;
        extrema_layer.emplace_back(extrema);
      }

      for (unsigned int ilayer = 0; ilayer < m_activeVolumesNumbersSegmented[iSys]; ilayer++) {
        dd4hep::DDSegmentation::CellID volumeId = 0;
        // Get VolumeID
        static_cast<dd4hep::DDSegmentation::BitFieldCoder>(*decoder)[m_fieldNamesSegmented[iSys]].set(
            volumeId, m_fieldValuesSegmented[iSys]);
        static_cast<dd4hep::DDSegmentation::BitFieldCoder>(*decoder)[m_activeFieldNamesSegmented[iSys]].set(volumeId,
                                                                                                            ilayer);
        static_cast<dd4hep::DDSegmentation::BitFieldCoder>(*decoder)["theta"].set(volumeId, 0);
        static_cast<dd4hep::DDSegmentation::BitFieldCoder>(*decoder)["module"].set(volumeId, 0);
        // Loop over segmentation cells to find crosstalk neighbours in ECAL
        for (int imodule = extrema_layer[ilayer][1].first; imodule <= extrema_layer[ilayer][1].second;
             imodule += moduleThetaSegmentation->mergedModules(ilayer)) {
          for (int itheta = extrema_layer[ilayer][2].first; itheta <= extrema_layer[ilayer][2].second;
               itheta += moduleThetaSegmentation->mergedThetaCells(ilayer)) {
            dd4hep::DDSegmentation::CellID cellId = volumeId;
            decoder->set(cellId, "module", imodule);
            decoder->set(cellId, "theta", itheta); // start from the minimum existing theta cell in this layer
            uint64_t id = cellId;
            map.insert(std::pair<uint64_t, std::vector<std::pair<uint64_t, double>>>(
                id, det::crosstalk::getNeighboursModuleThetaMerged(
                        *moduleThetaSegmentation, *decoder, {m_activeFieldNamesSegmented[iSys], "module", "theta"},
                        extrema_layer, id, m_xtalk_coef_radial, m_xtalk_coef_theta, m_xtalk_coef_diagonal,
			m_xtalk_coef_tower)));
          }
        }
      }
    }

    if (msgLevel() <= MSG::DEBUG) {
      std::vector<int> counter;
      counter.assign(40, 0);
      for (const auto& item : map) {
        counter[item.second.size()]++;
      }
      for (uint iCount = 0; iCount < counter.size(); iCount++) {
        if (counter[iCount] != 0) {
          info() << counter[iCount] << " cells have " << iCount << " neighbours" << endmsg;
        }
      }
    }
    info() << "total number of cells:  " << map.size() << endmsg;
  }

  if (msgLevel() <= MSG::DEBUG) {
    std::vector<int> counter;
    counter.assign(40, 0);
    for (const auto& item : map) {
      counter[item.second.size()]++;
    }
    for (uint iCount = 0; iCount < counter.size(); iCount++) {
      if (counter[iCount] != 0) {
        debug() << counter[iCount] << " cells have " << iCount << " neighbours" << endmsg;
      }
    }
  }
  // debug() << "cells with neighbours across Calo boundaries: " << count << endmsg;

  // Check if output directory exists
  std::string outDirPath = gSystem->DirName(m_outputFileName.c_str());
  if (!gSystem->OpenDirectory(outDirPath.c_str())) {
    error() << "Output directory \"" << outDirPath << "\" does not exists! Please create it." << endmsg;
    return StatusCode::FAILURE;
  }

  std::unique_ptr<TFile> outFile(TFile::Open(m_outputFileName.c_str(), "RECREATE"));
  outFile->cd();
  TTree tree("crosstalk_neighbours", "Tree with map of neighbours");
  uint64_t saveCellId;
  std::vector<uint64_t> saveNeighbours;
  std::vector<double> saveCrosstalks;
  tree.Branch("cellId", &saveCellId, "cellId/l");
  tree.Branch("list_crosstalk_neighbours", &saveNeighbours);
  tree.Branch("list_crosstalks", &saveCrosstalks);

  // Debug: save cell position
  std::vector<int> saveCellInfo;
  if (m_debugCellInfo){
    tree.Branch("CellInfo", &saveCellInfo);
  }

  int count_map = 0;
  for (const auto& item : map) {
    saveCellId = item.first;
    std::vector<std::pair<uint64_t, double>> temp_pair = item.second;
    std::vector<uint64_t> temp_neighbours;
    std::vector<double> temp_crosstalks;
    for (const auto& this_temp : temp_pair) {
      temp_neighbours.push_back(this_temp.first);
      temp_crosstalks.push_back(this_temp.second);
    }
    saveNeighbours = temp_neighbours;
    saveCrosstalks = temp_crosstalks;
    // Debug: save cell position
    if (m_debugCellInfo){
      for (uint iSys = 0; iSys < m_readoutNamesSegmented.size(); iSys++) {
        dd4hep::DDSegmentation::Segmentation* aSegmentation =
            m_geoSvc->getDetector()->readout(m_readoutNamesSegmented[iSys]).segmentation().segmentation();
        std::string segmentationType = aSegmentation->type();
        dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo* moduleThetaSegmentation = nullptr;
        auto decoder = m_geoSvc->getDetector()->readout(m_readoutNamesSegmented[iSys]).idSpec().decoder();
        if (segmentationType == "FCCSWGridModuleThetaMerged_k4geo") {
          moduleThetaSegmentation =
              dynamic_cast<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo*>(aSegmentation);
        }
        saveCellInfo = det::crosstalk::getCellIndices(*moduleThetaSegmentation, *decoder,
                                                      {m_activeFieldNamesSegmented[iSys], "module", "theta"}, saveCellId);
      }
    }
    tree.Fill();
    count_map++;
    if (!count_map % 1000)
      std::cout << "Number of cells: " << count_map << std::endl;
  }
  outFile->Write();
  outFile->Close();

  return StatusCode::SUCCESS;
}

StatusCode CreateFCCeeCaloXTalkNeighbours::finalize() { return Service::finalize(); }
