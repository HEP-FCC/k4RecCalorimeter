#include "CreateFCCeeCaloNeighbours.h"

// DD4hep
#include "DD4hep/Detector.h"

// k4FWCore
#include "k4Interface/IGeoSvc.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"
#include "detectorSegmentations/FCCSWEndcapTurbine_k4geo.h"
#include "detectorSegmentations/FCCSWGridModuleThetaMerged_k4geo.h"
#include "detectorSegmentations/FCCSWGridPhiTheta_k4geo.h"
#include "detectorSegmentations/FCCSWHCalPhiRow_k4geo.h"
#include "detectorSegmentations/FCCSWHCalPhiTheta_k4geo.h"
#include "detectorSegmentations/GridTheta_k4geo.h"

// ROOT
#include "TFile.h"
#include "TH1F.h"
#include "TSystem.h"
#include "TTree.h"

DECLARE_COMPONENT(CreateFCCeeCaloNeighbours)

CreateFCCeeCaloNeighbours::CreateFCCeeCaloNeighbours(const std::string& aName, ISvcLocator* aSL)
    : base_class(aName, aSL) {
  declareProperty("outputFileName", m_outputFileName, "Name of the output file");
}

CreateFCCeeCaloNeighbours::~CreateFCCeeCaloNeighbours() {}

StatusCode CreateFCCeeCaloNeighbours::initialize() {
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
  std::unordered_map<uint64_t, std::vector<uint64_t>> map;

  // will be used for volume connecting
  int eCalLastLayer;
  std::pair<int, int> extremaECalLastLayerPhi;
  std::pair<int, int> extremaECalLastLayerModule;
  std::pair<int, int> extremaECalLastLayerTheta;
  int ecalBarrelId = -1; // index of ecal barrel in segmentation list
  dd4hep::DDSegmentation::BitFieldCoder* decoderECalBarrel = nullptr;
  dd4hep::DDSegmentation::FCCSWHCalPhiTheta_k4geo* hcalBarrelSegmentation = nullptr;
  dd4hep::DDSegmentation::FCCSWHCalPhiTheta_k4geo* hcalEndcapSegmentation = nullptr;
  dd4hep::DDSegmentation::FCCSWHCalPhiRow_k4geo* hcalBarrelPhiRowSegmentation = nullptr;
  dd4hep::DDSegmentation::FCCSWHCalPhiRow_k4geo* hcalEndcapPhiRowSegmentation = nullptr;

  //  will be used for volume connecting
  std::pair<int, int> extremaHCalFirstLayerPhi;
  std::pair<int, int> extremaHCalFirstLayerTheta;
  int hcalBarrelId = -1; // index of hcal barrel in segmentation list
  int hcalEndcapId = -1; // index of hcal endcap in segmentation list
  dd4hep::DDSegmentation::BitFieldCoder* decoderHCalBarrel = nullptr;
  dd4hep::DDSegmentation::BitFieldCoder* decoderHCalEndcap = nullptr;
  dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo* ecalBarrelModuleThetaSegmentation = nullptr;
  dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo* ecalBarrelPhiThetaSegmentation = nullptr;
  dd4hep::DDSegmentation::FCCSWEndcapTurbine_k4geo* ecalEndcapTurbineSegmentation = nullptr;

  // precalculate lookup tables for ECAL barrel-endcap connection
  if (m_connectECal) {
    info() << "Will connect ECal barrel and endcap cells. Precalculating lookup tables needed for connection" << endmsg;
    StatusCode st = initialize_lookups();
    if (st != StatusCode::SUCCESS)
      return st;
  }

  for (uint iSys = 0; iSys < m_readoutNamesSegmented.size(); iSys++) {
    // Check if readout exists
    info() << "Creating neighbours for readout: " << m_readoutNamesSegmented[iSys] << endmsg;
    if (m_geoSvc->getDetector()->readouts().find(m_readoutNamesSegmented[iSys]) ==
        m_geoSvc->getDetector()->readouts().end()) {
      error() << "Readout <<" << m_readoutNamesSegmented[iSys] << ">> does not exist." << endmsg;
      return StatusCode::FAILURE;
    }

    // get segmentation
    dd4hep::DDSegmentation::Segmentation* aSegmentation =
        m_geoSvc->getDetector()->readout(m_readoutNamesSegmented[iSys]).segmentation().segmentation();
    if (aSegmentation == nullptr) {
      error() << "Segmentation does not exist." << endmsg;
      return StatusCode::FAILURE;
    }

    std::string segmentationType = aSegmentation->type();
    info() << "Segmentation: type " << segmentationType << endmsg;

    dd4hep::DDSegmentation::GridTheta_k4geo* segmentation = nullptr;
    dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo* phiThetaSegmentation = nullptr;
    dd4hep::DDSegmentation::FCCSWHCalPhiTheta_k4geo* hcalPhiThetaSegmentation = nullptr;
    dd4hep::DDSegmentation::FCCSWHCalPhiRow_k4geo* hcalPhiRowSegmentation = nullptr;
    dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo* moduleThetaSegmentation = nullptr;
    if (segmentationType == "FCCSWGridModuleThetaMerged_k4geo") {
      segmentation = dynamic_cast<dd4hep::DDSegmentation::GridTheta_k4geo*>(aSegmentation);
      moduleThetaSegmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo*>(aSegmentation);
    } else if (segmentationType == "FCCSWGridPhiTheta_k4geo") {
      segmentation = dynamic_cast<dd4hep::DDSegmentation::GridTheta_k4geo*>(aSegmentation);
      phiThetaSegmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo*>(aSegmentation);
    } else if (segmentationType == "FCCSWHCalPhiTheta_k4geo") {
      segmentation = dynamic_cast<dd4hep::DDSegmentation::GridTheta_k4geo*>(aSegmentation);
      hcalPhiThetaSegmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWHCalPhiTheta_k4geo*>(aSegmentation);
    } else if (segmentationType == "FCCSWHCalPhiRow_k4geo") {
      hcalPhiRowSegmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWHCalPhiRow_k4geo*>(aSegmentation);
    } else if (segmentationType == "FCCSWEndcapTurbine_k4geo")

    {
      ecalEndcapTurbineSegmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWEndcapTurbine_k4geo*>(aSegmentation);
    } else {
      error() << "Segmentation type not handled." << endmsg;
      return StatusCode::FAILURE;
    }

    if (((segmentation == nullptr && segmentationType != "FCCSWHCalPhiRow_k4geo") ||
         (phiThetaSegmentation == nullptr && hcalPhiThetaSegmentation == nullptr && hcalPhiRowSegmentation == nullptr &&
          moduleThetaSegmentation == nullptr)) &&
        (ecalEndcapTurbineSegmentation == nullptr)) {
      error() << "Unable to cast segmentation pointer!!!!" << endmsg;
      return StatusCode::FAILURE;
    }

    if ((segmentationType != "FCCSWHCalPhiRow_k4geo") && (segmentationType != "FCCSWEndcapTurbine_k4geo")) {
      info() << "Segmentation: size in Theta " << segmentation->gridSizeTheta() << endmsg;
      info() << "Segmentation: offset in Theta " << segmentation->offsetTheta() << endmsg;
    }
    if (segmentationType == "FCCSWGridModuleThetaMerged_k4geo") {
      info() << "Segmentation: bins in Module " << moduleThetaSegmentation->nModules() << endmsg;
    } else if (segmentationType == "FCCSWGridPhiTheta_k4geo") {
      info() << "Segmentation: size in Phi " << phiThetaSegmentation->gridSizePhi() << endmsg;
      info() << "Segmentation: offset in Phi " << phiThetaSegmentation->offsetPhi() << endmsg;
    } else if (segmentationType == "FCCSWHCalPhiTheta_k4geo") {
      info() << "Segmentation: size in Phi " << hcalPhiThetaSegmentation->gridSizePhi() << endmsg;
      info() << "Segmentation: offset in Phi " << hcalPhiThetaSegmentation->offsetPhi() << endmsg;
    } else if (segmentationType == "FCCSWEndcapTurbine_k4geo") {
      for (int iWheel = 0; iWheel < 3; iWheel++) {
        info() << "Segmentation: nModules for wheel " << iWheel << ": "
               << ecalEndcapTurbineSegmentation->nModules(iWheel) << endmsg;
        info() << "Segmentation: size in rho for wheel " << iWheel << ": "
               << ecalEndcapTurbineSegmentation->gridSizeRho(iWheel) << endmsg;
        info() << "Segmentation: size in z for wheel " << iWheel << ": "
               << ecalEndcapTurbineSegmentation->gridSizeZ(iWheel) << endmsg;
      }
    }

    // retrieve decoders and other info needed for volume (endcap-barrel) connection
    auto decoder = m_geoSvc->getDetector()->readout(m_readoutNamesSegmented[iSys]).idSpec().decoder();
    if (m_connectBarrels) {
      if (m_fieldNamesSegmented[iSys] == "system" && m_fieldValuesSegmented[iSys] == m_ecalBarrelSysId) {
        decoderECalBarrel = decoder;
        ecalBarrelId = iSys;
        if (segmentationType == "FCCSWGridModuleThetaMerged_k4geo") {
          ecalBarrelModuleThetaSegmentation = moduleThetaSegmentation;
        } else if (segmentationType == "FCCSWGridPhiTheta_k4geo") {
          ecalBarrelPhiThetaSegmentation = phiThetaSegmentation;
        } else {
          error() << "ECAL barrel segmentation type not handled for connectBarrels option." << endmsg;
          return StatusCode::FAILURE;
        }
      }

      if (m_fieldNamesSegmented[iSys] == "system" && m_fieldValuesSegmented[iSys] == m_hcalBarrelSysId) {
        decoderHCalBarrel = decoder;
        hcalBarrelId = iSys;
        if (segmentationType == "FCCSWHCalPhiTheta_k4geo") {
          hcalBarrelSegmentation = hcalPhiThetaSegmentation;
        } else if (segmentationType == "FCCSWHCalPhiRow_k4geo") {
          hcalBarrelPhiRowSegmentation = hcalPhiRowSegmentation;
        } else {
          error() << "HCAL barrel segmentation type not handled for connectBarrels option." << endmsg;
          return StatusCode::FAILURE;
        }
      }
    }

    // retrieve decoders and other info needed for HCal Barrel and Endcap connection
    if (m_connectHCal) {
      if (m_fieldNamesSegmented[iSys] == "system" && m_fieldValuesSegmented[iSys] == m_hcalBarrelSysId) {
        decoderHCalBarrel = decoder;
        hcalBarrelId = iSys;
        if (segmentationType == "FCCSWHCalPhiTheta_k4geo") {
          hcalBarrelSegmentation = hcalPhiThetaSegmentation;
        } else if (segmentationType == "FCCSWHCalPhiRow_k4geo") {
          hcalBarrelPhiRowSegmentation = hcalPhiRowSegmentation;
        } else {
          error() << "HCAL barrel segmentation type not handled to connect Barrel and Endcap." << endmsg;
          return StatusCode::FAILURE;
        }
      }
      if (m_fieldNamesSegmented[iSys] == "system" && m_fieldValuesSegmented[iSys] == m_hcalEndcapSysId) {
        decoderHCalEndcap = decoder;
        hcalEndcapId = iSys;
        if (segmentationType == "FCCSWHCalPhiTheta_k4geo") {
          hcalEndcapSegmentation = hcalPhiThetaSegmentation;
        } else if (segmentationType == "FCCSWHCalPhiRow_k4geo") {
          hcalEndcapPhiRowSegmentation = hcalPhiRowSegmentation;
        } else {
          error() << "HCAL Endcap segmentation type not handled to connect Barrel and Endcap." << endmsg;
          return StatusCode::FAILURE;
        }
      }
    }

    // Loop over all cells in the calorimeter and retrieve existing cellIDs and find neihbours
    if (segmentationType == "FCCSWGridPhiTheta_k4geo") {
      // Loop over active layers
      std::vector<std::pair<int, int>> extrema;
      extrema.push_back(std::make_pair(0, m_activeVolumesNumbersSegmented[iSys] - 1));
      extrema.push_back(std::make_pair(0, 0));
      extrema.push_back(std::make_pair(0, 0));
      for (unsigned int ilayer = 0; ilayer < m_activeVolumesNumbersSegmented[iSys]; ilayer++) {
        dd4hep::DDSegmentation::CellID volumeId = 0;
        // Get VolumeID
        // for ECAL OK (volume extends along z)
        // for HCAL not OK.. need to find HCalLayerVol..
        (*decoder)[m_fieldNamesSegmented[iSys]].set(volumeId, m_fieldValuesSegmented[iSys]);
        (*decoder)[m_activeFieldNamesSegmented[iSys]].set(volumeId, ilayer);
        (*decoder)["theta"].set(volumeId, 0);
        (*decoder)["phi"].set(volumeId, 0);
        // Get number of segmentation cells within the active volume
        std::array<uint, 3> numCells;
        // for vlume such as HCAL for which the numberOfCells function does not work
        // (because the G4 volume does not span the whole z of the detector)
        // use theta range passed via activeVolumesTheta option
        // if (m_fieldValuesSegmented[iSys] == m_hcalBarrelSysId)
        if (m_activeVolumesTheta[iSys].size() > 0) {
          int nPhiCells = phiThetaSegmentation->phiBins();
          double thetaCellSize = phiThetaSegmentation->gridSizeTheta();
          double thetaOffset = phiThetaSegmentation->offsetTheta();
          double thetaMin = m_activeVolumesTheta[iSys][ilayer];
          double thetaMax = TMath::Pi() - thetaMin;
          // debug
          debug() << "thetaMin, thetaMax = " << thetaMin << " " << thetaMax << endmsg;
          uint minThetaID = int(floor((thetaMin + 0.5 * thetaCellSize - thetaOffset) / thetaCellSize));
          uint maxThetaID = int(floor((thetaMax + 0.5 * thetaCellSize - thetaOffset) / thetaCellSize));
          uint nThetaCells = 1 + maxThetaID - minThetaID;
          numCells[0] = nPhiCells;
          numCells[1] = nThetaCells;
          numCells[2] = minThetaID;
        } else {
          numCells = det::utils::numberOfCells(volumeId, *phiThetaSegmentation);
        }
        // extrema 1: 0, ID of last cell in phi
        extrema[1] = std::make_pair(0, numCells[0] - 1);
        // extrema[2]: min theta ID, n theta cells
        extrema[2] = std::make_pair(numCells[2], numCells[1] + numCells[2] - 1);

        // for layer 0 of HCal barrel,  will be used for volume connecting
        if (m_fieldNamesSegmented[iSys] == "system" && m_fieldValuesSegmented[iSys] == m_hcalBarrelSysId &&
            ilayer == 0) {
          extremaHCalFirstLayerPhi = std::make_pair(0, numCells[0] - 1);
          extremaHCalFirstLayerTheta = std::make_pair(numCells[2], numCells[1] + numCells[2] - 1);
        }
        // for layer N-1 of ECal barrel,  will be used for volume connecting
        if (m_fieldNamesSegmented[iSys] == "system" && m_fieldValuesSegmented[iSys] == m_ecalBarrelSysId &&
            ilayer == (m_activeVolumesNumbersSegmented[iSys] - 1)) {
          eCalLastLayer = m_activeVolumesNumbersSegmented[iSys] - 1;
          extremaECalLastLayerPhi = std::make_pair(0, (numCells[0] - 1));
          extremaECalLastLayerTheta = std::make_pair(numCells[2], numCells[2] + (numCells[1] - 1));
        }

        debug() << "Layer: " << ilayer << endmsg;
        debug() << "Extrema[0]: " << extrema[0].first << " , " << extrema[0].second << endmsg;
        debug() << "Extrema[1]: " << extrema[1].first << " , " << extrema[1].second << endmsg;
        debug() << "Extrema[2]: " << extrema[2].first << " , " << extrema[2].second << endmsg;
        debug() << "Number of segmentation cells in phi, in theta, and min theta ID, : " << numCells << endmsg;
        // Loop over segmentation cells
        for (unsigned int iphi = 0; iphi < numCells[0]; iphi++) {
          for (unsigned int itheta = 0; itheta < numCells[1]; itheta++) {
            dd4hep::DDSegmentation::CellID cellId = volumeId;
            decoder->set(cellId, "phi", iphi);
            decoder->set(cellId, "theta",
                         itheta + numCells[2]); // start from the minimum existing theta cell in this layer
            uint64_t id = cellId;
            map.insert(std::pair<uint64_t, std::vector<uint64_t>>(
                id, det::utils::neighbours(*decoder, {m_activeFieldNamesSegmented[iSys], "phi", "theta"}, extrema, id,
                                           {false, true, false}, m_includeDiagonalCells)));
          }
        }
      }
    }

    // Loop over all cells in the HCal and retrieve existing cellIDs and find neighbours
    else if (segmentationType == "FCCSWHCalPhiTheta_k4geo") {
      // Loop over active layers
      for (unsigned int ilayer = 0; ilayer < m_activeVolumesNumbersSegmented[iSys]; ilayer++) {
        dd4hep::DDSegmentation::CellID volumeId = 0;
        // Get VolumeID
        (*decoder)[m_fieldNamesSegmented[iSys]].set(volumeId, m_fieldValuesSegmented[iSys]);
        (*decoder)[m_activeFieldNamesSegmented[iSys]].set(volumeId, ilayer);
        (*decoder)["theta"].set(volumeId, 0);
        (*decoder)["phi"].set(volumeId, 0);
        // Get number of segmentation cells within the active volume
        std::array<uint, 3> numCells;
        std::vector<int> thetaBins(hcalPhiThetaSegmentation->thetaBins(ilayer));
        numCells[0] = hcalPhiThetaSegmentation->phiBins();

        // Barrel
        if (m_fieldNamesSegmented[iSys] == "system" && m_fieldValuesSegmented[iSys] == m_hcalBarrelSysId) {
          // number of cells in the layer
          numCells[1] = thetaBins.size();
          // ID of first cell theta bin (smallest theta)
          if (numCells[1] > 0)
            numCells[2] = thetaBins[0];
          else
            error() << "Can not get number of cells from the segmentation object: " << segmentationType << endmsg;

          // for layer 0 of HCal barrel,  will be used for volume connecting
          if (ilayer == 0) {
            extremaHCalFirstLayerPhi = std::make_pair(0, numCells[0] - 1);
            extremaHCalFirstLayerTheta = std::make_pair(numCells[2], numCells[1] + numCells[2] - 1);
          }

          debug() << "Layer: " << ilayer << endmsg;
          debug() << "Number of segmentation cells in phi, in theta, and min theta ID, for HCal barrel : " << numCells
                  << endmsg;
          // Loop over segmentation cells
          for (unsigned int iphi = 0; iphi < numCells[0]; iphi++) {
            for (unsigned int itheta = 0; itheta < numCells[1]; itheta++) {
              dd4hep::DDSegmentation::CellID cellId = volumeId;
              decoder->set(cellId, "phi", iphi);
              decoder->set(cellId, "theta",
                           itheta + numCells[2]); // start from the minimum existing theta cell in this layer
              uint64_t id = cellId;
              map.insert(std::pair<uint64_t, std::vector<uint64_t>>(
                  id, hcalPhiThetaSegmentation->neighbours(id, m_includeDiagonalCellsHCal)));
            }
          }
        } // Barrel

        // For endcap, determine neighbours separately for positive- and negative-z part
        if (m_fieldNamesSegmented[iSys] == "system" && m_fieldValuesSegmented[iSys] == m_hcalEndcapSysId) {
          // number of cells in the layer
          numCells[1] =
              thetaBins.size() /
              2; // thetaBins contains theta IDs for both positive- and negative-z parts, hence divide size by 2

          // ID of first cell theta bin (smallest theta) in the positive-z part of the endcap
          if (numCells[1] > 0)
            numCells[2] = thetaBins[0];
          else
            error() << "Can not get number of cells from the segmentation object: " << segmentationType << endmsg;

          debug() << "Layer: " << ilayer << endmsg;
          debug() << "Number of segmentation cells in phi, in theta, and min theta ID, for positive-z HCal endcap : "
                  << numCells << endmsg;
          // Loop over segmentation cells
          for (unsigned int iphi = 0; iphi < numCells[0]; iphi++) {
            for (unsigned int itheta = 0; itheta < numCells[1]; itheta++) {
              dd4hep::DDSegmentation::CellID cellId = volumeId;
              decoder->set(cellId, "phi", iphi);
              decoder->set(cellId, "theta",
                           itheta + numCells[2]); // start from the minimum existing theta cell in this layer
              uint64_t id = cellId;
              map.insert(std::pair<uint64_t, std::vector<uint64_t>>(
                  id, hcalPhiThetaSegmentation->neighbours(id, m_includeDiagonalCellsHCal)));
            }
          }

          // ID of first cell theta bin (smallest theta) in the negative-z part of the endcap
          if (numCells[1] > 0)
            numCells[2] = thetaBins[thetaBins.size() / 2];

          debug() << "Layer: " << ilayer << endmsg;
          debug() << "Number of segmentation cells in phi, in theta, and min theta ID, for negative-z HCal endcap : "
                  << numCells << endmsg;
          // Loop over segmentation cells
          for (unsigned int iphi = 0; iphi < numCells[0]; iphi++) {
            for (unsigned int itheta = 0; itheta < numCells[1]; itheta++) {
              dd4hep::DDSegmentation::CellID cellId = volumeId;
              decoder->set(cellId, "phi", iphi);
              decoder->set(cellId, "theta",
                           itheta + numCells[2]); // start from the minimum existing theta cell in this layer
              uint64_t id = cellId;
              map.insert(std::pair<uint64_t, std::vector<uint64_t>>(
                  id, hcalPhiThetaSegmentation->neighbours(id, m_includeDiagonalCellsHCal)));
            }
          }
        } // Endcap
      }
    }

    // Loop over all cells in the HCal and retrieve existing cellIDs and find neighbours
    else if (segmentationType == "FCCSWHCalPhiRow_k4geo") {
      // Loop over active layers
      for (unsigned int ilayer = 0; ilayer < m_activeVolumesNumbersSegmented[iSys]; ilayer++) {
        dd4hep::DDSegmentation::CellID volumeId = 0;
        // Get VolumeID
        (*decoder)[m_fieldNamesSegmented[iSys]].set(volumeId, m_fieldValuesSegmented[iSys]);
        (*decoder)[m_activeFieldNamesSegmented[iSys]].set(volumeId, ilayer);
        (*decoder)["row"].set(volumeId, 0);
        (*decoder)["phi"].set(volumeId, 0);
        // Get number of segmentation cells within the active volume
        std::array<int, 3> numCells;
        std::vector<int> cellIndexes(hcalPhiRowSegmentation->cellIndexes(ilayer));
        numCells[0] = hcalPhiRowSegmentation->phiBins();

        // Barrel
        if (m_fieldNamesSegmented[iSys] == "system" && m_fieldValuesSegmented[iSys] == m_hcalBarrelSysId) {
          // number of cells in the layer
          numCells[1] = cellIndexes.size();
          // minimum cell index
          if (numCells[1] > 0)
            numCells[2] = cellIndexes.front();
          else
            error() << "Can not get number of cells from the segmentation object: " << segmentationType << endmsg;

          // for layer 0 of HCal barrel,  will be used for volume connecting
          if (ilayer == 0) {
            extremaHCalFirstLayerPhi = std::make_pair(0, numCells[0] - 1);
            extremaHCalFirstLayerTheta = std::make_pair(cellIndexes.front(), cellIndexes.back());
          }

          debug() << "Layer: " << ilayer << endmsg;
          debug() << "Number of segmentation cells in phi, in row, and min cell index, for HCal barrel : " << numCells
                  << endmsg;
          // Loop over segmentation cells
          for (int iphi = 0; iphi < numCells[0]; iphi++) {
            for (int irow = 0; irow < numCells[1]; irow++) {
              dd4hep::DDSegmentation::CellID cellId = volumeId;
              decoder->set(cellId, "phi", iphi);
              decoder->set(cellId, "row",
                           irow + numCells[2]); // start from the minimum existing cell index in this layer
              uint64_t id = cellId;
              map.insert(std::pair<uint64_t, std::vector<uint64_t>>(id, hcalPhiRowSegmentation->neighbours(id)));
            }
          }
        } // Barrel

        // For endcap, determine neighbours separately for positive- and negative-z part
        if (m_fieldNamesSegmented[iSys] == "system" && m_fieldValuesSegmented[iSys] == m_hcalEndcapSysId) {
          // number of cells in the layer
          numCells[1] =
              cellIndexes.size() /
              2; // cellIndexes contains cell IDs for both positive- and negative-z parts, hence divide size by 2

          // minimum cell index in the positive-z part of the endcap
          if (numCells[1] > 0)
            numCells[2] = cellIndexes.front();
          else
            error() << "Can not get number of cells from the segmentation object: " << segmentationType << endmsg;

          debug() << "Layer: " << ilayer << endmsg;
          debug() << "Number of segmentation cells in phi, in row, and min cell index, for positive-z HCal endcap : "
                  << numCells << endmsg;
          // Loop over segmentation cells
          for (int iphi = 0; iphi < numCells[0]; iphi++) {
            for (int irow = 0; irow < numCells[1]; irow++) {
              dd4hep::DDSegmentation::CellID cellId = volumeId;
              decoder->set(cellId, "phi", iphi);
              decoder->set(cellId, "row",
                           irow + numCells[2]); // start from the minimum existing cell index in this layer
              unsigned int pseudoLayer = hcalPhiRowSegmentation->definePseudoLayer(cellId);
              decoder->set(cellId, "pseudoLayer", pseudoLayer);
              uint64_t id = cellId;
              map.insert(std::pair<uint64_t, std::vector<uint64_t>>(id, hcalPhiRowSegmentation->neighbours(id)));
            }
          }

          // minimum cell index in the negative-z part of the endcap
          if (numCells[1] > 0)
            numCells[2] = cellIndexes.back();

          debug() << "Layer: " << ilayer << endmsg;
          debug() << "Number of segmentation cells in phi, in row, and min cell index, for negative-z HCal endcap : "
                  << numCells << endmsg;
          // Loop over segmentation cells
          for (int iphi = 0; iphi < numCells[0]; iphi++) {
            for (int irow = 0; irow < numCells[1]; irow++) {
              dd4hep::DDSegmentation::CellID cellId = volumeId;
              decoder->set(cellId, "phi", iphi);
              decoder->set(cellId, "row", irow + numCells[2]); // start from the minimum cell index in this layer
              unsigned int pseudoLayer = hcalPhiRowSegmentation->definePseudoLayer(cellId);
              decoder->set(cellId, "pseudoLayer", pseudoLayer);
              uint64_t id = cellId;
              map.insert(std::pair<uint64_t, std::vector<uint64_t>>(id, hcalPhiRowSegmentation->neighbours(id)));
            }
          }
        } // Endcap
      }
    }

    else if (segmentationType == "FCCSWGridModuleThetaMerged_k4geo") {
      std::vector<std::pair<int, int>> extrema;
      // extrema[0]: min layer, n layers
      extrema.push_back(std::make_pair(0, m_activeVolumesNumbersSegmented[iSys] - 1));
      extrema.push_back(std::make_pair(0, 0));
      extrema.push_back(std::make_pair(0, 0));
      for (unsigned int ilayer = 0; ilayer < m_activeVolumesNumbersSegmented[iSys]; ilayer++) {
        dd4hep::DDSegmentation::CellID volumeId = 0;
        // Get VolumeID
        (*decoder)[m_fieldNamesSegmented[iSys]].set(volumeId, m_fieldValuesSegmented[iSys]);
        (*decoder)[m_activeFieldNamesSegmented[iSys]].set(volumeId, ilayer);
        (*decoder)["theta"].set(volumeId, 0);
        (*decoder)["module"].set(volumeId, 0);
        // Get number of segmentation cells within the active volume
        // numberOfCells: return Array of the number of cells in (module, theta) and the minimum theta ID.
        auto numCells = det::utils::numberOfCells(volumeId, *moduleThetaSegmentation);
        // extrema 1: min module number (0), max module number
        extrema[1] = std::make_pair(0, (numCells[0] - 1) * moduleThetaSegmentation->mergedModules(ilayer));
        // extrema[2]: min theta ID, n (merged) theta cells
        extrema[2] = std::make_pair(numCells[2], numCells[2] + (numCells[1] - 1) *
                                                                   moduleThetaSegmentation->mergedThetaCells(ilayer));
        // for layer N-1 of ECal barrel,  will be used for volume connecting
        if (ilayer == (m_activeVolumesNumbersSegmented[iSys] - 1) && m_fieldNamesSegmented[iSys] == "system" &&
            m_fieldValuesSegmented[iSys] == m_ecalBarrelSysId) {
          eCalLastLayer = m_activeVolumesNumbersSegmented[iSys] - 1;
          extremaECalLastLayerModule =
              std::make_pair(0, (numCells[0] - 1) * moduleThetaSegmentation->mergedModules(eCalLastLayer));
          extremaECalLastLayerTheta = std::make_pair(
              numCells[2], numCells[2] + (numCells[1] - 1) * moduleThetaSegmentation->mergedThetaCells(eCalLastLayer));
        }
        debug() << "Layer: " << ilayer << endmsg;
        debug() << "Extrema[0]: " << extrema[0].first << " , " << extrema[0].second << endmsg;
        debug() << "Extrema[1]: " << extrema[1].first << " , " << extrema[1].second << endmsg;
        debug() << "Extrema[2]: " << extrema[2].first << " , " << extrema[2].second << endmsg;
        debug() << "Number of segmentation cells in (module,theta): " << numCells << endmsg;
        // Loop over segmentation cells to find neighbours in ECAL
        for (int imodule = extrema[1].first; imodule <= extrema[1].second;
             imodule += moduleThetaSegmentation->mergedModules(ilayer)) {
          for (int itheta = extrema[2].first; itheta <= extrema[2].second;
               itheta += moduleThetaSegmentation->mergedThetaCells(ilayer)) {
            dd4hep::DDSegmentation::CellID cellId = volumeId;
            decoder->set(cellId, "module", imodule);
            decoder->set(cellId, "theta", itheta); // start from the minimum existing theta cell in this layer
            uint64_t id = cellId;
            auto neighboursList = det::utils::neighbours_ModuleThetaMerged(
                *moduleThetaSegmentation, *decoder, {m_activeFieldNamesSegmented[iSys], "module", "theta"}, extrema, id, m_includeDiagonalCells);

            // check if we are on an edge in theta, and add neighbours in the
            // endcap calorimeter if so

            if (((itheta == extrema[2].first) || (itheta == extrema[2].second)) && m_connectECal) {
              // find barrel cell position
              double eCalBarrelTheta = m_EMB_theta_lookup[ilayer];
              if (itheta == extrema[2].second)
                eCalBarrelTheta = TMath::Pi() - eCalBarrelTheta;
              double eCalBarrelPhi = m_EMB_phi_lookup[ilayer][imodule];

              for (uint iSys2 = 0; iSys2 < m_readoutNamesSegmented.size(); iSys2++) {
                // get segmentation
                dd4hep::DDSegmentation::Segmentation* aSegmentation2 =
                    m_geoSvc->getDetector()->readout(m_readoutNamesSegmented[iSys2]).segmentation().segmentation();
                if (aSegmentation2 == nullptr) {
                  error() << "Segmentation does not exist." << endmsg;
                  return StatusCode::FAILURE;
                }
	      
                std::string segmentationType2 = aSegmentation2->type();
                if (segmentationType2 == "FCCSWEndcapTurbine_k4geo") {
                  auto endcapDecoder =
                      m_geoSvc->getDetector()->readout(m_readoutNamesSegmented[iSys2]).idSpec().decoder();
                  if (ecalEndcapTurbineSegmentation == nullptr) {
                    ecalEndcapTurbineSegmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWEndcapTurbine_k4geo*>(aSegmentation2);
                  }
 
                  int iWheel = 2;
                  int iSide;
                  if (itheta == extrema[2].first) {
                    iSide = 1;
                  } else {
                    iSide = -1;
                  }

                  unsigned numECCellsRho = ecalEndcapTurbineSegmentation->numCellsRho(iWheel);
                  // Loop over segmentation cells
                  dd4hep::DDSegmentation::CellID endcapCellId = 0;
		
                  (*endcapDecoder)[m_fieldNamesSegmented[iSys2]].set(endcapCellId, m_ecalEndcapSysId);
                  (*endcapDecoder)["z"].set(endcapCellId, 0);
                  (*endcapDecoder)["wheel"].set(endcapCellId, iWheel);
                  (*endcapDecoder)["side"].set(endcapCellId, iSide);

                  unsigned iMatchRho = 999;
                  float minDelTheta = 1.0;

                  for (unsigned iECrho = 0; iECrho < numECCellsRho; iECrho++) {
                    double endcapTheta = m_EMEC_theta_lookup[iECrho];
                    if (itheta == extrema[2].second)
                      endcapTheta = TMath::Pi() - endcapTheta;
                    float delTheta = TMath::Abs(endcapTheta - eCalBarrelTheta);
                    if (delTheta < minDelTheta) {
                      minDelTheta = delTheta;
                      iMatchRho = iECrho;
                    }
                  }
                  if (minDelTheta > 0.5)
                    continue;
                  unsigned iECmodule;
                  if (iSide == -1) {
                    iECmodule = m_EMEC_h_module_vs_phi_neg[iMatchRho]->GetBinContent(
                        m_EMEC_h_module_vs_phi_neg[iMatchRho]->FindBin(eCalBarrelPhi));
                  } else {
                    iECmodule = m_EMEC_h_module_vs_phi_pos[iMatchRho]->GetBinContent(
                        m_EMEC_h_module_vs_phi_pos[iMatchRho]->FindBin(eCalBarrelPhi));
                  }
                  debug() << "Adding endcap cell to barrel list with module " << iECmodule << " for side " << iSide
                          << " and barrel phi " << eCalBarrelPhi << endmsg;
                  (*endcapDecoder)["module"].set(endcapCellId, iECmodule);
                  (*endcapDecoder)["rho"].set(endcapCellId, iMatchRho);
                  (*endcapDecoder)[m_activeFieldNamesSegmented[iSys]].set(endcapCellId,
                                                ecalEndcapTurbineSegmentation->expLayer(iWheel, iMatchRho, 0));
                  neighboursList.push_back(endcapCellId);
                }
              }
            }
            map.insert(std::pair<uint64_t, std::vector<uint64_t>>(id, neighboursList));
          }
        }
      }
    } else if (segmentationType == "FCCSWEndcapTurbine_k4geo") {
      // Loop over active layers
      std::vector<std::pair<int, int>> extrema;
      extrema.push_back(std::make_pair(0, 0)); // modules (set per wheel)
      extrema.push_back(std::make_pair(0, 0)); // rho (set per wheel)
      extrema.push_back(std::make_pair(0, 0)); // z (set per wheel)

      unsigned layerOffset[3];
      layerOffset[0] = 0;
      layerOffset[1] =
          ecalEndcapTurbineSegmentation->numCellsRhoCalib(0) * ecalEndcapTurbineSegmentation->numCellsZCalib(0);
      layerOffset[2] = layerOffset[1] + ecalEndcapTurbineSegmentation->numCellsRhoCalib(1) *
                                            ecalEndcapTurbineSegmentation->numCellsZCalib(1);
      for (int iSide = -1; iSide < 2; iSide += 2) {
        for (unsigned int iWheel = 0; iWheel < 3; iWheel++) {
          dd4hep::DDSegmentation::CellID volumeId = 0;
          // Get VolumeID
          // for ECAL OK (volume extends along z)
          (*decoder)[m_fieldNamesSegmented[iSys]].set(volumeId, m_fieldValuesSegmented[iSys]);
          (*decoder)["side"].set(volumeId, iSide);
          (*decoder)["module"].set(volumeId, 0);
          (*decoder)["rho"].set(volumeId, 0);
          (*decoder)["z"].set(volumeId, 0);

          unsigned numModules = ecalEndcapTurbineSegmentation->nModules(iWheel);
          unsigned numCellsRho = ecalEndcapTurbineSegmentation->numCellsRho(iWheel);
          unsigned numCellsZ = ecalEndcapTurbineSegmentation->numCellsZ(iWheel);
          unsigned numCellsRhoCalib = ecalEndcapTurbineSegmentation->numCellsRhoCalib(iWheel);
          unsigned numCellsZCalib = ecalEndcapTurbineSegmentation->numCellsZCalib(iWheel);
          // extrema 0: 0, ID of last module
          extrema[0] = std::make_pair(0, numModules - 1);
          // extrema[1]: 0, ID of last rho cell
          extrema[1] = std::make_pair(0, numCellsRho - 1);
          // extrema[2]: 0, ID of last z cell
          extrema[2] = std::make_pair(0, numCellsZ - 1);

          debug() << "Wheel: " << iWheel << endmsg;
          debug() << "Extrema[0]: " << extrema[0].first << " , " << extrema[0].second << endmsg;
          debug() << "Extrema[1]: " << extrema[1].first << " , " << extrema[1].second << endmsg;
          debug() << "Extrema[2]: " << extrema[2].first << " , " << extrema[2].second << endmsg;

          // Set up some variables that will be useful when finding neighbors
          // in the EM barrel.  These will be calculated the first time through
          // the loop
          float minRhoEMBarrel = -1.;

          // Loop over segmentation cells
          for (unsigned imodule = 0; imodule < numModules; imodule++) {
            for (unsigned irho = 0; irho < numCellsRho; irho++) {
              for (unsigned iz = 0; iz < numCellsZ; iz++) {
                // check if we're at the boundary between wheels
                bool atInnerBoundary = false, atOuterBoundary = false, nextToEMBarrel = false;
                if (iWheel > 0 && irho == 0)
                  atInnerBoundary = true;
                if (iWheel < 2 && irho == numCellsRho - 1)
                  atOuterBoundary = true;
                dd4hep::DDSegmentation::CellID cellId = volumeId;
                decoder->set(cellId, "wheel", iWheel);
                decoder->set(cellId, "module", imodule);
                decoder->set(cellId, "rho", irho);
                decoder->set(cellId, "z", iz);

                double endcapRho = ecalEndcapTurbineSegmentation->rho(cellId);
                if (iWheel == 2 && iz == 0 && ((endcapRho > minRhoEMBarrel) || (minRhoEMBarrel < 0.))) {
                  nextToEMBarrel = true;
                }

                double endcapPhi = TMath::ATan2(ecalEndcapTurbineSegmentation->position(cellId).y(),
                                                ecalEndcapTurbineSegmentation->position(cellId).x());

                unsigned iLayerZ = iz / (numCellsZ / numCellsZCalib);
                unsigned iLayerRho = irho / (numCellsRho / numCellsRhoCalib);
                unsigned iLayer = layerOffset[iWheel] + iLayerRho * numCellsZCalib + iLayerZ;
                decoder->set(cellId, m_activeFieldNamesSegmented[iSys], iLayer);
                uint64_t id = cellId;
                debug() << "Mapping cell " << cellId << " " << std::hex << cellId << std::dec << endmsg;
                auto neighboursList = det::utils::neighbours(*decoder, {"module", "rho", "z"}, extrema, id,
                                                             {true, false, false}, m_includeDiagonalCells);
                // now correct the layer index for the neighbours, since the
                // rho and z indices change
                unsigned idx = 0;
                for (auto nCell : neighboursList) {
                  unsigned lirho = decoder->get(nCell, "rho");
                  unsigned liz = decoder->get(nCell, "z");
                  unsigned correctLayer = layerOffset[iWheel] +
                                          lirho / (numCellsRho / numCellsRhoCalib) * numCellsZCalib +
                                          liz / (numCellsZ / numCellsZCalib);
                  decoder->set(nCell, m_activeFieldNamesSegmented[iSys], correctLayer);
                  neighboursList[idx] = nCell;
                  idx++;
                }
                // if we're at a boundary between wheels, add the
                // appropriate cells in the neighboring wheel
                if (atInnerBoundary || atOuterBoundary) {
                  unsigned otherWheel, otherRho, newLayerOffset;
                  if (atInnerBoundary) {
                    otherWheel = iWheel - 1;
                    otherRho = ecalEndcapTurbineSegmentation->numCellsRho(otherWheel);
                  } else {
                    otherWheel = iWheel + 1;
                    otherRho = 0;
                  }
                  newLayerOffset = layerOffset[otherWheel];
                  unsigned numModulesotherWheel = ecalEndcapTurbineSegmentation->nModules(otherWheel);
                  unsigned numCellsZotherWheel = ecalEndcapTurbineSegmentation->numCellsZ(otherWheel);
                  unsigned numCellsRhootherWheel = ecalEndcapTurbineSegmentation->numCellsRho(otherWheel);
                  unsigned numCellsZCalibotherWheel = ecalEndcapTurbineSegmentation->numCellsZCalib(otherWheel);
                  unsigned numCellsRhoCalibotherWheel = ecalEndcapTurbineSegmentation->numCellsRhoCalib(otherWheel);
                  uint64_t newZ;
                  int newModule;
                  if (numCellsZ > numCellsZotherWheel) {
                    newZ = decoder->get(cellId, "z") / ((1.0 * numCellsZ) / numCellsZotherWheel);
                  } else {
                    newZ = decoder->get(cellId, "z") * ((1.0 * numCellsZotherWheel) / numCellsZ);
                  }
                  // calculate offset in module index due to differences in
                  // blade angle in different wheels
                  double rho, z;
                  if (atInnerBoundary) {
                    rho = ecalEndcapTurbineSegmentation->offsetRho(iWheel);
                  } else {
                    rho = ecalEndcapTurbineSegmentation->offsetRho(otherWheel);
                  }
                  z = ((int)(iz - numCellsZ / 2)) * ecalEndcapTurbineSegmentation->gridSizeZ(iWheel);
                  double alpha = ecalEndcapTurbineSegmentation->bladeAngle(iWheel);
                  double alphaotherWheel = ecalEndcapTurbineSegmentation->bladeAngle(otherWheel);
                  int moduleOffset = numModulesotherWheel * (z / (2 * TMath::Pi() * rho)) *
                                     (1. / TMath::Tan(alpha) - 1. / TMath::Tan(alphaotherWheel));

                  if (numModules > numModulesotherWheel) {
                    newModule =
                        decoder->get(cellId, "module") / ((1.0 * numModules) / numModulesotherWheel) + moduleOffset;
                  } else {
                    newModule =
                        decoder->get(cellId, "module") * ((1.0 * numModulesotherWheel) / numModules) + moduleOffset;
                  }
                  if (newModule < 0)
                    newModule = numModulesotherWheel + newModule;
                  unsigned newLayer =
                      newLayerOffset +
                      otherRho / (numCellsRhootherWheel / numCellsRhoCalibotherWheel) * numCellsZCalibotherWheel +
                      newZ / (numCellsZotherWheel / numCellsZCalibotherWheel);
                  for (int ibmod = 0; ibmod < 2; ibmod++) {
                    for (int izmod = 0; izmod < 2; izmod++) {
                      uint64_t newCellId = cellId;
                      decoder->set(newCellId, "wheel", otherWheel);
                      decoder->set(newCellId, "rho", otherRho);
                      decoder->set(newCellId, "z", newZ + izmod / 1);
                      decoder->set(newCellId, "module", newModule + ibmod / 1);
                      decoder->set(newCellId, m_activeFieldNamesSegmented[iSys], newLayer);
                      neighboursList.push_back(newCellId);
                    }
                  }
                }
                // add neighboring cells in the EM barrel calorimeter, if relevant
                if (nextToEMBarrel && m_connectECal) {

                  double endcapTheta = m_EMEC_theta_lookup[irho];
                  if (iSide == -1)
                    endcapTheta = TMath::Pi() - endcapTheta;

                  for (uint iSys2 = 0; iSys2 < m_readoutNamesSegmented.size(); iSys2++) {
                    // get segmentation
                    dd4hep::DDSegmentation::Segmentation* aSegmentation2 =
                        m_geoSvc->getDetector()->readout(m_readoutNamesSegmented[iSys2]).segmentation().segmentation();
                    if (aSegmentation2 == nullptr) {
                      error() << "Segmentation does not exist." << endmsg;
                      return StatusCode::FAILURE;
                    }

                    std::string segmentationType2 = aSegmentation2->type();
                    if (segmentationType2 == "FCCSWGridModuleThetaMerged_k4geo") {
                      moduleThetaSegmentation =
                          dynamic_cast<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo*>(aSegmentation2);
                      auto barrelDecoder =
                          m_geoSvc->getDetector()->readout(m_readoutNamesSegmented[iSys2]).idSpec().decoder();
                      m_geoSvc->getDetector()->readouts().find(m_readoutNamesSegmented[iSys2]);
                      // Loop over all relevant cells in the barrel and find position
                      unsigned iBarrelTheta = 0;
                      unsigned nLayersBarrel = m_activeVolumesNumbersSegmented[iSys2];
                      auto numCells = det::utils::numberOfCells(volumeId, *moduleThetaSegmentation);
                      unsigned iMatchLayer = 0;
                      float minDelTheta = 1.0;

                      for (unsigned iBarrelLayer = 0; iBarrelLayer < nLayersBarrel; iBarrelLayer++) {
                        double eCalBarrelTheta = m_EMB_theta_lookup[iBarrelLayer];
                        if (iSide == -1)
                          eCalBarrelTheta = TMath::Pi() - eCalBarrelTheta;
                        float delTheta = TMath::Abs(endcapTheta - eCalBarrelTheta);
                        if (delTheta < minDelTheta) {
                          minDelTheta = delTheta;
                          iMatchLayer = iBarrelLayer;
                        }
                      }
                      if (minDelTheta > 0.5)
                        continue;

                      unsigned iBarrelModule = m_EMB_h_module_vs_phi[iMatchLayer]->GetBinContent(
                          m_EMB_h_module_vs_phi[iMatchLayer]->FindBin(endcapPhi));
                      // set volumeID correctly
                      dd4hep::DDSegmentation::CellID barrelVolumeId = 0;

                      (*barrelDecoder)[m_fieldNamesSegmented[iSys2]].set(barrelVolumeId, m_ecalBarrelSysId);
                      (*barrelDecoder)["module"].set(barrelVolumeId, iBarrelModule);
                      (*barrelDecoder)[m_activeFieldNamesSegmented[iSys2]].set(barrelVolumeId, iMatchLayer);

                      numCells = det::utils::numberOfCells(barrelVolumeId, *moduleThetaSegmentation);

                      if (iSide > 0) {
                        iBarrelTheta = numCells[2];
                      } else {
                        iBarrelTheta =
                            numCells[2] + (numCells[1] - 1) * moduleThetaSegmentation->mergedThetaCells(iMatchLayer);
                      }
                      (*barrelDecoder)["theta"].set(barrelVolumeId, iBarrelTheta);
                      debug() << "We have found a neighbor for cellId  " << cellId << " with module " << iBarrelModule
                              << " from phi " << endcapPhi << " and bin "
                              << m_EMB_h_module_vs_phi[iMatchLayer]->FindBin(endcapPhi) << "!" << endmsg;
                      neighboursList.push_back(barrelVolumeId);
                      break;
                    }
                  }
                }
                map.insert(std::pair<uint64_t, std::vector<uint64_t>>(id, neighboursList));
              }
            }
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
    info() << "total number of cells in map so far:  " << map.size() << endmsg;
  }

  //////////////////////////////////////////////////
  ///   connection HCAL Barrel + HCAL Endcap     ///
  //////////////////////////////////////////////////
  // FCCSWHCalPhiTheta segmentation for Barrel while FCCSWHCalPhiTheta/FCCSWHCalPhiRow for Endcap
  if (m_connectHCal && hcalBarrelSegmentation && (hcalEndcapSegmentation || hcalEndcapPhiRowSegmentation)) {
    info() << "Linking the HCal Barrel and HCal Endcap" << endmsg;

    // Get the part1 and part2 endcap layer indexes that will be used to connect Barrel cells to the Endcap cells
    std::vector<int> minLayerIdEndcap = {0, 0};
    std::vector<int> maxLayerIdEndcap = {0, 0};
    std::vector<std::pair<uint, uint>> minMaxLayerId((hcalEndcapSegmentation) ? hcalEndcapSegmentation->getMinMaxLayerId() :
                                                                                hcalEndcapPhiRowSegmentation->getMinMaxLayerId());
    if (minMaxLayerId.empty()) {
      error() << "hcalEndcap segmentation is not configured correctly, cannot link Endcap and Barrel!" << endmsg;
      return StatusCode::FAILURE;
    }

    // Loop over active layers
    for (unsigned int ilayer = 0; ilayer < m_activeVolumesNumbersSegmented[hcalBarrelId]; ilayer++) {
      debug() << "HCal Barrel Layer: " << ilayer << endmsg;
      std::array<dd4hep::DDSegmentation::CellID, 2> barrelCellID = {0, 0};
      (*decoderHCalBarrel)[m_fieldNamesSegmented[hcalBarrelId]].set(barrelCellID[0],
                                                                    m_fieldValuesSegmented[hcalBarrelId]);
      (*decoderHCalBarrel)[m_activeFieldNamesSegmented[hcalBarrelId]].set(barrelCellID[0], ilayer);
      (*decoderHCalBarrel)["theta"].set(barrelCellID[0], 0);
      (*decoderHCalBarrel)["phi"].set(barrelCellID[0], 0);
      (*decoderHCalBarrel)[m_fieldNamesSegmented[hcalBarrelId]].set(barrelCellID[1],
                                                                    m_fieldValuesSegmented[hcalBarrelId]);
      (*decoderHCalBarrel)[m_activeFieldNamesSegmented[hcalBarrelId]].set(barrelCellID[1], ilayer);
      (*decoderHCalBarrel)["theta"].set(barrelCellID[1], 0);
      (*decoderHCalBarrel)["phi"].set(barrelCellID[1], 0);

      std::vector<int> thetaBins(hcalBarrelSegmentation->thetaBins(ilayer));
      int minCellId = thetaBins.front(); // first cell (lowest theta bin) in the barrel layer
      int maxCellId = thetaBins.back();  // last cell (highest theta bin) in the barrel layer
      int nPhiBins = hcalBarrelSegmentation->phiBins();
      for (int iphi = 0; iphi < nPhiBins; iphi++) {
        // set the final bits for the first cell in the barrel layer
        decoderHCalBarrel->set(barrelCellID[0], "phi", iphi);
        decoderHCalBarrel->set(barrelCellID[0], "theta", minCellId);
        // set the final bits for the last cell in the barrel layer
        decoderHCalBarrel->set(barrelCellID[1], "phi", iphi);
        decoderHCalBarrel->set(barrelCellID[1], "theta", maxCellId);

        // Number of neighbours for the first and the last cell in the barrel layer
        std::array<uint64_t, 2> numberOfNeighbours = {map[barrelCellID[0]].size(), map[barrelCellID[1]].size()};

        // min and max polar angle for the first cell in the barrel layer
        std::array<double, 2> minCellIdTheta = hcalBarrelSegmentation->cellTheta(barrelCellID[0]);
        // min and max polar angle for the last cell in the barrel layer
        std::array<double, 2> maxCellIdTheta = hcalBarrelSegmentation->cellTheta(barrelCellID[1]);

        // part1 min and max layer index
        minLayerIdEndcap[0] = minMaxLayerId[0].first;
        maxLayerIdEndcap[0] = minMaxLayerId[0].second;
        // part2 min and max layer index
        minLayerIdEndcap[1] = minMaxLayerId[1].first;
        maxLayerIdEndcap[1] = minMaxLayerId[1].second;

        // create the endcap cellID that should become the neighbour for the barrel cell
        dd4hep::DDSegmentation::CellID endcapCellID = 0;
        (*decoderHCalEndcap)[m_fieldNamesSegmented[hcalEndcapId]].set(endcapCellID,
                                                                      m_fieldValuesSegmented[hcalEndcapId]);
        (*decoderHCalEndcap)["phi"].set(endcapCellID, iphi);

        // loop over the part1 layers
        for (int part1LayerId = minLayerIdEndcap[0]; part1LayerId <= maxLayerIdEndcap[0]; part1LayerId++) {
          // set the layer bit field for part1 layer
          (*decoderHCalEndcap)[m_activeFieldNamesSegmented[hcalEndcapId]].set(endcapCellID, part1LayerId);

          // get cell indexes in the given endcap layer
          std::vector<int> endcapCellIndexes( (hcalEndcapSegmentation) ? hcalEndcapSegmentation->thetaBins(part1LayerId) :
                                                                         hcalEndcapPhiRowSegmentation->cellIndexes(part1LayerId) );

          // consider all cells in the first layer of part1 endcap
          if (part1LayerId == minLayerIdEndcap[0]) {
            for (auto idx : endcapCellIndexes) {
              if(hcalEndcapSegmentation) // phi-theta segmentation
              {
                decoderHCalEndcap->set(endcapCellID, "theta", idx);
              }
              else // phi-row segmentation
              {
                decoderHCalEndcap->set(endcapCellID, "row", idx);
                // set the pseudo-layer bit field if the phi-row segmentation is used for the Endcap
                unsigned int pseudoLayer = hcalEndcapPhiRowSegmentation->definePseudoLayer(endcapCellID);
                decoderHCalEndcap->set(endcapCellID, "pseudoLayer", pseudoLayer);
              }

              // get thetaMin and thetaMax of the endcap cell
              std::array<double, 2> cTheta = (hcalEndcapSegmentation) ? hcalEndcapSegmentation->cellTheta(endcapCellID) :
                                                                        hcalEndcapPhiRowSegmentation->cellTheta(endcapCellID);
              // add in the neighbours list if it overlaps with first barrel cell in theta
              if ((cTheta[0] <= minCellIdTheta[0] && cTheta[1] > minCellIdTheta[0]) ||
                  (cTheta[0] > minCellIdTheta[0] && cTheta[0] < minCellIdTheta[1])) {
                map[barrelCellID[0]].push_back(
                    endcapCellID); // add to the neighbours list of the first cell in the barrel layer
                map[endcapCellID].push_back(
                    barrelCellID[0]); // add the first barrel cell to the list of the endcap cell neighbours
                // add the next and previous phi modules as well
                // previous: if the current is 0 then previous is the last bin (id = phiBins - 1) else current - 1
                (*decoderHCalEndcap)["phi"].set(endcapCellID, (iphi == 0) ? (nPhiBins - 1) : (iphi - 1));
                map[barrelCellID[0]].push_back(
                    endcapCellID); // add to the neighbours list of the first cell in the barrel layer
                map[endcapCellID].push_back(
                    barrelCellID[0]); // add the first barrel cell to the list of the endcap cell neighbours
                // next: if the current is the last bin (id = phiBins - 1) then the next is the first bin (id = 0) else
                // current + 1
                (*decoderHCalEndcap)["phi"].set(endcapCellID, (iphi == (nPhiBins - 1)) ? 0 : (iphi + 1));
                map[barrelCellID[0]].push_back(
                    endcapCellID); // add to the neighbours list of the first cell in the barrel layer
                map[endcapCellID].push_back(
                    barrelCellID[0]); // add the first barrel cell to the list of the endcap cell neighbours
                // restore back the iphi value
                (*decoderHCalEndcap)["phi"].set(endcapCellID, iphi);
              }
              // add in the neighbours list if it overlaps with the last barrel cell in theta
              if ((cTheta[0] <= maxCellIdTheta[0] && cTheta[1] > maxCellIdTheta[0]) ||
                  (cTheta[0] > maxCellIdTheta[0] && cTheta[0] < maxCellIdTheta[1])) {
                map[barrelCellID[1]].push_back(
                    endcapCellID); // add to the neighbours list of the last cell in the barrel layer
                map[endcapCellID].push_back(
                    barrelCellID[1]); // add the last barrel cell to the list of the endcap cell neighbours
                // add the next and previous phi modules as well
                // previous: if the current is 0 then previous is the last bin (id = phiBins - 1) else current - 1
                (*decoderHCalEndcap)["phi"].set(endcapCellID, (iphi == 0) ? (nPhiBins - 1) : (iphi - 1));
                map[barrelCellID[1]].push_back(
                    endcapCellID); // add to the neighbours list of the last cell in the barrel layer
                map[endcapCellID].push_back(
                    barrelCellID[1]); // add the last barrel cell to the list of the endcap cell neighbours
                // next: if the current is the last bin (id = phiBins - 1) then the next is the first bin (id = 0) else
                // current + 1
                (*decoderHCalEndcap)["phi"].set(endcapCellID, (iphi == (nPhiBins - 1)) ? 0 : (iphi + 1));
                map[barrelCellID[1]].push_back(
                    endcapCellID); // add to the neighbours list of the last cell in the barrel layer
                map[endcapCellID].push_back(
                    barrelCellID[1]); // add the last barrel cell to the list of the endcap cell neighbours
                // restore back the iphi value
                (*decoderHCalEndcap)["phi"].set(endcapCellID, iphi);
              }
            }
          } else // consider only the first cell
          {
            // consider the the first cell (highest theta) from the positive-z part endcap layer
            // to be added in the neighbours list of the first barrel cell (lowest theta bin)
            if(hcalEndcapSegmentation) // phi-theta segmentation
            {
              decoderHCalEndcap->set(endcapCellID, "theta",
                endcapCellIndexes[endcapCellIndexes.size() / 2 - 1]);
            }
            else // phi-row segmentation
            {
              decoderHCalEndcap->set(endcapCellID, "row",
                endcapCellIndexes[0]);
              // set the pseudo-layer bit field if the phi-row segmentation is used for the Endcap
              unsigned int pseudoLayer = hcalEndcapPhiRowSegmentation->definePseudoLayer(endcapCellID);
              decoderHCalEndcap->set(endcapCellID, "pseudoLayer", pseudoLayer);
            }

            // get thetaMin and thetaMax of the endcap cell
            std::array<double, 2> cTheta = (hcalEndcapSegmentation) ? hcalEndcapSegmentation->cellTheta(endcapCellID) :
                                                                      hcalEndcapPhiRowSegmentation->cellTheta(endcapCellID);
            // add in the neighbours list if it overlaps in theta with the first barrel cell (lowest theta bin)
            if ((cTheta[0] <= minCellIdTheta[0] && cTheta[1] > minCellIdTheta[0]) ||
                (cTheta[0] > minCellIdTheta[0] && cTheta[0] < minCellIdTheta[1])) {
              map[barrelCellID[0]].push_back(
                  endcapCellID); // add to the neighbours list of the first (lowest theta bin) cell in the barrel layer
              map[endcapCellID].push_back(
                  barrelCellID[0]); // add the first barrel cell to the list of the endcap cell neighbours
              // add the next and previous phi modules as well
              // previous: if the current is 0 then previous is the last bin (id = phiBins - 1) else current - 1
              (*decoderHCalEndcap)["phi"].set(endcapCellID, (iphi == 0) ? (nPhiBins - 1) : (iphi - 1));
              map[barrelCellID[0]].push_back(
                  endcapCellID); // add to the neighbours list of the first cell in the barrel layer
              map[endcapCellID].push_back(
                  barrelCellID[0]); // add the last barrel cell to the list of the endcap cell neighbours
              // next: if the current is the last bin (id = phiBins - 1) then the next is the first bin (id = 0) else
              // current + 1
              (*decoderHCalEndcap)["phi"].set(endcapCellID, (iphi == (nPhiBins - 1)) ? 0 : (iphi + 1));
              map[barrelCellID[0]].push_back(
                  endcapCellID); // add to the neighbours list of the first cell in the barrel layer
              map[endcapCellID].push_back(
                  barrelCellID[0]); // add the last barrel cell to the list of the endcap cell neighbours
              // restore back the iphi value
              (*decoderHCalEndcap)["phi"].set(endcapCellID, iphi);
            }

            // consider the the first cell (lowest theta) from the negative-z part endcap layer
            // to be added in the neighbours list of the last barrel cell (highest theta bin)
            if(hcalEndcapSegmentation) // phi-theta segmentation
            {
              decoderHCalEndcap->set(endcapCellID, "theta",
                endcapCellIndexes[endcapCellIndexes.size() / 2]);
            }
            else // phi-row segmentation
            {
              decoderHCalEndcap->set(endcapCellID, "row",
                endcapCellIndexes[endcapCellIndexes.size() / 2]);
              // set the pseudo-layer bit field if the phi-row segmentation is used for the Endcap
              unsigned int pseudoLayer = hcalEndcapPhiRowSegmentation->definePseudoLayer(endcapCellID);
              decoderHCalEndcap->set(endcapCellID, "pseudoLayer", pseudoLayer);
            }

            // get thetaMin and thetaMax of the endcap cell
            cTheta = (hcalEndcapSegmentation) ? hcalEndcapSegmentation->cellTheta(endcapCellID) :
                                                                      hcalEndcapPhiRowSegmentation->cellTheta(endcapCellID);
            // add in the neighbours list if it overlaps in theta with the last barrel cell (highest theta bin)
            if ((cTheta[0] <= maxCellIdTheta[0] && cTheta[1] > maxCellIdTheta[0]) ||
                (cTheta[0] > maxCellIdTheta[0] && cTheta[0] < maxCellIdTheta[1])) {
              map[barrelCellID[1]].push_back(
                  endcapCellID); // add to the neighbours list of the last cell (highest theta bin) in the barrel layer
              map[endcapCellID].push_back(
                  barrelCellID[1]); // add the last barrel cell to the list of the endcap cell neighbours
              // add the next and previous phi modules as well
              // previous: if the current is 0 then previous is the last bin (id = phiBins - 1) else current - 1
              (*decoderHCalEndcap)["phi"].set(endcapCellID, (iphi == 0) ? (nPhiBins - 1) : (iphi - 1));
              map[barrelCellID[1]].push_back(
                  endcapCellID); // add to the neighbours list of the last cell in the barrel layer
              map[endcapCellID].push_back(
                  barrelCellID[1]); // add the last barrel cell to the list of the endcap cell neighbours
              // next: if the current is the last bin (id = phiBins - 1) then the next is the first bin (id = 0) else
              // current + 1
              (*decoderHCalEndcap)["phi"].set(endcapCellID, (iphi == (nPhiBins - 1)) ? 0 : (iphi + 1));
              map[barrelCellID[1]].push_back(
                  endcapCellID); // add to the neighbours list of the last cell in the barrel layer
              map[endcapCellID].push_back(
                  barrelCellID[1]); // add the last barrel cell to the list of the endcap cell neighbours
              // restore back the iphi value
              (*decoderHCalEndcap)["phi"].set(endcapCellID, iphi);
            }
          }
        }

        // loop over the part2 layers
        // consider only the first cell of each layer
        for (int part2LayerId = minLayerIdEndcap[1]; part2LayerId <= maxLayerIdEndcap[1]; part2LayerId++) {
          // set the layer bit field for part2 layer
          (*decoderHCalEndcap)[m_activeFieldNamesSegmented[hcalEndcapId]].set(endcapCellID, part2LayerId);

          // get cell indexes in the given endcap layer
          std::vector<int> endcapCellIndexes( (hcalEndcapSegmentation) ? hcalEndcapSegmentation->thetaBins(part2LayerId) :
                                                                         hcalEndcapPhiRowSegmentation->cellIndexes(part2LayerId) );

          // consider the the first cell (highest theta) from the positive-z part endcap layer
          // to be added in the neighbours list of the first barrel cell (lowest theta bin)
          if(hcalEndcapSegmentation) // phi-theta segmentation
          {
            decoderHCalEndcap->set(endcapCellID, "theta",
              endcapCellIndexes[endcapCellIndexes.size() / 2 - 1]);
          }
          else // phi-row segmentation
          {
            decoderHCalEndcap->set(endcapCellID, "row",
              endcapCellIndexes[0]);
            // set the pseudo-layer bit field if the phi-row segmentation is used for the Endcap
            unsigned int pseudoLayer = hcalEndcapPhiRowSegmentation->definePseudoLayer(endcapCellID);
            decoderHCalEndcap->set(endcapCellID, "pseudoLayer", pseudoLayer);
          }

          // we do not want to add the part2 layer cells that have part1 cells in the neighbours list
          // try to find if there is a part1 cell in the list of neighbours
          bool hasPart1Neighbour = false;
          // loop over the neighbours of the first cell
          for (auto nCellid : map[endcapCellID]) {
            int nLayerId = decoderHCalEndcap->get(nCellid, "layer");
            int systemId = decoderHCalEndcap->get(nCellid, "system");
            // check if the neighbour cell is in the endcap part1 layer
            if (nLayerId >= minLayerIdEndcap[0] && nLayerId <= maxLayerIdEndcap[0] && systemId == hcalEndcapId) {
              hasPart1Neighbour = true;
              break;
            }
          }

          // get thetaMin and thetaMax of the endcap cell
          std::array<double, 2> cTheta = (hcalEndcapSegmentation) ? hcalEndcapSegmentation->cellTheta(endcapCellID) :
                                                                    hcalEndcapPhiRowSegmentation->cellTheta(endcapCellID);
          // add in the neighbours list if it overlaps in theta with the first barrel cell (lowest theta bin)
          if ((cTheta[0] <= minCellIdTheta[0] && cTheta[1] > minCellIdTheta[0]) ||
              (cTheta[0] > minCellIdTheta[0] && cTheta[0] < minCellIdTheta[1])) {
            map[barrelCellID[0]].push_back(
                endcapCellID); // add to the neighbours list of the first cell (lowest theta bin) in the barrel layer
            map[endcapCellID].push_back(
                barrelCellID[0]); // add the first barrel cell to the list of the endcap cell neighbours
            // add the next and previous phi modules as well
            // previous: if the current is 0 then previous is the last bin (id = phiBins - 1) else current - 1
            (*decoderHCalEndcap)["phi"].set(endcapCellID, (iphi == 0) ? (nPhiBins - 1) : (iphi - 1));
            map[barrelCellID[0]].push_back(
                endcapCellID); // add to the neighbours list of the first cell in the barrel layer
            map[endcapCellID].push_back(
                barrelCellID[0]); // add the last barrel cell to the list of the endcap cell neighbours
            // next: if the current is the last bin (id = phiBins - 1) then the next is the first bin (id = 0) else
            // current + 1
            (*decoderHCalEndcap)["phi"].set(endcapCellID, (iphi == (nPhiBins - 1)) ? 0 : (iphi + 1));
            map[barrelCellID[0]].push_back(
                endcapCellID); // add to the neighbours list of the first cell in the barrel layer
            map[endcapCellID].push_back(
                barrelCellID[0]); // add the last barrel cell to the list of the endcap cell neighbours
            // restore back the iphi value
            (*decoderHCalEndcap)["phi"].set(endcapCellID, iphi);
          }

          // consider the the first cell (lowest theta) from the negative-z part endcap layer
          // to be added in the neighbours list of the last barrel cell (highest theta bin)
          if(hcalEndcapSegmentation) // phi-theta segmentation
          {
            decoderHCalEndcap->set(endcapCellID, "theta",
              endcapCellIndexes[endcapCellIndexes.size() / 2]);
          }
          else // phi-row segmentation
          {
            decoderHCalEndcap->set(endcapCellID, "row",
              endcapCellIndexes[endcapCellIndexes.size() / 2]);
              // set the pseudo-layer bit field if the phi-row segmentation is used for the Endcap
              unsigned int pseudoLayer = hcalEndcapPhiRowSegmentation->definePseudoLayer(endcapCellID);
              decoderHCalEndcap->set(endcapCellID, "pseudoLayer", pseudoLayer);
          }

          // get thetaMin and thetaMax of the endcap cell
          cTheta = (hcalEndcapSegmentation) ? hcalEndcapSegmentation->cellTheta(endcapCellID) :
                                                                      hcalEndcapPhiRowSegmentation->cellTheta(endcapCellID);
          // add in the neighbours list if it overlaps in theta with the last barrel cell (highest theta bin)
          if ((cTheta[0] <= maxCellIdTheta[0] && cTheta[1] > maxCellIdTheta[0]) ||
              (cTheta[0] > maxCellIdTheta[0] && cTheta[0] < maxCellIdTheta[1])) {
            map[barrelCellID[1]].push_back(
                endcapCellID); // add to the neighbours list of the last cell in the barrel layer
            map[endcapCellID].push_back(
                barrelCellID[1]); // add the last barrel cell to the list of the endcap cell neighbours
            // add the next and previous phi modules as well
            // previous: if the current is 0 then previous is the last bin (id = phiBins - 1) else current - 1
            (*decoderHCalEndcap)["phi"].set(endcapCellID, (iphi == 0) ? (nPhiBins - 1) : (iphi - 1));
            map[barrelCellID[1]].push_back(
                endcapCellID); // add to the neighbours list of the last cell in the barrel layer
            map[endcapCellID].push_back(
                barrelCellID[1]); // add the last barrel cell to the list of the endcap cell neighbours
            // next: if the current is the last bin (id = phiBins - 1) then the next is the first bin (id = 0) else
            // current + 1
            (*decoderHCalEndcap)["phi"].set(endcapCellID, (iphi == (nPhiBins - 1)) ? 0 : (iphi + 1));
            map[barrelCellID[1]].push_back(
                endcapCellID); // add to the neighbours list of the last cell in the barrel layer
            map[endcapCellID].push_back(
                barrelCellID[1]); // add the last barrel cell to the list of the endcap cell neighbours
            // restore back the iphi value
            (*decoderHCalEndcap)["phi"].set(endcapCellID, iphi);
          }

          // if the cell in this layer has a neighbour in the part1 layers, then stop here the part2 layers loop
          if (hasPart1Neighbour)
            break;
        }
        debug() << "Layer: " << ilayer << " iphi: " << iphi << " First Barrel CellID: " << barrelCellID[0] << ": "
                << map[barrelCellID[0]].size() - numberOfNeighbours[0]
                << " neighbours added after connecting Barrel+Endcap." << endmsg;
        debug() << "Layer: " << ilayer << " iphi: " << iphi << " Last Barrel CellID: " << barrelCellID[1] << ": "
                << map[barrelCellID[1]].size() - numberOfNeighbours[1]
                << " neighbours added after connecting Barrel+Endcap." << endmsg;
      } // iphi loop
    } // loop over the barrel layers
  }
  // FCCSWHCalPhiRow segmentation for both Barrel and Endcap
  else if (m_connectHCal && hcalBarrelPhiRowSegmentation && hcalEndcapPhiRowSegmentation) {
    info() << "Linking the HCal Barrel and HCal Endcap" << endmsg;

    // Get the part1 and part2 endcap layer indexes that will be used to connect Barrel cells to the Endcap cells
    std::vector<int> minLayerIdEndcap = {0, 0};
    std::vector<int> maxLayerIdEndcap = {0, 0};
    std::vector<std::pair<uint, uint>> minMaxLayerId(hcalEndcapPhiRowSegmentation->getMinMaxLayerId());
    if (minMaxLayerId.empty()) {
      error() << "hcalEndcap segmentation is not configured correctly, cannot link Endcap and Barrel!" << endmsg;
      return StatusCode::FAILURE;
    }

    // Loop over active layers
    for (unsigned int ilayer = 0; ilayer < m_activeVolumesNumbersSegmented[hcalBarrelId]; ilayer++) {
      debug() << "HCal Barrel Layer: " << ilayer << endmsg;
      std::array<dd4hep::DDSegmentation::CellID, 2> barrelCellID = {0, 0};
      (*decoderHCalBarrel)[m_fieldNamesSegmented[hcalBarrelId]].set(barrelCellID[0],
                                                                    m_fieldValuesSegmented[hcalBarrelId]);
      (*decoderHCalBarrel)[m_activeFieldNamesSegmented[hcalBarrelId]].set(barrelCellID[0], ilayer);
      (*decoderHCalBarrel)["row"].set(barrelCellID[0], 0);
      (*decoderHCalBarrel)["phi"].set(barrelCellID[0], 0);
      (*decoderHCalBarrel)[m_fieldNamesSegmented[hcalBarrelId]].set(barrelCellID[1],
                                                                    m_fieldValuesSegmented[hcalBarrelId]);
      (*decoderHCalBarrel)[m_activeFieldNamesSegmented[hcalBarrelId]].set(barrelCellID[1], ilayer);
      (*decoderHCalBarrel)["row"].set(barrelCellID[1], 0);
      (*decoderHCalBarrel)["phi"].set(barrelCellID[1], 0);

      std::vector<int> cellIndexes(hcalBarrelPhiRowSegmentation->cellIndexes(ilayer));
      int minCellId = cellIndexes.front(); // first cell in the barrel layer
      int maxCellId = cellIndexes.back();  // last cell in the barrel layer
      int nPhiBins = hcalBarrelPhiRowSegmentation->phiBins();
      for (int iphi = 0; iphi < nPhiBins; iphi++) {
        // set the final bits for the first cell in the barrel layer
        decoderHCalBarrel->set(barrelCellID[0], "phi", iphi);
        decoderHCalBarrel->set(barrelCellID[0], "row", minCellId);
        // set the final bits for the last cell in the barrel layer
        decoderHCalBarrel->set(barrelCellID[1], "phi", iphi);
        decoderHCalBarrel->set(barrelCellID[1], "row", maxCellId);

        // Number of neighbours for the first and the last cell in the barrel layer
        std::array<uint64_t, 2> numberOfNeighbours = {map[barrelCellID[0]].size(), map[barrelCellID[1]].size()};

        // min and max polar angle for the first cell in the barrel layer
        std::array<double, 2> minCellIdTheta = hcalBarrelPhiRowSegmentation->cellTheta(barrelCellID[0]);
        // min and max polar angle for the last cell in the barrel layer
        std::array<double, 2> maxCellIdTheta = hcalBarrelPhiRowSegmentation->cellTheta(barrelCellID[1]);

        // part1 min and max layer index
        minLayerIdEndcap[0] = minMaxLayerId[0].first;
        maxLayerIdEndcap[0] = minMaxLayerId[0].second;
        // part2 min and max layer index
        minLayerIdEndcap[1] = minMaxLayerId[1].first;
        maxLayerIdEndcap[1] = minMaxLayerId[1].second;

        // create the endcap cellID that should become the neighbour for the barrel cell
        dd4hep::DDSegmentation::CellID endcapCellID = 0;
        (*decoderHCalEndcap)[m_fieldNamesSegmented[hcalEndcapId]].set(endcapCellID,
                                                                      m_fieldValuesSegmented[hcalEndcapId]);
        (*decoderHCalEndcap)["phi"].set(endcapCellID, iphi);

        // loop over the part1 layers
        for (int part1LayerId = minLayerIdEndcap[0]; part1LayerId <= maxLayerIdEndcap[0]; part1LayerId++) {
          // set the layer bit field for part1 layer
          (*decoderHCalEndcap)[m_activeFieldNamesSegmented[hcalEndcapId]].set(endcapCellID, part1LayerId);

          // get cell indexes in the given endcap layer
          std::vector<int> endcapCellIndexes(hcalEndcapPhiRowSegmentation->cellIndexes(part1LayerId));

          // consider all cells in the first layer of part1 endcap
          if (part1LayerId == minLayerIdEndcap[0]) {
            for (auto idx : endcapCellIndexes) {
              decoderHCalEndcap->set(endcapCellID, "row", idx);
              // get thetaMin and thetaMax of the endcap cell
              std::array<double, 2> cTheta = hcalEndcapPhiRowSegmentation->cellTheta(endcapCellID);
              // add in the neighbours list if it overlaps with first barrel cell in theta
              if ((cTheta[0] <= minCellIdTheta[0] && cTheta[1] > minCellIdTheta[0]) ||
                  (cTheta[0] > minCellIdTheta[0] && cTheta[0] < minCellIdTheta[1])) {
                map[barrelCellID[0]].push_back(
                    endcapCellID); // add to the neighbours list of the first cell in the barrel layer
                map[endcapCellID].push_back(
                    barrelCellID[0]); // add the first barrel cell to the list of the endcap cell neighbours
                // add the next and previous phi modules as well
                // previous: if the current is 0 then previous is the last bin (id = phiBins - 1) else current - 1
                (*decoderHCalEndcap)["phi"].set(endcapCellID, (iphi == 0) ? (nPhiBins - 1) : (iphi - 1));
                map[barrelCellID[0]].push_back(
                    endcapCellID); // add to the neighbours list of the first cell in the barrel layer
                map[endcapCellID].push_back(
                    barrelCellID[0]); // add the last barrel cell to the list of the endcap cell neighbours
                // next: if the current is the last bin (id = phiBins - 1) then the next is the first bin (id = 0) else
                // current + 1
                (*decoderHCalEndcap)["phi"].set(endcapCellID, (iphi == (nPhiBins - 1)) ? 0 : (iphi + 1));
                map[barrelCellID[0]].push_back(
                    endcapCellID); // add to the neighbours list of the first cell in the barrel layer
                map[endcapCellID].push_back(
                    barrelCellID[0]); // add the last barrel cell to the list of the endcap cell neighbours
                // restore back the iphi value
                (*decoderHCalEndcap)["phi"].set(endcapCellID, iphi);
              }
              // add in the neighbours list if it overlaps with first barrel cell in theta
              if ((cTheta[0] <= maxCellIdTheta[0] && cTheta[1] > maxCellIdTheta[0]) ||
                  (cTheta[0] > maxCellIdTheta[0] && cTheta[0] < maxCellIdTheta[1])) {
                map[barrelCellID[1]].push_back(
                    endcapCellID); // add to the neighbours list of the last cell in the barrel layer
                map[endcapCellID].push_back(
                    barrelCellID[1]); // add the last barrel cell to the list of the endcap cell neighbours
                // add the next and previous phi modules as well
                // previous: if the current is 0 then previous is the last bin (id = phiBins - 1) else current - 1
                (*decoderHCalEndcap)["phi"].set(endcapCellID, (iphi == 0) ? (nPhiBins - 1) : (iphi - 1));
                map[barrelCellID[1]].push_back(
                    endcapCellID); // add to the neighbours list of the last cell in the barrel layer
                map[endcapCellID].push_back(
                    barrelCellID[1]); // add the last barrel cell to the list of the endcap cell neighbours
                // next: if the current is the last bin (id = phiBins - 1) then the next is the first bin (id = 0) else
                // current + 1
                (*decoderHCalEndcap)["phi"].set(endcapCellID, (iphi == (nPhiBins - 1)) ? 0 : (iphi + 1));
                map[barrelCellID[1]].push_back(
                    endcapCellID); // add to the neighbours list of the last cell in the barrel layer
                map[endcapCellID].push_back(
                    barrelCellID[1]); // add the last barrel cell to the list of the endcap cell neighbours
                // restore back the iphi value
                (*decoderHCalEndcap)["phi"].set(endcapCellID, iphi);
              }
            }
          } else // consider only the first cell
          {
            std::array<double, 2> cTheta;
            decoderHCalEndcap->set(
                endcapCellID, "row",
                endcapCellIndexes[endcapCellIndexes.size() / 2]); // first cell in the negative-z part endcap layer
            // add in the neighbours list if it overlaps with first barrel cell in theta
            cTheta =
                hcalEndcapPhiRowSegmentation->cellTheta(endcapCellID); // get thetaMin and thetaMax of the endcap cell
            if ((cTheta[0] <= minCellIdTheta[0] && cTheta[1] > minCellIdTheta[0]) ||
                (cTheta[0] > minCellIdTheta[0] && cTheta[0] < minCellIdTheta[1])) {
              map[barrelCellID[0]].push_back(
                  endcapCellID); // add to the neighbours list of the first cell in the barrel layer
              map[endcapCellID].push_back(
                  barrelCellID[0]); // add the first barrel cell to the list of the endcap cell neighbours
              // add the next and previous phi modules as well
              // previous: if the current is 0 then previous is the last bin (id = phiBins - 1) else current - 1
              (*decoderHCalEndcap)["phi"].set(endcapCellID, (iphi == 0) ? (nPhiBins - 1) : (iphi - 1));
              map[barrelCellID[0]].push_back(
                  endcapCellID); // add to the neighbours list of the first cell in the barrel layer
              map[endcapCellID].push_back(
                  barrelCellID[0]); // add the last barrel cell to the list of the endcap cell neighbours
              // next: if the current is the last bin (id = phiBins - 1) then the next is the first bin (id = 0) else
              // current + 1
              (*decoderHCalEndcap)["phi"].set(endcapCellID, (iphi == (nPhiBins - 1)) ? 0 : (iphi + 1));
              map[barrelCellID[0]].push_back(
                  endcapCellID); // add to the neighbours list of the first cell in the barrel layer
              map[endcapCellID].push_back(
                  barrelCellID[0]); // add the last barrel cell to the list of the endcap cell neighbours
              // restore back the iphi value
              (*decoderHCalEndcap)["phi"].set(endcapCellID, iphi);
            }
            decoderHCalEndcap->set(endcapCellID, "row",
                                   endcapCellIndexes[0]); // first cell in the positive-z part endcap layer
            // add in the neighbours list if it overlaps with the last barrel cell in theta
            cTheta =
                hcalEndcapPhiRowSegmentation->cellTheta(endcapCellID); // get thetaMin and thetaMax of the endcap cell
            if ((cTheta[0] <= maxCellIdTheta[0] && cTheta[1] > maxCellIdTheta[0]) ||
                (cTheta[0] > maxCellIdTheta[0] && cTheta[0] < maxCellIdTheta[1])) {
              map[barrelCellID[1]].push_back(
                  endcapCellID); // add to the neighbours list of the last cell in the barrel layer
              map[endcapCellID].push_back(
                  barrelCellID[1]); // add the last barrel cell to the list of the endcap cell neighbours
              // add the next and previous phi modules as well
              // previous: if the current is 0 then previous is the last bin (id = phiBins - 1) else current - 1
              (*decoderHCalEndcap)["phi"].set(endcapCellID, (iphi == 0) ? (nPhiBins - 1) : (iphi - 1));
              map[barrelCellID[1]].push_back(
                  endcapCellID); // add to the neighbours list of the last cell in the barrel layer
              map[endcapCellID].push_back(
                  barrelCellID[1]); // add the last barrel cell to the list of the endcap cell neighbours
              // next: if the current is the last bin (id = phiBins - 1) then the next is the first bin (id = 0) else
              // current + 1
              (*decoderHCalEndcap)["phi"].set(endcapCellID, (iphi == (nPhiBins - 1)) ? 0 : (iphi + 1));
              map[barrelCellID[1]].push_back(
                  endcapCellID); // add to the neighbours list of the last cell in the barrel layer
              map[endcapCellID].push_back(
                  barrelCellID[1]); // add the last barrel cell to the list of the endcap cell neighbours
              // restore back the iphi value
              (*decoderHCalEndcap)["phi"].set(endcapCellID, iphi);
            }
          }
        }

        // loop over the part2 layers
        // consider only the first cell of each layer
        for (int part2LayerId = minLayerIdEndcap[1]; part2LayerId <= maxLayerIdEndcap[1]; part2LayerId++) {
          // set the layer bit field for part2 layer
          (*decoderHCalEndcap)[m_activeFieldNamesSegmented[hcalEndcapId]].set(endcapCellID, part2LayerId);

          // get cell indexes in the given endcap layer
          std::vector<int> endcapCellIndexes(hcalEndcapPhiRowSegmentation->cellIndexes(part2LayerId));

          std::array<double, 2> cTheta;
          decoderHCalEndcap->set(
              endcapCellID, "row",
              endcapCellIndexes[endcapCellIndexes.size() / 2]); // first cell in the negative-z part endcap layer

          // we do not want to add the part2 layer cells that have part1 cells in the neighbours list
          // try to find if there is a part1 cell in the list of neighbours
          bool hasPart1Neighbour = false;
          // loop over the neighbours of the first cell
          for (auto nCellid : map[endcapCellID]) {
            int nLayerId = decoderHCalEndcap->get(nCellid, "layer");
            int systemId = decoderHCalEndcap->get(nCellid, "system");
            // check if the neighbour cell is in the endcap part1 layer
            if (nLayerId >= minLayerIdEndcap[0] && nLayerId <= maxLayerIdEndcap[0] && systemId == hcalEndcapId) {
              hasPart1Neighbour = true;
              break;
            }
          }

          // add in the neighbours list if it overlaps with first barrel cell in theta
          cTheta =
              hcalEndcapPhiRowSegmentation->cellTheta(endcapCellID); // get thetaMin and thetaMax of the endcap cell
          if ((cTheta[0] <= minCellIdTheta[0] && cTheta[1] > minCellIdTheta[0]) ||
              (cTheta[0] > minCellIdTheta[0] && cTheta[0] < minCellIdTheta[1])) {
            map[barrelCellID[0]].push_back(
                endcapCellID); // add to the neighbours list of the first cell in the barrel layer
            map[endcapCellID].push_back(
                barrelCellID[0]); // add the first barrel cell to the list of the endcap cell neighbours
            // add the next and previous phi modules as well
            // previous: if the current is 0 then previous is the last bin (id = phiBins - 1) else current - 1
            (*decoderHCalEndcap)["phi"].set(endcapCellID, (iphi == 0) ? (nPhiBins - 1) : (iphi - 1));
            map[barrelCellID[0]].push_back(
                endcapCellID); // add to the neighbours list of the first cell in the barrel layer
            map[endcapCellID].push_back(
                barrelCellID[0]); // add the last barrel cell to the list of the endcap cell neighbours
            // next: if the current is the last bin (id = phiBins - 1) then the next is the first bin (id = 0) else
            // current + 1
            (*decoderHCalEndcap)["phi"].set(endcapCellID, (iphi == (nPhiBins - 1)) ? 0 : (iphi + 1));
            map[barrelCellID[0]].push_back(
                endcapCellID); // add to the neighbours list of the first cell in the barrel layer
            map[endcapCellID].push_back(
                barrelCellID[0]); // add the last barrel cell to the list of the endcap cell neighbours
            // restore back the iphi value
            (*decoderHCalEndcap)["phi"].set(endcapCellID, iphi);
          }
          decoderHCalEndcap->set(endcapCellID, "row",
                                 endcapCellIndexes[0]); // first cell in the positive-z part endcap layer
          // add in the neighbours list if it overlaps with the last barrel cell in theta
          cTheta =
              hcalEndcapPhiRowSegmentation->cellTheta(endcapCellID); // get thetaMin and thetaMax of the endcap cell
          if ((cTheta[0] <= maxCellIdTheta[0] && cTheta[1] > maxCellIdTheta[0]) ||
              (cTheta[0] > maxCellIdTheta[0] && cTheta[0] < maxCellIdTheta[1])) {
            map[barrelCellID[1]].push_back(
                endcapCellID); // add to the neighbours list of the last cell in the barrel layer
            map[endcapCellID].push_back(
                barrelCellID[1]); // add the last barrel cell to the list of the endcap cell neighbours
            // add the next and previous phi modules as well
            // previous: if the current is 0 then previous is the last bin (id = phiBins - 1) else current - 1
            (*decoderHCalEndcap)["phi"].set(endcapCellID, (iphi == 0) ? (nPhiBins - 1) : (iphi - 1));
            map[barrelCellID[1]].push_back(
                endcapCellID); // add to the neighbours list of the last cell in the barrel layer
            map[endcapCellID].push_back(
                barrelCellID[1]); // add the last barrel cell to the list of the endcap cell neighbours
            // next: if the current is the last bin (id = phiBins - 1) then the next is the first bin (id = 0) else
            // current + 1
            (*decoderHCalEndcap)["phi"].set(endcapCellID, (iphi == (nPhiBins - 1)) ? 0 : (iphi + 1));
            map[barrelCellID[1]].push_back(
                endcapCellID); // add to the neighbours list of the last cell in the barrel layer
            map[endcapCellID].push_back(
                barrelCellID[1]); // add the last barrel cell to the list of the endcap cell neighbours
            // restore back the iphi value
            (*decoderHCalEndcap)["phi"].set(endcapCellID, iphi);
          }

          // if the cell in this layer has a neighbour in the part1 layers, then stop here the part2 layers loop
          if (hasPart1Neighbour)
            break;
        }
        debug() << "Layer: " << ilayer << " iphi: " << iphi << " First Barrel CellID: " << barrelCellID[0] << ": "
                << map[barrelCellID[0]].size() - numberOfNeighbours[0]
                << " neighbours added after connecting Barrel+Endcap." << endmsg;
        debug() << "Layer: " << ilayer << " iphi: " << iphi << " Last Barrel CellID: " << barrelCellID[1] << ": "
                << map[barrelCellID[1]].size() - numberOfNeighbours[1]
                << " neighbours added after connecting Barrel+Endcap." << endmsg;
      } // iphi loop
    } // loop over the barrel layers
  }
  else if (m_connectHCal)
  {
    error() << "Unable to connect the HCal Barrel and HCal Endcap! segmentations are not handled correctly. -->\n"
            << "Please provide one of the following:\n"
            << "  PhiTheta segmentation readout for both HCal Barrel and HCal Endcap;\n"
            << "  PhiRow segmentation readout for both HCal Barrel and HCal Endcap;\n"
            << "  PhiTheta segmentation readout for HCal Barrel while PhiRow segmentation for HCal Endcap." << endmsg;
    return StatusCode::FAILURE;
  }

  //////////////////////////////////////////////////
  ///      BARREL: connection ECAL + HCAL        ///
  /////////////////////////////////////////////////

  std::set<dd4hep::DDSegmentation::CellID> linkedCells;

  if (m_connectBarrels) {
    info() << "Linking the ECal and HCal Barrels" << endmsg;
    // first check if ECAL barrel and HCal barrel are configured
    if (decoderECalBarrel == nullptr || decoderHCalBarrel == nullptr) {
      error() << "ECAL barrel and/or HCal barrel are not configured correctly, cannot link barrels!" << endmsg;
      return StatusCode::FAILURE;
    }

    // check which type of segmentation is available for HCal barrel
    bool thetaSegHCal = false;
    if (hcalBarrelSegmentation && hcalBarrelPhiRowSegmentation) {
      error() << "Please provide either PhiTheta or PhiRow segmentation readout for HCal..." << endmsg;
      return StatusCode::FAILURE;
    } else if (hcalBarrelSegmentation)
      thetaSegHCal = true;

    // print how many cells in each dimensions will be matched
    info() << "ECAL layer " << eCalLastLayer << " is a neighbour of HCAL layer 0." << endmsg;
    if (ecalBarrelModuleThetaSegmentation) {
      info() << "ECAL modules " << extremaECalLastLayerModule.first << " - " << extremaECalLastLayerModule.second;
    } else {
      info() << "ECAL phi cells " << extremaECalLastLayerPhi.first << " - " << extremaECalLastLayerPhi.second;
    }
    info() << " will be matched to HCAL phi cells " << extremaHCalFirstLayerPhi.first << " - "
           << extremaHCalFirstLayerPhi.second << endmsg;
    info() << "ECAL theta cells " << extremaECalLastLayerTheta.first << " - " << extremaECalLastLayerTheta.second
           << " will be matched to HCAL theta cells " << extremaHCalFirstLayerTheta.first << " - "
           << extremaHCalFirstLayerTheta.second << endmsg;

    dd4hep::DDSegmentation::CellID volumeIdECal = 0;
    dd4hep::DDSegmentation::CellID volumeIdHCal = 0;
    (*decoderECalBarrel)["system"].set(volumeIdECal, m_ecalBarrelSysId);
    (*decoderECalBarrel)[m_activeFieldNamesSegmented[ecalBarrelId]].set(volumeIdECal, eCalLastLayer);
    (*decoderHCalBarrel)["system"].set(volumeIdHCal, m_hcalBarrelSysId);
    (*decoderHCalBarrel)[m_activeFieldNamesSegmented[hcalBarrelId]].set(volumeIdHCal, 0);

    // retrieve needed parameters
    double hCalPhiSize =
        (thetaSegHCal) ? hcalBarrelSegmentation->gridSizePhi() : hcalBarrelPhiRowSegmentation->gridSizePhi();
    double eCalThetaOffset, eCalThetaSize;
    double eCalPhiOffset, eCalPhiSize;
    int eCalModules(-1);
    int eCalPhiBins(-1);
    if (ecalBarrelModuleThetaSegmentation) {
      eCalThetaOffset = ecalBarrelModuleThetaSegmentation->offsetTheta();
      eCalThetaSize = ecalBarrelModuleThetaSegmentation->gridSizeTheta();
      eCalModules = ecalBarrelModuleThetaSegmentation->nModules();
      // for ECAL barrel with module readout, need to find out phi position of module 0 of ECal barrel last layer
      // get it from volume manager
      // for a normal phi grid it would suffice to use the segmentation phi(cellId) method
      dd4hep::VolumeManager volman = m_geoSvc->getDetector()->volumeManager();
      auto detelement = volman.lookupDetElement(volumeIdECal);
      double inLocal[] = {0, 0, 0};
      const auto outGlobal = detelement.nominal().localToWorld(inLocal);
      eCalPhiOffset = std::atan2(outGlobal.Y(), outGlobal.X());
      info() << "ECAL modules : " << eCalModules << " , phi of module 0 in last ECAL layer : " << eCalPhiOffset
             << endmsg;
      eCalPhiSize = TMath::TwoPi() / eCalModules;
    } else {
      eCalThetaOffset = ecalBarrelPhiThetaSegmentation->offsetTheta();
      eCalThetaSize = ecalBarrelPhiThetaSegmentation->gridSizeTheta();
      eCalPhiOffset = ecalBarrelPhiThetaSegmentation->offsetPhi();
      eCalPhiSize = ecalBarrelPhiThetaSegmentation->gridSizePhi();
      eCalPhiBins = ecalBarrelPhiThetaSegmentation->phiBins();
    }

    dd4hep::DDSegmentation::CellID cellIdECal = volumeIdECal;
    dd4hep::DDSegmentation::CellID cellIdHCal = volumeIdHCal;

    // Loop over HCAL segmentation cells to find neighbours in ECAL
    // Loop on phi
    for (int iphi = extremaHCalFirstLayerPhi.first; iphi <= extremaHCalFirstLayerPhi.second; iphi++) {
      // determine phi extent of HCAL cell
      (*decoderHCalBarrel)["phi"].set(cellIdHCal, iphi);
      double phi =
          (thetaSegHCal) ? hcalBarrelSegmentation->phi(cellIdHCal) : hcalBarrelPhiRowSegmentation->phi(cellIdHCal);
      double phiMin = phi - 0.5 * hCalPhiSize;
      double phiMax = phi + 0.5 * hCalPhiSize;

      // find ECAL barrel modules corresponding to this phi range
      int minModuleECal(-1), maxModuleECal(-1), minPhiECal(-1), maxPhiECal(-1);
      if (ecalBarrelModuleThetaSegmentation) {
        minModuleECal = int(floor((phiMin - eCalPhiOffset + 0.5 * eCalPhiSize) / eCalPhiSize));
        if (minModuleECal < 0)
          minModuleECal +=
              eCalModules; // need module to be >=0 so that subtracting module%mergedModules gives the correct result
        minModuleECal -= minModuleECal % ecalBarrelModuleThetaSegmentation->mergedModules(eCalLastLayer);
        maxModuleECal = int(floor((phiMax - eCalPhiOffset + 0.5 * eCalPhiSize) / eCalPhiSize));
        if (maxModuleECal < 0)
          maxModuleECal += eCalModules;
        maxModuleECal -= maxModuleECal % ecalBarrelModuleThetaSegmentation->mergedModules(eCalLastLayer);
        // due to ciclic behaviour of modules, one could have e.g min = 1534 and max = 2 (for N=1536)
        if (maxModuleECal < minModuleECal)
          maxModuleECal += eCalModules;
      } else {
        minPhiECal = int(floor((phiMin - eCalPhiOffset + 0.5 * eCalPhiSize) / eCalPhiSize));
        maxPhiECal = int(floor((phiMax - eCalPhiOffset + 0.5 * eCalPhiSize) / eCalPhiSize));
      }

      if (thetaSegHCal) {
        // get all cells (theta bins) from the first layer of HCal barrel
        std::vector<int> hcalThetaBins(hcalBarrelSegmentation->thetaBins(0));

        for (auto itheta : hcalThetaBins) {
          (*decoderHCalBarrel)["theta"].set(cellIdHCal, itheta);
          // min and max polar angle of the hcal cell
          std::array<double, 2> hcalCellTheta = hcalBarrelSegmentation->cellTheta(cellIdHCal);
          double theta = hcalBarrelSegmentation->theta(cellIdHCal);
          double thetaMin = hcalCellTheta[0];
          double thetaMax = hcalCellTheta[1];

          // find ECAL barrel theta bins corresponding to this theta range
          int minThetaECal = int(floor((thetaMin - eCalThetaOffset + 0.5 * eCalThetaSize) / eCalThetaSize));
          int maxThetaECal = int(floor((thetaMax - eCalThetaOffset + 0.5 * eCalThetaSize) / eCalThetaSize));
          if (ecalBarrelModuleThetaSegmentation) {
            minThetaECal -= minThetaECal % ecalBarrelModuleThetaSegmentation->mergedThetaCells(eCalLastLayer);
            maxThetaECal -= maxThetaECal % ecalBarrelModuleThetaSegmentation->mergedThetaCells(eCalLastLayer);
          }
          // print out some debug info
          debug() << "HCAL cell: cellId = " << cellIdHCal << " phi bin = " << iphi << " theta bin = " << itheta
                  << endmsg;
          debug() << "HCAL cell: phi = " << phi << " theta = " << theta << endmsg;
          debug() << "HCAL cell: thetaMin = " << thetaMin << " thetaMax = " << thetaMax << endmsg;
          debug() << "HCAL cell: phiMin = " << phiMin << " phiMax = " << phiMax << endmsg;
          debug() << "HCAL cell: theta bin of neighbours in ECAL layer max = " << minThetaECal << " - " << maxThetaECal
                  << endmsg;
          if (ecalBarrelModuleThetaSegmentation) {
            debug() << "HCAL cell: module of neighbours in ECAL layer max = " << minModuleECal << " - " << maxModuleECal
                    << endmsg;
          } else {
            debug() << "HCAL cell: phi bin of neighbours in ECAL layer max = " << minPhiECal << " - " << maxPhiECal
                    << endmsg;
          }

          // add the neighbours if within the acceptance of the layer
          int thetaStep = (ecalBarrelModuleThetaSegmentation)
                              ? ecalBarrelModuleThetaSegmentation->mergedThetaCells(eCalLastLayer)
                              : 1;
          for (int ithetaECal = minThetaECal; ithetaECal <= maxThetaECal; ithetaECal += thetaStep) {
            if (ithetaECal < extremaECalLastLayerTheta.first || ithetaECal > extremaECalLastLayerTheta.second) {
              debug() << "theta bin " << ithetaECal << " out of range, skipping" << endmsg;
              continue;
            }
            (*decoderECalBarrel)["theta"].set(cellIdECal, ithetaECal);
            if (ecalBarrelModuleThetaSegmentation) {
              for (int iModule = minModuleECal; iModule <= maxModuleECal;
                   iModule += ecalBarrelModuleThetaSegmentation->mergedModules(eCalLastLayer)) {
                int module = iModule;
                if (module < 0)
                  module += eCalModules;
                if (module >= eCalModules)
                  module -= eCalModules;
                if (module < extremaECalLastLayerModule.first || module > extremaECalLastLayerModule.second) {
                  warning() << "module " << module << " out of range, skipping" << endmsg;
                  continue;
                }
                (*decoderECalBarrel)["module"].set(cellIdECal, module);
                debug() << "HCAL cell: neighbour in ECAL has cellID = " << cellIdECal << endmsg;
                map.find((uint64_t)cellIdHCal)->second.push_back((uint64_t)cellIdECal);
                map.find((uint64_t)cellIdECal)->second.push_back((uint64_t)cellIdHCal);
                linkedCells.insert(cellIdHCal);
                linkedCells.insert(cellIdECal);
              }
            } else {
              for (int iphiECal = minPhiECal; iphiECal <= maxPhiECal; iphiECal++) {
                int phiBinECal = iphiECal;
                // bring phi bin in 0..nbins-1 range
                if (phiBinECal < 0)
                  phiBinECal += eCalPhiBins;
                if (phiBinECal >= eCalPhiBins)
                  phiBinECal -= eCalPhiBins;
                if (phiBinECal < extremaECalLastLayerPhi.first || phiBinECal > extremaECalLastLayerPhi.second) {
                  warning() << "phi bin " << phiBinECal << " out of range, skipping" << endmsg;
                  continue;
                }
                (*decoderECalBarrel)["phi"].set(cellIdECal, phiBinECal);
                debug() << "HCAL cell: neighbour in ECAL has cellID = " << cellIdECal << endmsg;
                map.find((uint64_t)cellIdHCal)->second.push_back((uint64_t)cellIdECal);
                map.find((uint64_t)cellIdECal)->second.push_back((uint64_t)cellIdHCal);
                linkedCells.insert(cellIdHCal);
                linkedCells.insert(cellIdECal);
              }
            }
          }
        }
      } else // HCal phi-row segmentaion
      {
        // get all cell indexes from the first layer of HCal barrel
        std::vector<int> hcalCellIndexes(hcalBarrelPhiRowSegmentation->cellIndexes(0));

        // find which HCal cells overlap in theta with the ECal cell
        for (auto idx : hcalCellIndexes) {
          (*decoderHCalBarrel)["row"].set(cellIdHCal, idx);
          // min and max polar angle of the hcal cell
          std::array<double, 2> hcalCellTheta = hcalBarrelPhiRowSegmentation->cellTheta(cellIdHCal);
          double theta = hcalBarrelPhiRowSegmentation->theta(cellIdHCal);
          double thetaMin = hcalCellTheta[0];
          double thetaMax = hcalCellTheta[1];

          // find ECAL barrel theta bins corresponding to this theta range
          int minThetaECal = int(floor((thetaMin - eCalThetaOffset + 0.5 * eCalThetaSize) / eCalThetaSize));
          int maxThetaECal = int(floor((thetaMax - eCalThetaOffset + 0.5 * eCalThetaSize) / eCalThetaSize));
          if (ecalBarrelModuleThetaSegmentation) {
            minThetaECal -= minThetaECal % ecalBarrelModuleThetaSegmentation->mergedThetaCells(eCalLastLayer);
            maxThetaECal -= maxThetaECal % ecalBarrelModuleThetaSegmentation->mergedThetaCells(eCalLastLayer);
          }
          // print out some debug info
          debug() << "HCAL cell: cellId = " << cellIdHCal << " phi bin = " << iphi << " row/cell index = " << idx
                  << endmsg;
          debug() << "HCAL cell: phi = " << phi << " theta = " << theta << endmsg;
          debug() << "HCAL cell: thetaMin = " << thetaMin << " thetaMax = " << thetaMax << endmsg;
          debug() << "HCAL cell: phiMin = " << phiMin << " phiMax = " << phiMax << endmsg;
          debug() << "HCAL cell: theta bin of neighbours in ECAL layer max = " << minThetaECal << " - " << maxThetaECal
                  << endmsg;
          if (ecalBarrelModuleThetaSegmentation) {
            debug() << "HCAL cell: module of neighbours in ECAL layer max = " << minModuleECal << " - " << maxModuleECal
                    << endmsg;
          } else {
            debug() << "HCAL cell: phi bin of neighbours in ECAL layer max = " << minPhiECal << " - " << maxPhiECal
                    << endmsg;
          }

          // add the neighbours if within the acceptance of the layer
          int thetaStep = (ecalBarrelModuleThetaSegmentation)
                              ? ecalBarrelModuleThetaSegmentation->mergedThetaCells(eCalLastLayer)
                              : 1;
          for (int ithetaECal = minThetaECal; ithetaECal <= maxThetaECal; ithetaECal += thetaStep) {
            if (ithetaECal < extremaECalLastLayerTheta.first || ithetaECal > extremaECalLastLayerTheta.second) {
              debug() << "theta bin " << ithetaECal << " out of range, skipping" << endmsg;
              continue;
            }
            (*decoderECalBarrel)["theta"].set(cellIdECal, ithetaECal);
            if (ecalBarrelModuleThetaSegmentation) {
              for (int iModule = minModuleECal; iModule <= maxModuleECal;
                   iModule += ecalBarrelModuleThetaSegmentation->mergedModules(eCalLastLayer)) {
                int module = iModule;
                if (module < 0)
                  module += eCalModules;
                if (module >= eCalModules)
                  module -= eCalModules;
                if (module < extremaECalLastLayerModule.first || module > extremaECalLastLayerModule.second) {
                  warning() << "module " << module << " out of range, skipping" << endmsg;
                  continue;
                }
                (*decoderECalBarrel)["module"].set(cellIdECal, module);
                debug() << "HCAL cell: neighbour in ECAL has cellID = " << cellIdECal << endmsg;
                map.find((uint64_t)cellIdHCal)->second.push_back((uint64_t)cellIdECal);
                map.find((uint64_t)cellIdECal)->second.push_back((uint64_t)cellIdHCal);
                linkedCells.insert(cellIdHCal);
                linkedCells.insert(cellIdECal);
              }
            } else {
              for (int iphiECal = minPhiECal; iphiECal <= maxPhiECal; iphiECal++) {
                int phiBinECal = iphiECal;
                // bring phi bin in 0..nbins-1 range
                if (phiBinECal < 0)
                  phiBinECal += eCalPhiBins;
                if (phiBinECal >= eCalPhiBins)
                  phiBinECal -= eCalPhiBins;
                if (phiBinECal < extremaECalLastLayerPhi.first || phiBinECal > extremaECalLastLayerPhi.second) {
                  warning() << "phi bin " << phiBinECal << " out of range, skipping" << endmsg;
                  continue;
                }
                (*decoderECalBarrel)["phi"].set(cellIdECal, phiBinECal);
                debug() << "HCAL cell: neighbour in ECAL has cellID = " << cellIdECal << endmsg;
                map.find((uint64_t)cellIdHCal)->second.push_back((uint64_t)cellIdECal);
                map.find((uint64_t)cellIdECal)->second.push_back((uint64_t)cellIdHCal);
                linkedCells.insert(cellIdHCal);
                linkedCells.insert(cellIdECal);
              }
            }
          }
        }
      }
    }
  } // end of connectBarrels

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
  debug() << "cells with neighbours across Calo boundaries: " << linkedCells.size() << endmsg;

  // Check if output directory exists
  std::string outDirPath = gSystem->DirName(m_outputFileName.c_str());
  if (!gSystem->OpenDirectory(outDirPath.c_str())) {
    error() << "Output directory \"" << outDirPath << "\" does not exists! Please create it." << endmsg;
    return StatusCode::FAILURE;
  }

  std::unique_ptr<TFile> outFile(TFile::Open(m_outputFileName.c_str(), "RECREATE"));
  outFile->cd();
  TTree tree("neighbours", "Tree with map of neighbours");
  uint64_t saveCellId;
  std::vector<uint64_t> saveNeighbours;
  tree.Branch("cellId", &saveCellId, "cellId/l");
  tree.Branch("neighbours", &saveNeighbours);
  for (const auto& item : map) {
    saveCellId = item.first;
    saveNeighbours = item.second;
    tree.Fill();
  }
  outFile->Write();
  outFile->Close();

  return StatusCode::SUCCESS;
}

StatusCode CreateFCCeeCaloNeighbours::initialize_lookups() {

  for (uint iSys = 0; iSys < m_readoutNamesSegmented.size(); iSys++) {
    // get segmentation
    dd4hep::DDSegmentation::Segmentation* aSegmentation =
        m_geoSvc->getDetector()->readout(m_readoutNamesSegmented[iSys]).segmentation().segmentation();
    if (aSegmentation == nullptr) {
      error() << "Segmentation does not exist." << endmsg;
      return StatusCode::FAILURE;
    }

    std::string segmentationType = aSegmentation->type();
    if (segmentationType == "FCCSWEndcapTurbine_k4geo") {
      // initialize lookups for endcap cells
      info() << "Creating lookup tables for ECal endcap segmentation of type " << segmentationType << endmsg;
      auto ecalEndcapTurbineSegmentation =
          dynamic_cast<dd4hep::DDSegmentation::FCCSWEndcapTurbine_k4geo*>(aSegmentation);
      auto endcapDecoder = m_geoSvc->getDetector()->readout(m_readoutNamesSegmented[iSys]).idSpec().decoder();
      unsigned layerOffset[3];
      layerOffset[0] = 0;
      layerOffset[1] =
          ecalEndcapTurbineSegmentation->numCellsRhoCalib(0) * ecalEndcapTurbineSegmentation->numCellsZCalib(0);
      layerOffset[2] = layerOffset[1] + ecalEndcapTurbineSegmentation->numCellsRhoCalib(1) *
                                            ecalEndcapTurbineSegmentation->numCellsZCalib(1);
      int iWheel = 2;

      unsigned numECModules = ecalEndcapTurbineSegmentation->nModules(iWheel);
      unsigned numECCellsRho = ecalEndcapTurbineSegmentation->numCellsRho(iWheel);
      info() << "wheel = " << iWheel << ", modules = " << numECModules << ", rho cells = " << numECCellsRho << endmsg;

      // GM: dont think this is actually needed since we do push_back later..
      m_EMEC_h_module_vs_phi_pos.reserve(numECCellsRho);
      m_EMEC_h_module_vs_phi_neg.reserve(numECCellsRho);

      // unused vectors
      //m_EMEC_pos_phi_lookup.reserve(numECCellsRho);
      //m_EMEC_neg_phi_lookup.reserve(numECCellsRho);

      // start sizing the vectors appropriately
      //for (unsigned iECRho = 0; iECRhoe < numECCellsRho; iECRho++) {
      //  m_EMEC_pos_phi_lookup[iECRho].reserve(numECModules);
      //  m_EMEC_neg_phi_lookup[iECRho].reserve(numECModules);
      //}

      // Loop over segmentation cells
      unsigned iECz = 0;

      dd4hep::DDSegmentation::CellID endcapCellId = 0;

      (*endcapDecoder)[m_fieldNamesSegmented[iSys]].set(endcapCellId, m_fieldValuesSegmented[iSys]);
      (*endcapDecoder)["z"].set(endcapCellId, iECz);
      (*endcapDecoder)["wheel"].set(endcapCellId, iWheel);

      for (int iSide = -1; iSide < 2; iSide += 2) {
        (*endcapDecoder)["side"].set(endcapCellId, iSide);

        for (unsigned iECrho = 0; iECrho < numECCellsRho; iECrho++) {
          TH1F* h_neg = 0;
          TH1F* h_pos = 0;
          if (iSide ==-1) {
            std::string histname_neg = "h_modulevphi_neg_" + std::to_string(iECrho);
            h_neg = new TH1F(histname_neg.c_str(), histname_neg.c_str(), numECModules, -TMath::Pi(), TMath::Pi());
            m_EMEC_h_module_vs_phi_neg.push_back(h_neg);
          }
          else {
            std::string histname_pos = "h_modulevphi_pos_" + std::to_string(iECrho);
            h_pos = new TH1F(histname_pos.c_str(), histname_pos.c_str(), numECModules, -TMath::Pi(), TMath::Pi());
            m_EMEC_h_module_vs_phi_pos.push_back(h_pos);
          }
          for (unsigned int iECmodule = 0; iECmodule < numECModules; iECmodule++) {
            (*endcapDecoder)["module"].set(endcapCellId, iECmodule);
            (*endcapDecoder)["rho"].set(endcapCellId, iECrho);
            (*endcapDecoder)[m_activeFieldNamesSegmented[iSys]].set(endcapCellId, ecalEndcapTurbineSegmentation->expLayer(iWheel, iECrho, iECz));

            double endcapPhi = TMath::ATan2(ecalEndcapTurbineSegmentation->position(endcapCellId).y(),
                                            ecalEndcapTurbineSegmentation->position(endcapCellId).x());
            double endcapRho = ecalEndcapTurbineSegmentation->rho(endcapCellId);
            double endcapZ = ecalEndcapTurbineSegmentation->z(endcapCellId);
            double endcapTheta = TMath::ATan2(endcapRho, endcapZ);

            if (iSide == -1) {
              //m_EMEC_neg_phi_lookup[iECrho].push_back(endcapPhi);
              h_neg->SetBinContent(h_neg->FindBin(endcapPhi), iECmodule);
            } else {
              //m_EMEC_pos_phi_lookup[iECrho].push_back(endcapPhi);
              h_pos->SetBinContent(h_pos->FindBin(endcapPhi), iECmodule);
            }
            if (iECmodule == 0 && iSide == 1)
              m_EMEC_theta_lookup.push_back(endcapTheta);
          }
        }
      }
    }

    else if (segmentationType == "FCCSWGridModuleThetaMerged_k4geo") {
      // initialize lookups for barrel cells
      info() << "Creating lookup tables for ECal barrel segmentation of type " << segmentationType << endmsg;
      auto moduleThetaSegmentation =
          dynamic_cast<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo*>(aSegmentation);
      auto barrelDecoder = m_geoSvc->getDetector()->readout(m_readoutNamesSegmented[iSys]).idSpec().decoder();
      m_geoSvc->getDetector()->readouts().find(m_readoutNamesSegmented[iSys]);
      // Loop over all relevant cells in the barrel and find position
      unsigned nBarrelLayers = m_activeVolumesNumbersSegmented[iSys];
      m_EMB_phi_lookup.reserve(nBarrelLayers);
      auto numCells = det::utils::numberOfCells(m_ecalBarrelSysId, *moduleThetaSegmentation);

      for (unsigned iBarrelLayer = 0; iBarrelLayer < nBarrelLayers; iBarrelLayer++) {
        std::vector<Float_t> phis, thetas;
        unsigned nModulesBarrel = (numCells[0] - 1) * moduleThetaSegmentation->mergedModules(iBarrelLayer);
        std::string histname = "h_EMB_modulevphi_" + std::to_string(iBarrelLayer);
        TH1F* h_phi = new TH1F(histname.c_str(), histname.c_str(), nModulesBarrel, -TMath::Pi(), TMath::Pi());
        for (unsigned iBarrelModule = 0; iBarrelModule < nModulesBarrel; iBarrelModule++) {
          // set volumeID correctly
          dd4hep::DDSegmentation::CellID barrelVolumeId = 0;
          (*barrelDecoder)[m_fieldNamesSegmented[iSys]].set(barrelVolumeId, m_ecalBarrelSysId);
          (*barrelDecoder)["module"].set(barrelVolumeId, iBarrelModule);
          (*barrelDecoder)[m_activeFieldNamesSegmented[iSys]].set(barrelVolumeId, iBarrelLayer);

          numCells = det::utils::numberOfCells(barrelVolumeId, *moduleThetaSegmentation);

          (*barrelDecoder)["theta"].set(barrelVolumeId, numCells[2]);

          // for ECAL barrel with module readout, need to find out phi position of module 0 of ECal barrel last layer
          // get it from volume manager
          // for a normal phi grid it would suffice to use the segmentation phi(cellId) method
          dd4hep::VolumeManager volman = m_geoSvc->getDetector()->volumeManager();

          auto detelement = volman.lookupDetElement(barrelVolumeId);
          double inLocal[] = {0, 0, 0};
          const auto outGlobal = detelement.nominal().localToWorld(inLocal);
          double eCalBarrelPhiOffset = std::atan2(outGlobal.Y(), outGlobal.X());
          double eCalBarrelPhi = eCalBarrelPhiOffset;
          double eCalBarrelTheta = moduleThetaSegmentation->theta(barrelVolumeId);

          h_phi->SetBinContent(h_phi->FindBin(eCalBarrelPhi), iBarrelModule);
          phis.push_back(eCalBarrelPhi);

          if (iBarrelModule == 0)
            m_EMB_theta_lookup.push_back(eCalBarrelTheta);
        }
        m_EMB_phi_lookup.push_back(phis);
        m_EMB_h_module_vs_phi.push_back(h_phi);
      }
    }
  }

  debug() << "Lookup test " << m_EMEC_theta_lookup[9] << endmsg;
  return StatusCode::SUCCESS;
}

// StatusCode CreateFCCeeCaloNeighbours::finalize() { return Service::finalize(); }

StatusCode CreateFCCeeCaloNeighbours::finalize() {
    // for debug, to be removed
    /*
    TFile* f=new TFile("lookup.root", "RECREATE");
    for (unsigned int i=0; i<m_EMB_h_module_vs_phi.size(); i++) {
        m_EMB_h_module_vs_phi[i]->Write();
    }
    for (unsigned int i=0; i<m_EMEC_h_module_vs_phi_pos.size(); i++) {
        m_EMEC_h_module_vs_phi_pos[i]->Write();
    }
    for (unsigned int i=0; i<m_EMEC_h_module_vs_phi_neg.size(); i++) {
        m_EMEC_h_module_vs_phi_neg[i]->Write();
    }
    std::cout << "m_EMB_phi_lookup" << std::endl;
    for (unsigned int i=0; i<m_EMB_phi_lookup.size(); i++) {
        for (unsigned int j=0; j<m_EMB_phi_lookup[i].size(); j++) {
            std::cout << m_EMB_phi_lookup[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "m_EMB_theta_lookup" << std::endl;
    for (unsigned int i=0; i<m_EMB_theta_lookup.size(); i++) {
        std::cout << m_EMB_theta_lookup[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "m_EMEC_theta_lookup" << std::endl;
    for (unsigned int i=0; i<m_EMEC_theta_lookup.size(); i++) {
        std::cout << m_EMEC_theta_lookup[i] << " ";
    }
    std::cout << std::endl;
    */
    
    /*
    std::cout << "m_EMEC_pos_phi_lookup" << std::endl;
    for (unsigned int i=0; i<m_EMEC_pos_phi_lookup.size(); i++) {
        for (unsigned int j=0; j<m_EMEC_pos_phi_lookup[i].size(); j++) {
            std::cout << m_EMEC_pos_phi_lookup[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "m_EMEC_neg_phi_lookup" << std::endl;
    for (unsigned int i=0; i<m_EMEC_neg_phi_lookup.size(); i++) {
        for (unsigned int j=0; j<m_EMEC_neg_phi_lookup[i].size(); j++) {
            std::cout << m_EMEC_neg_phi_lookup[i][j] << " ";
        }
        std::cout << std::endl;
    }
    f->Close();
    */
    return Service::finalize();
}
