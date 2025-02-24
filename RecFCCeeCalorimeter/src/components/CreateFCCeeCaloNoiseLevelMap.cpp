#include "CreateFCCeeCaloNoiseLevelMap.h"

#include "DD4hep/Detector.h"
#include "detectorCommon/DetUtils_k4geo.h"
#include "k4Interface/IGeoSvc.h"
#include "detectorSegmentations/FCCSWGridModuleThetaMerged_k4geo.h"
#include "detectorSegmentations/FCCSWEndcapTurbine_k4geo.h"
#include "detectorSegmentations/FCCSWHCalPhiTheta_k4geo.h"
#include "detectorSegmentations/FCCSWHCalPhiRow_k4geo.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"

DECLARE_COMPONENT(CreateFCCeeCaloNoiseLevelMap)

CreateFCCeeCaloNoiseLevelMap::CreateFCCeeCaloNoiseLevelMap(const std::string &aName, ISvcLocator *aSL)
    : base_class(aName, aSL)
{
  declareProperty("ECalBarrelNoiseTool", m_ecalBarrelNoiseTool, "Handle for the cell noise tool of Barrel ECal");
  declareProperty("ECalEndcapNoiseTool", m_ecalEndcapNoiseTool, "Handle for the cell noise tool of Endcap ECal");
  declareProperty("HCalBarrelNoiseTool", m_hcalBarrelNoiseTool, "Handle for the cell noise tool of Barrel HCal");
  declareProperty("HCalEndcapNoiseTool", m_hcalEndcapNoiseTool, "Handle for the cell noise tool of Endcap HCal");
  declareProperty("outputFileName", m_outputFileName, "Name of the output file");
}

CreateFCCeeCaloNoiseLevelMap::~CreateFCCeeCaloNoiseLevelMap() {}

StatusCode CreateFCCeeCaloNoiseLevelMap::initialize()
{
  // Initialize necessary Gaudi components
  if (Service::initialize().isFailure())
  {
    return StatusCode::FAILURE;
  }
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc)
  {
    error() << "Unable to locate Geometry Service. "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }

  if (!m_ecalBarrelNoiseTool.empty())
  {
    if (!m_ecalBarrelNoiseTool.retrieve())
    {
      error() << "Unable to retrieve the ECAL barrel noise tool!!!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  if (!m_ecalEndcapNoiseTool.empty())
  {
    if (!m_ecalEndcapNoiseTool.retrieve())
    {
      error() << "Unable to retrieve the ECAL endcap noise tool!!!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  if (!m_hcalBarrelNoiseTool.empty())
  {
    if (!m_hcalBarrelNoiseTool.retrieve())
    {
      error() << "Unable to retrieve the HCAL barrel noise tool!!!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  if (!m_hcalEndcapNoiseTool.empty())
  {
    if (!m_hcalEndcapNoiseTool.retrieve())
    {
      error() << "Unable to retrieve the HCAL Endcap noise tool!!!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  std::unordered_map<uint64_t, std::pair<double, double>> map;

  for (uint iSys = 0; iSys < m_readoutNamesSegmented.size(); iSys++)
  {
    // Check if readouts exist
    info() << "Readout: " << m_readoutNamesSegmented[iSys] << endmsg;
    if (m_geoSvc->getDetector()->readouts().find(m_readoutNamesSegmented[iSys]) == m_geoSvc->getDetector()->readouts().end())
    {
      error() << "Readout <<" << m_readoutNamesSegmented[iSys] << ">> does not exist." << endmsg;
      return StatusCode::FAILURE;
    }

    // get segmentation
    dd4hep::DDSegmentation::Segmentation *aSegmentation = m_geoSvc->getDetector()->readout(m_readoutNamesSegmented[iSys]).segmentation().segmentation();
    std::string segmentationType = aSegmentation->type();
    info() << "Segmentation type : " << segmentationType << endmsg;
    
    dd4hep::DDSegmentation::GridTheta_k4geo *segmentation = nullptr;
    dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo *phiThetaSegmentation = nullptr;
    dd4hep::DDSegmentation::FCCSWHCalPhiTheta_k4geo *hcalPhiThetaSegmentation = nullptr;
    dd4hep::DDSegmentation::FCCSWHCalPhiRow_k4geo *hcalPhiRowSegmentation = nullptr;
    dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo *moduleThetaSegmentation = nullptr;
    dd4hep::DDSegmentation::FCCSWEndcapTurbine_k4geo *endcapTurbineSegmentation = nullptr;
    if (segmentationType == "FCCSWGridModuleThetaMerged_k4geo")
    {
      segmentation = dynamic_cast<dd4hep::DDSegmentation::GridTheta_k4geo *>(aSegmentation);
      moduleThetaSegmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo *>(aSegmentation);
    }
    else if (segmentationType == "FCCSWGridPhiTheta_k4geo")
    {
      segmentation = dynamic_cast<dd4hep::DDSegmentation::GridTheta_k4geo *>(aSegmentation);
      phiThetaSegmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo *>(aSegmentation);
    }
    else if (segmentationType == "FCCSWEndcapTurbine_k4geo")
    {
      endcapTurbineSegmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWEndcapTurbine_k4geo *>(aSegmentation);
    }
    else if (segmentationType == "FCCSWHCalPhiTheta_k4geo")
    {
      segmentation = dynamic_cast<dd4hep::DDSegmentation::GridTheta_k4geo *>(aSegmentation);
      hcalPhiThetaSegmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWHCalPhiTheta_k4geo *>(aSegmentation);
    }
    else if (segmentationType == "FCCSWHCalPhiRow_k4geo")
    {
      segmentation = dynamic_cast<dd4hep::DDSegmentation::GridTheta_k4geo *>(aSegmentation);
      hcalPhiRowSegmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWHCalPhiRow_k4geo *>(aSegmentation);
    }
    else
    {
      error() << "Segmentation type not handled." << endmsg;
      return StatusCode::FAILURE;
    }

    if(segmentation && segmentationType != "FCCSWHCalPhiRow_k4geo")
    {
      info() << "Segmentation: size in Theta " << segmentation->gridSizeTheta() << endmsg;
      info() << "Segmentation: offset in Theta " << segmentation->offsetTheta() << endmsg;
    }

    if (segmentationType == "FCCSWGridModuleThetaMerged_k4geo")
    {
      info() << "Segmentation: bins in Module " << moduleThetaSegmentation->nModules() << endmsg;
    }
    else if (segmentationType == "FCCSWGridPhiTheta_k4geo")
    {
      info() << "Segmentation: size in Phi " << phiThetaSegmentation->gridSizePhi() << endmsg;
      info() << "Segmentation: offset in Phi " << phiThetaSegmentation->offsetPhi() << endmsg;
    }
    else if (segmentationType == "FCCSWHCalPhiTheta_k4geo")
    {
      info() << "Segmentation: size in Phi " << hcalPhiThetaSegmentation->gridSizePhi() << endmsg;
      info() << "Segmentation: offset in Phi " << hcalPhiThetaSegmentation->offsetPhi() << endmsg;
    }

    // Loop over all cells in the calorimeter and retrieve existing cellIDs
    auto decoder = m_geoSvc->getDetector()->readout(m_readoutNamesSegmented[iSys]).idSpec().decoder();

    if (segmentationType == "FCCSWGridPhiTheta_k4geo" || segmentationType == "FCCSWGridModuleThetaMerged_k4geo")
    {
      // Loop over active layers
      std::vector<std::pair<int, int>> extrema;
      // extrema[0]: (0, nLayers-1)
      extrema.push_back(std::make_pair(0, m_activeVolumesNumbersSegmented[iSys] - 1));
      // extrema[1] and extrema[2] will contain theta and phi/module ranges, will be set later
      // as they (theta) depend on layer
      extrema.push_back(std::make_pair(0, 0));
      extrema.push_back(std::make_pair(0, 0));
      for (unsigned int ilayer = 0; ilayer < m_activeVolumesNumbersSegmented[iSys]; ilayer++)
      {
        if (segmentationType == "FCCSWGridPhiTheta_k4geo")
        {
          // for ECAL, G4 volumes extend along the full z so
          // they can be used - once obtained from the cellID - to
          // determine the full theta range of the layer.
          // For HCal the G4 volume containing the cell spans only
          // a subrange of the full z range of the HCal layer - one
          // would need to find the HCalLayerVol volume to determine
          // the theta range from it.
          // Thus for HCAL we use the theta range passed via activeVolumesTheta

          // Set system and layer of volume ID
          dd4hep::DDSegmentation::CellID volumeId = 0;
          (*decoder)[m_fieldNamesSegmented[iSys]].set(volumeId, m_fieldValuesSegmented[iSys]);
          (*decoder)[m_activeFieldNamesSegmented[iSys]].set(volumeId, ilayer);
          (*decoder)["theta"].set(volumeId, 0);
          (*decoder)["phi"].set(volumeId, 0);
          // Get number of segmentation cells within the active volume
          std::array<uint, 3> numCells;
          // for volumes (such as Allegro's HCAL) for which theta range
          // is passed via activeVolumesTheta option then use that instead of
          // numberOfCells helper function
          if (m_activeVolumesTheta[iSys].size()>0)
          {
            int nPhiCells = TMath::TwoPi() / phiThetaSegmentation->gridSizePhi();
            double thetaCellSize = phiThetaSegmentation->gridSizeTheta();
            double thetaOffset = phiThetaSegmentation->offsetTheta();
            double thetaMin = m_activeVolumesTheta[iSys][ilayer];
            double thetaMax = TMath::Pi() - thetaMin;
            uint minThetaID = int(floor((thetaMin + 0.5 * thetaCellSize - thetaOffset) / thetaCellSize));
            uint maxThetaID = int(floor((thetaMax + 0.5 * thetaCellSize - thetaOffset) / thetaCellSize));
            uint nThetaCells = 1 + maxThetaID - minThetaID;
            numCells[0] = nPhiCells;
            numCells[1] = nThetaCells;
            numCells[2] = minThetaID;
          }
          else
          {
            numCells = det::utils::numberOfCells(volumeId, *phiThetaSegmentation);
          }
          // extrema[1]: 0, ID of last cell in phi
          extrema[1] = std::make_pair(0, numCells[0] - 1);
          // extrema[2]: min theta ID, ID of last cell in theta
          extrema[2] = std::make_pair(numCells[2], numCells[1] + numCells[2] - 1);
          debug() << "Layer: " << ilayer << endmsg;
          debug() << "Extrema[0]: " << extrema[0].first << " , " << extrema[0].second << endmsg;
          debug() << "Extrema[1]: " << extrema[1].first << " , " << extrema[1].second << endmsg;
          debug() << "Extrema[2]: " << extrema[2].first << " , " << extrema[2].second << endmsg;
          debug() << "Number of segmentation cells in phi, in theta, and min theta ID, : " << numCells << endmsg;
          // Loop over segmentation cells
          for (unsigned int iphi = 0; iphi < numCells[0]; iphi++)
          {
            for (unsigned int itheta = 0; itheta < numCells[1]; itheta++)
            {
              dd4hep::DDSegmentation::CellID cellId = volumeId;
              decoder->set(cellId, "phi", iphi);
              decoder->set(cellId, "theta", itheta + numCells[2]); // start from the minimum existing theta cell in this layer
              uint64_t id = cellId;
              double noiseRMS = 0.;
              double noiseOffset = 0.;
              if (m_fieldValuesSegmented[iSys] == m_ecalBarrelSysId)
              {
                noiseRMS = m_ecalBarrelNoiseTool->getNoiseRMSPerCell(id);
                noiseOffset = m_ecalBarrelNoiseTool->getNoiseOffsetPerCell(id);
              }
              else if (m_fieldValuesSegmented[iSys] == m_hcalBarrelSysId)
              {
                noiseRMS = m_hcalBarrelNoiseTool->getNoiseRMSPerCell(id);
                noiseOffset = m_hcalBarrelNoiseTool->getNoiseOffsetPerCell(id);
              }
              else
              {
                warning() << "Unexpected system value for phi-theta readout " << m_fieldValuesSegmented[iSys]
                          << ", setting noise RMS and offset to 0.0" << endmsg;
              }
              map.insert(std::pair<uint64_t, std::pair<double, double>>(id, std::make_pair(noiseRMS, noiseOffset)));
            }
          }
        }
        else if (segmentationType == "FCCSWGridModuleThetaMerged_k4geo")
        {
          // Set system and layer of volume ID
          dd4hep::DDSegmentation::CellID volumeId = 0;
          (*decoder)[m_fieldNamesSegmented[iSys]].set(volumeId, m_fieldValuesSegmented[iSys]);
          (*decoder)[m_activeFieldNamesSegmented[iSys]].set(volumeId, ilayer);
          (*decoder)["theta"].set(volumeId, 0);
          (*decoder)["module"].set(volumeId, 0);
          // Get number of segmentation cells within the active volume, for given layer
          auto numCells = det::utils::numberOfCells(volumeId, *moduleThetaSegmentation);
          // extrema[1]: Range of module ID (0, max module ID)
          extrema[1] = std::make_pair(0, (numCells[0] - 1) * moduleThetaSegmentation->mergedModules(ilayer));
          // extrema[2]: Range of theta ID (min theta ID, max theta ID
          extrema[2] = std::make_pair(numCells[2], numCells[2] + (numCells[1] - 1) * moduleThetaSegmentation->mergedThetaCells(ilayer));
          debug() << "Layer: " << ilayer << endmsg;
          debug() << "Number of segmentation cells in (module, theta): " << numCells << endmsg;
          // Loop over segmentation cells
          for (unsigned int imodule = 0; imodule < numCells[0]; imodule++)
          {
            for (unsigned int itheta = 0; itheta < numCells[1]; itheta++)
            {
              dd4hep::DDSegmentation::CellID cellId = volumeId;
              decoder->set(cellId, "module", imodule * moduleThetaSegmentation->mergedModules(ilayer));
              decoder->set(cellId, "theta", numCells[2] + itheta * moduleThetaSegmentation->mergedThetaCells(ilayer)); // start from the minimum existing theta cell in this layer
              uint64_t id = cellId;
              double noiseRMS = 0.;
              double noiseOffset = 0.;
              if (m_fieldValuesSegmented[iSys] == m_ecalBarrelSysId)
              {
                noiseRMS = m_ecalBarrelNoiseTool->getNoiseRMSPerCell(id);
                noiseOffset = m_ecalBarrelNoiseTool->getNoiseOffsetPerCell(id);
              }
              else
              {
                warning() << "Unexpected system value for module-theta readout " << m_fieldValuesSegmented[iSys]
                          << ", setting noise RMS and offset to 0.0" << endmsg;
              }
              map.insert(std::pair<uint64_t, std::pair<double, double>>(id, std::make_pair(noiseRMS, noiseOffset)));
            }
          }
        }
      }
    }
    else if (segmentationType == "FCCSWEndcapTurbine_k4geo")
    {
      unsigned layerOffset[3];
      layerOffset[0] = 0;
      layerOffset[1] = endcapTurbineSegmentation->numCellsRhoCalib(0)*endcapTurbineSegmentation->numCellsZCalib(0);
      layerOffset[2] = layerOffset[1] + endcapTurbineSegmentation->numCellsRhoCalib(1)*endcapTurbineSegmentation->numCellsZCalib(1);
      // Loop over all cells in the calorimeter and retrieve existing cellIDs
      // Loop over active layers
      std::vector<std::pair<int, int>> extrema;
      extrema.push_back(std::make_pair(0, 2));  // wheel
      extrema.push_back(std::make_pair(0, 0)); // modules (set per wheel)
      extrema.push_back(std::make_pair(0, 0)); // rho (set per wheel)
      extrema.push_back(std::make_pair(0, 0)); // z (set per wheel)

      for (int iSide = -1; iSide < 2; iSide+=2) {
        for (unsigned int iWheel = 0; iWheel < 3; iWheel++) {
          dd4hep::DDSegmentation::CellID volumeId = 0;
          // Get VolumeID
          (*decoder)[m_fieldNamesSegmented[iSys]].set(volumeId, m_fieldValuesSegmented[iSys]);
          (*decoder)["side"].set(volumeId, iSide);
          (*decoder)["wheel"].set(volumeId, iWheel);
          (*decoder)["module"].set(volumeId, 0);
          (*decoder)["rho"].set(volumeId, 0);
          (*decoder)["z"].set(volumeId, 0);

          int numModules = endcapTurbineSegmentation->nModules(iWheel);
          int numCellsRho = endcapTurbineSegmentation->numCellsRho(iWheel);
          int numCellsZ = endcapTurbineSegmentation->numCellsZ(iWheel);
          int numCellsRhoCalib = endcapTurbineSegmentation->numCellsRhoCalib(iWheel);
          int numCellsZCalib = endcapTurbineSegmentation->numCellsZCalib(iWheel);

          // extrema 1: 0, ID of last module
          extrema[1] = std::make_pair(0, numModules - 1);
          // extrema[2]: 0, ID of last rho cell
          extrema[2] = std::make_pair(0, numCellsRho - 1);
          // extrema[3]: 0, ID of last z cell
          extrema[3] = std::make_pair(0, numCellsZ - 1);

          // Loop over segmentation cells
          for (int imodule = 0; imodule < numModules; imodule++) {
            for (int irho = 0; irho < numCellsRho; irho++) {
              for (int iz = 0; iz < numCellsZ; iz++) {
                dd4hep::DDSegmentation::CellID cellId = volumeId;
                decoder->set(cellId, "module", imodule);
                decoder->set(cellId, "rho", irho);
                decoder->set(cellId, "z", iz);
                unsigned iLayerZ = iz/(numCellsZ/numCellsZCalib);
                unsigned iLayerRho = irho/(numCellsRho/numCellsRhoCalib);
                unsigned iLayer = layerOffset[iWheel] + iLayerRho*numCellsZCalib + iLayerZ;
                decoder->set(cellId, "layer", iLayer);
                uint64_t id = cellId;
                if (iSide==1 && iWheel == 2 && imodule == 113 && irho == 19 && iz == 1) {
                  debug() << "in test cell, iLayer = " << iLayer << " and cell ID = " << id << endmsg;
                }
                double noiseRMS = m_ecalEndcapNoiseTool->getNoiseRMSPerCell(id);
                double noiseOffset = m_ecalEndcapNoiseTool->getNoiseOffsetPerCell(id);
                map.insert( std::pair<uint64_t, std::pair<double, double> >(id, std::make_pair(noiseRMS, noiseOffset) ) );
              } // end loop over z
            } // end loop over rho
          } // end loop over module
        } // end loop over wheel
      } // end looper over side
    } // end if (segmentationType == "FCCSWEndcapTurbine_k4geo")

    if(segmentationType == "FCCSWHCalPhiTheta_k4geo" || segmentationType == "FCCSWHCalPhiRow_k4geo")
    {
      for (unsigned int ilayer = 0; ilayer < m_activeVolumesNumbersSegmented[iSys]; ilayer++)
      {
        // HCal phi-theta segmentation
        if (segmentationType == "FCCSWHCalPhiTheta_k4geo")
        {
          // Set system and layer of volume ID
          dd4hep::DDSegmentation::CellID volumeId = 0;
          (*decoder)[m_fieldNamesSegmented[iSys]].set(volumeId, m_fieldValuesSegmented[iSys]);
          (*decoder)[m_activeFieldNamesSegmented[iSys]].set(volumeId, ilayer);
          (*decoder)["theta"].set(volumeId, 0);
          (*decoder)["phi"].set(volumeId, 0);

          // Get number of segmentation cells within the active volume
          std::array<uint, 3> numCells;
          std::vector<int> thetaBins(hcalPhiThetaSegmentation->thetaBins(ilayer));
          numCells[0] = hcalPhiThetaSegmentation->phiBins();

          // Barrel
          if(m_fieldNamesSegmented[iSys] == "system" && m_fieldValuesSegmented[iSys] == m_hcalBarrelSysId)
          {
            // number of cells in the layer
            numCells[1] = thetaBins.size();
            // ID of first cell theta bin (smallest theta)
            if(numCells[1] > 0) numCells[2] = thetaBins[0];
            else
            {
              error() << "Can not get number of cells from the segmentation object: " << segmentationType << endmsg;
              return StatusCode::FAILURE;
            }
            debug() << "Layer: " << ilayer << endmsg;
            debug() << "Number of segmentation cells in phi, in theta, and min theta ID, for HCal barrel : " << numCells << endmsg;
            // Loop over segmentation cells
            for (unsigned int iphi = 0; iphi < numCells[0]; iphi++)
            {
              for (unsigned int itheta = 0; itheta < numCells[1]; itheta++)
              {
                dd4hep::DDSegmentation::CellID cellId = volumeId;
                decoder->set(cellId, "phi", iphi);
                decoder->set(cellId, "theta", itheta + numCells[2]); // start from the minimum existing theta cell in this layer
                uint64_t id = cellId;
                double noiseRMS = m_hcalBarrelNoiseTool->getNoiseRMSPerCell(id);
                double noiseOffset = m_hcalBarrelNoiseTool->getNoiseOffsetPerCell(id);
                map.insert(std::pair<uint64_t, std::pair<double, double>>(id, std::make_pair(noiseRMS, noiseOffset)));
              }
            }
          } // barrel

          // For endcap, determine noise separately for positive- and negative-z part
          if(m_fieldNamesSegmented[iSys] == "system" && m_fieldValuesSegmented[iSys] == m_hcalEndcapSysId)
          {
            // number of cells in the layer
            numCells[1] = thetaBins.size()/2; // thetaBins contains theta IDs for both positive- and negative-z parts, hence divide size by 2

            // ID of first cell theta bin (smallest theta) in the positive-z part of the endcap
            if(numCells[1] > 0) numCells[2] = thetaBins[0];
            else
            {
              error() << "Can not get number of cells from the segmentation object: " << segmentationType << endmsg;
              return StatusCode::FAILURE;
            }

            debug() << "Layer: " << ilayer << endmsg;
            debug() << "Number of segmentation cells in phi, in theta, and min theta ID, for positive-z HCal endcap : " << numCells << endmsg;
            // Loop over segmentation cells
            for (unsigned int iphi = 0; iphi < numCells[0]; iphi++)
            {
              for (unsigned int itheta = 0; itheta < numCells[1]; itheta++)
              {
                dd4hep::DDSegmentation::CellID cellId = volumeId;
                decoder->set(cellId, "phi", iphi);
                decoder->set(cellId, "theta", itheta + numCells[2]); // start from the minimum existing theta cell in this layer
                uint64_t id = cellId;
                double noiseRMS = m_hcalEndcapNoiseTool->getNoiseRMSPerCell(id);
                double noiseOffset = m_hcalEndcapNoiseTool->getNoiseOffsetPerCell(id);
                map.insert(std::pair<uint64_t, std::pair<double, double>>(id, std::make_pair(noiseRMS, noiseOffset)));
              }
            }

            // ID of first cell theta bin (smallest theta) in the negative-z part of the endcap
            numCells[2] = thetaBins[thetaBins.size()/2];

            debug() << "Layer: " << ilayer << endmsg;
            debug() << "Number of segmentation cells in phi, in theta, and min theta ID, for negative-z HCal endcap : " << numCells << endmsg;
            // Loop over segmentation cells
            for (unsigned int iphi = 0; iphi < numCells[0]; iphi++)
            {
              for (unsigned int itheta = 0; itheta < numCells[1]; itheta++)
              {
                dd4hep::DDSegmentation::CellID cellId = volumeId;
                decoder->set(cellId, "phi", iphi);
                decoder->set(cellId, "theta", itheta + numCells[2]); // start from the minimum existing theta cell in this layer
                uint64_t id = cellId;
                double noiseRMS = m_hcalEndcapNoiseTool->getNoiseRMSPerCell(id);
                double noiseOffset = m_hcalEndcapNoiseTool->getNoiseOffsetPerCell(id);
                map.insert(std::pair<uint64_t, std::pair<double, double>>(id, std::make_pair(noiseRMS, noiseOffset)));
              }
            }
          } // Endcap
        }
        // HCal phi-row segmentation
        else if (segmentationType == "FCCSWHCalPhiRow_k4geo")
        {
          // Set system and layer of volume ID
          dd4hep::DDSegmentation::CellID volumeId = 0;
          (*decoder)[m_fieldNamesSegmented[iSys]].set(volumeId, m_fieldValuesSegmented[iSys]);
          (*decoder)[m_activeFieldNamesSegmented[iSys]].set(volumeId, ilayer);
          (*decoder)["row"].set(volumeId, 0);
          (*decoder)["phi"].set(volumeId, 0);

          // Get number of segmentation cells within the active volume
          std::array<int, 3> numCells;
          std::vector<int> cellIndexes(hcalPhiRowSegmentation->cellIndexes(ilayer));
          numCells[0] = hcalPhiRowSegmentation->phiBins();

          // Barrel
          if(m_fieldNamesSegmented[iSys] == "system" && m_fieldValuesSegmented[iSys] == m_hcalBarrelSysId)
          {
            // number of cells in the layer
            numCells[1] = cellIndexes.size();
            // minimum cell index
            if(numCells[1] > 0) numCells[2] = cellIndexes.front();
            else
            {
              error() << "Can not get number of cells from the segmentation object: " << segmentationType << endmsg;
              return StatusCode::FAILURE;
            }

            debug() << "Layer: " << ilayer << endmsg;
            debug() << "Number of segmentation cells in phi, in row, and min cell index, for HCal barrel : " << numCells << endmsg;
            // Loop over segmentation cells
            for (int iphi = 0; iphi < numCells[0]; iphi++)
            {
              for (int irow = 0; irow < numCells[1]; irow++)
              {
                dd4hep::DDSegmentation::CellID cellId = volumeId;
                decoder->set(cellId, "phi", iphi);
                decoder->set(cellId, "row", irow + numCells[2]); // start from the minimum existing cell index in this layer
                uint64_t id = cellId;
                double noiseRMS = m_hcalBarrelNoiseTool->getNoiseRMSPerCell(id);
                double noiseOffset = m_hcalBarrelNoiseTool->getNoiseOffsetPerCell(id);
                map.insert(std::pair<uint64_t, std::pair<double, double>>(id, std::make_pair(noiseRMS, noiseOffset)));
              }
            }
          } // barrel

          // For endcap, determine noise separately for positive- and negative-z part
          if(m_fieldNamesSegmented[iSys] == "system" && m_fieldValuesSegmented[iSys] == m_hcalEndcapSysId)
          {
            // number of cells in the layer
            numCells[1] = cellIndexes.size()/2; // cellIndexes contains cell IDs for both positive- and negative-z parts, hence divide size by 2

            // minimum cell index in the positive-z part of the endcap
            if(numCells[1] > 0) numCells[2] = cellIndexes.front();
            else
            {
              error() << "Can not get number of cells from the segmentation object: " << segmentationType << endmsg;
              return StatusCode::FAILURE;
            }
            debug() << "Layer: " << ilayer << endmsg;
            debug() << "Number of segmentation cells in phi, in row, and min cell index, for positive-z HCal endcap : " << numCells << endmsg;
            // Loop over segmentation cells
            for (int iphi = 0; iphi < numCells[0]; iphi++)
            {
              for (int irow = 0; irow < numCells[1]; irow++)
              {
                dd4hep::DDSegmentation::CellID cellId = volumeId;
                decoder->set(cellId, "phi", iphi);
                decoder->set(cellId, "row", irow + numCells[2]); // start from the minimum existing cell index in this layer
                uint64_t id = cellId;
                double noiseRMS = m_hcalEndcapNoiseTool->getNoiseRMSPerCell(id);
                double noiseOffset = m_hcalEndcapNoiseTool->getNoiseOffsetPerCell(id);
                map.insert(std::pair<uint64_t, std::pair<double, double>>(id, std::make_pair(noiseRMS, noiseOffset)));
              }
            }

            // minimum cell index in the negative-z part of the endcap
            numCells[2] = cellIndexes.back();

            debug() << "Layer: " << ilayer << endmsg;
            debug() << "Number of segmentation cells in phi, in row, and min cell index, for negative-z HCal endcap : " << numCells << endmsg;
            // Loop over segmentation cells
            for (int iphi = 0; iphi < numCells[0]; iphi++)
            {
              for (int irow = 0; irow < numCells[1]; irow++)
              {
                dd4hep::DDSegmentation::CellID cellId = volumeId;
                decoder->set(cellId, "phi", iphi);
                decoder->set(cellId, "row", irow + numCells[2]); // start from the minimum cell index in this layer
                uint64_t id = cellId;
                double noiseRMS = m_hcalEndcapNoiseTool->getNoiseRMSPerCell(id);
                double noiseOffset = m_hcalEndcapNoiseTool->getNoiseOffsetPerCell(id);
                map.insert(std::pair<uint64_t, std::pair<double, double>>(id, std::make_pair(noiseRMS, noiseOffset)));
              }
            }
          } // Endcap
        } // hcal phi-row segmentation
      } // layer loop
    } // hcal segmentations
  }

  // Check if output directory exists
  std::string outDirPath = gSystem->DirName(m_outputFileName.c_str());
  if (!gSystem->OpenDirectory(outDirPath.c_str())) {
    error() << "Output directory \"" << outDirPath
            << "\" does not exists! Please create it." << endmsg;
    return StatusCode::FAILURE;
  }

  std::unique_ptr<TFile> outFile(TFile::Open(m_outputFileName.c_str(), "RECREATE"));
  outFile->cd();
  TTree tree("noisyCells", "Tree with map of noise per cell");
  uint64_t saveCellId;
  double saveNoiseRMS;
  double saveNoiseOffset;
  tree.Branch("cellId", &saveCellId, "cellId/l");
  tree.Branch("noiseLevel", &saveNoiseRMS); // would be better to call it noiseRMS
  tree.Branch("noiseOffset", &saveNoiseOffset);
  for (const auto &item : map)
  {
    saveCellId = item.first;
    saveNoiseRMS = item.second.first;
    saveNoiseOffset = item.second.second;
    tree.Fill();
  }
  outFile->Write();
  outFile->Close();

  return StatusCode::SUCCESS;
}

StatusCode CreateFCCeeCaloNoiseLevelMap::finalize() { return Service::finalize(); }
