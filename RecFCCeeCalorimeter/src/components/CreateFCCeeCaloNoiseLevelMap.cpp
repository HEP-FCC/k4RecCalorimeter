#include "CreateFCCeeCaloNoiseLevelMap.h"

#include "DD4hep/Detector.h"
#include "detectorCommon/DetUtils_k4geo.h"
#include "k4Interface/IGeoSvc.h"
#include "detectorSegmentations/FCCSWGridModuleThetaMerged_k4geo.h"

#include "TFile.h"
#include "TTree.h"

DECLARE_COMPONENT(CreateFCCeeCaloNoiseLevelMap)

CreateFCCeeCaloNoiseLevelMap::CreateFCCeeCaloNoiseLevelMap(const std::string &aName, ISvcLocator *aSL)
    : base_class(aName, aSL)
{
  declareProperty("ECalBarrelNoiseTool", m_ecalBarrelNoiseTool, "Handle for the cell noise tool of Barrel ECal");
  declareProperty("HCalBarrelNoiseTool", m_hcalBarrelNoiseTool, "Handle for the cell noise tool of Barrel HCal");
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

  if (!m_ecalBarrelNoiseTool.retrieve())
  {
    error() << "Unable to retrieve the ECAL noise tool!!!" << endmsg;
    return StatusCode::FAILURE;
  }

  if (!m_hcalBarrelNoiseTool.empty())
  {
    if (!m_hcalBarrelNoiseTool.retrieve())
    {
      error() << "Unable to retrieve the HCAL noise tool!!!" << endmsg;
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
    dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo *moduleThetaSegmentation = nullptr;
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
    else
    {
      error() << "Segmentation type not handled." << endmsg;
      return StatusCode::FAILURE;
    }
    
    if (segmentation == nullptr || (phiThetaSegmentation == nullptr && moduleThetaSegmentation == nullptr))
    {
      error() << "Unable to cast segmentation pointer!!!!" << endmsg;
      return StatusCode::FAILURE;
    }

    info() << "Segmentation: size in Theta " << segmentation->gridSizeTheta() << endmsg;
    info() << "Segmentation: offset in Theta " << segmentation->offsetTheta() << endmsg;
    if (segmentationType == "FCCSWGridModuleThetaMerged_k4geo")
    {
      info() << "Segmentation: bins in Module " << moduleThetaSegmentation->nModules() << endmsg;
    }
    else if (segmentationType == "FCCSWGridPhiTheta_k4geo")
    {
      info() << "Segmentation: size in Phi " << phiThetaSegmentation->gridSizePhi() << endmsg;
      info() << "Segmentation: offset in Phi " << phiThetaSegmentation->offsetPhi() << endmsg;
    }

    // Loop over all cells in the calorimeter and retrieve existing cellIDs
    // Loop over active layers
    auto decoder = m_geoSvc->getDetector()->readout(m_readoutNamesSegmented[iSys]).idSpec().decoder();
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
            double noise = 0.;
            double noiseOffset = 0.;
            if (m_fieldValuesSegmented[iSys] == m_ecalBarrelSysId)
            {
              noise = m_ecalBarrelNoiseTool->getNoiseConstantPerCell(id);
              noiseOffset = m_ecalBarrelNoiseTool->getNoiseOffsetPerCell(id);
            }
            else if (m_fieldValuesSegmented[iSys] == m_hcalBarrelSysId)
            {
              noise = m_hcalBarrelNoiseTool->getNoiseConstantPerCell(id);
              noiseOffset = m_hcalBarrelNoiseTool->getNoiseOffsetPerCell(id);
            }
            else
            {
              warning() << "Unexpected system value for phi-theta readout " << m_fieldValuesSegmented[iSys]
                        << ", setting noise to 0.0" << endmsg;
            }
            map.insert(std::pair<uint64_t, std::pair<double, double>>(id, std::make_pair(noise, noiseOffset)));
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
            double noise = 0.;
            double noiseOffset = 0.;
            if (m_fieldValuesSegmented[iSys] == m_ecalBarrelSysId)
            {
              noise = m_ecalBarrelNoiseTool->getNoiseConstantPerCell(id);
              noiseOffset = m_ecalBarrelNoiseTool->getNoiseOffsetPerCell(id);
            }
            else
            {
              warning() << "Unexpected system value for module-theta readout " << m_fieldValuesSegmented[iSys]
                        << ", setting noise to 0.0" << endmsg;
            }
            map.insert(std::pair<uint64_t, std::pair<double, double>>(id, std::make_pair(noise, noiseOffset)));
          }
        }
      }
    }
  }

  std::unique_ptr<TFile> file(TFile::Open(m_outputFileName.c_str(), "RECREATE"));
  file->cd();
  TTree tree("noisyCells", "Tree with map of noise per cell");
  uint64_t saveCellId;
  double saveNoiseLevel;
  double saveNoiseOffset;
  tree.Branch("cellId", &saveCellId, "cellId/l");
  tree.Branch("noiseLevel", &saveNoiseLevel);
  tree.Branch("noiseOffset", &saveNoiseOffset);
  for (const auto &item : map)
  {
    saveCellId = item.first;
    saveNoiseLevel = item.second.first;
    saveNoiseOffset = item.second.second;
    tree.Fill();
  }
  file->Write();
  file->Close();

  return StatusCode::SUCCESS;
}

StatusCode CreateFCCeeCaloNoiseLevelMap::finalize() { return Service::finalize(); }
