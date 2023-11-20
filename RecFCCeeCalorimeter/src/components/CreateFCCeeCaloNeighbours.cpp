#include "CreateFCCeeCaloNeighbours.h"

// DD4hep
#include "DD4hep/Detector.h"

// k4FWCore
#include "k4Interface/IGeoSvc.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"
#include "detectorSegmentations/GridTheta_k4geo.h"
#include "detectorSegmentations/FCCSWGridPhiTheta_k4geo.h"
#include "detectorSegmentations/FCCSWGridModuleThetaMerged_k4geo.h"

// ROOT
#include "TFile.h"
#include "TTree.h"

DECLARE_COMPONENT(CreateFCCeeCaloNeighbours)

CreateFCCeeCaloNeighbours::CreateFCCeeCaloNeighbours(const std::string &aName, ISvcLocator *aSL)
    : base_class(aName, aSL)
{
  declareProperty("outputFileName", m_outputFileName, "Name of the output file");
}

CreateFCCeeCaloNeighbours::~CreateFCCeeCaloNeighbours() {}

StatusCode CreateFCCeeCaloNeighbours::initialize()
{
  // Initialize necessary Gaudi components
  if (Service::initialize().isFailure())
  {
    error() << "Unable to initialize Service()" << endmsg;
    return StatusCode::FAILURE;
  }
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc)
  {
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
  // double eCalThetaOffset = 0;
  // double eCalThetaSize = 0;
  // double eCalPhiOffset = 0;
  // double eCalModuleSize = 0;
  // double hCalThetaOffset = 0;
  // double hCalThetaSize = 0;
  // double hCalPhiOffset = 0;
  dd4hep::DDSegmentation::BitFieldCoder *decoderECalBarrel = nullptr;
  //  will be used for volume connecting
  std::pair<int, int> extremaHCalFirstLayerPhi;
  std::pair<int, int> extremaHCalFirstLayerTheta;
  std::pair<int, int> extremaHCalFirstLayerZ;
  dd4hep::DDSegmentation::BitFieldCoder *decoderHCalBarrel = nullptr;

  for (uint iSys = 0; iSys < m_readoutNamesSegmented.size(); iSys++)
  {
    // Check if readouts exist
    info() << "Readout: " << m_readoutNamesSegmented[iSys] << endmsg;
    if (m_geoSvc->getDetector()->readouts().find(m_readoutNamesSegmented[iSys]) == m_geoSvc->getDetector()->readouts().end())
    {
      error() << "Readout <<" << m_readoutNamesSegmented[iSys] << ">> does not exist." << endmsg;
      return StatusCode::FAILURE;
    }

    // get Theta-based segmentation
    dd4hep::DDSegmentation::Segmentation *aSegmentation = m_geoSvc->getDetector()->readout(m_readoutNamesSegmented[iSys]).segmentation().segmentation();
    if (aSegmentation == nullptr)
    {
      error() << "Segmentation does not exist." << endmsg;
      return StatusCode::FAILURE;
    }

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

    // retrieve decoders and other info needed for volume (ECal-HCal) connection
    auto decoder = m_geoSvc->getDetector()->readout(m_readoutNamesSegmented[iSys]).idSpec().decoder();
    // hardcoding values (4,8) for ECal and HCal barrel systems looks quite error prone..
    // to be fixed in the future
    if (m_fieldNamesSegmented[iSys] == "system" && m_fieldValuesSegmented[iSys] == 4)
    {
      decoderECalBarrel = decoder;
      // eCalThetaSize = segmentation->gridSizeTheta();
      // eCalModuleSize = 2 * M_PI / segmentation->nModules();
      // eCalModuleSize = 2 * M_PI / segmentation->phiBins();
      // eCalThetaOffset = segmentation->offsetTheta();
      // eCalPhiOffset = segmentation->offsetPhi();
    }
    if (m_fieldNamesSegmented[iSys] == "system" && m_fieldValuesSegmented[iSys] == 8)
    {
      decoderHCalBarrel = decoder;
      // hCalThetaSize = segmentation->gridSizeTheta();
      // hCalThetaOffset = segmentation->offsetTheta();
      // hCalPhiOffset = segmentation->offsetPhi();
    }

    if (segmentationType == "FCCSWGridPhiTheta_k4geo")
    {
      // GM: code copied from RecFCChhCalorimeter, just replacing Eta->Theta
      //     did not review the code itself

      // Loop over all cells in the calorimeter and retrieve existing cellIDs
      // Loop over active layers
      std::vector<std::pair<int, int>> extrema;
      extrema.push_back(std::make_pair(0, m_activeVolumesNumbersSegmented[iSys] - 1));
      extrema.push_back(std::make_pair(0, 0));
      extrema.push_back(std::make_pair(0, 0));
      for (unsigned int ilayer = 0; ilayer < m_activeVolumesNumbersSegmented[iSys]; ilayer++)
      {
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
        // for HCAL use theta range passed via activeVolumesTheta option
        if (m_fieldValuesSegmented[iSys]==8) {
          int nPhiCells = TMath::TwoPi()/phiThetaSegmentation->gridSizePhi();
          double thetaCellSize = phiThetaSegmentation->gridSizeTheta();
          double thetaOffset = phiThetaSegmentation->offsetTheta();
          double thetaMin = m_activeVolumesTheta[ilayer];
          double thetaMax = TMath::Pi() - thetaMin;
          // debug
          std::cout << "thetaMin, thetaMax = " << thetaMin << " " << thetaMax << std::endl;
          uint minThetaID = int(floor((thetaMin + 0.5 * thetaCellSize - thetaOffset) / thetaCellSize));
          uint maxThetaID = int(floor((thetaMax + 0.5 * thetaCellSize - thetaOffset) / thetaCellSize));
          uint nThetaCells = 1 + maxThetaID - minThetaID;
          numCells[0] = nPhiCells;
          numCells[1] = nThetaCells;
          numCells[2] = minThetaID;
        }
        else {
          numCells = det::utils::numberOfCells(volumeId, *phiThetaSegmentation);
        }
        // extrema 1: 0, ID of last cell in phi
        extrema[1] = std::make_pair(0, numCells[0] - 1);
        // extrema[2]: min theta ID, n theta cells     
        extrema[2] = std::make_pair(numCells[2], numCells[1] + numCells[2] - 1);

        // for layer N-1 of ECal barrel,  will be used for volume connecting
        /*
        if (
            m_fieldNamesSegmented[iSys] == "system" &&
            m_fieldValuesSegmented[iSys] == 4 &&
            ilayer == (m_activeVolumesNumbersSegmented[iSys] - 1)
        )
        {
          eCalLastLayer = m_activeVolumesNumbersSegmented[iSys] - 1;
          // not really needed, same for all layers... unless we really restrict theta for each layer 
          // to the physical volume (not sure it's really needed, we can maybe allow for non-physical 
          // cells at edge in neighbour map, they simply won´t be added to the cluster since they do
          // not exist in the list of cells)
          extremaECalLastLayerPhi = std::make_pair(0, numCells[0] - 1);
          extremaECalLastLayerTheta = std::make_pair(numCells[2], numCells[1] + numCells[2] - 1);
        }
        else if (
          m_fieldNamesSegmented[iSys] == "system" &&
          m_fieldValuesSegmented[iSys] == 8 && 
          m_readoutNamesSegmented[iSys] == "HCalBarrelReadout" // why is this needed? because there's also the HCAL endcap?
        )
        {
          uint cellsTheta = ceil((2 * m_activeVolumesTheta[ilayer] - segmentation->gridSizeTheta()) / 2 / segmentation->gridSizeTheta()) * 2 + 1; // ceil( 2*m_activeVolumesRadii[ilayer] / segmentation->gridSizeEta()) ;
          uint minThetaID = int(floor((-m_activeVolumesTheta[ilayer] + 0.5 * segmentation->gridSizeTheta() - segmentation->offsetTheta()) / segmentation->gridSizeTheta()));
          numCells[1] = cellsTheta;
          numCells[2] = minThetaID;
          // for layer 0 of HCal barrel,  will be used for volume connecting
          if (ilayer == 0)
          {
            extremaHCalFirstLayerPhi = std::make_pair(0, numCells[0] - 1);
            extremaHCalFirstLayerTheta = std::make_pair(numCells[2], numCells[1] + numCells[2] - 1);
            extrema[2] = std::make_pair(numCells[2], numCells[1] + numCells[2] - 1);
          }
        }
        */
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
	    if (ilayer==7 && iphi==95) debug() << itheta+numCells[2] << " " << cellId << endmsg;
            map.insert(std::pair<uint64_t, std::vector<uint64_t>>(
                id, det::utils::neighbours(*decoder, {m_activeFieldNamesSegmented[iSys], "phi", "theta"}, extrema,
                                           id, {false, true, false}, m_includeDiagonalCells)));
          }
        }
      }
    }

    else if (segmentationType == "FCCSWGridModuleThetaMerged_k4geo")
    {
      std::vector<std::pair<int, int>> extrema;
      // extrema[0]: min layer, n layers
      extrema.push_back(std::make_pair(0, m_activeVolumesNumbersSegmented[iSys] - 1));
      extrema.push_back(std::make_pair(0, 0));
      extrema.push_back(std::make_pair(0, 0));
      for (unsigned int ilayer = 0; ilayer < m_activeVolumesNumbersSegmented[iSys]; ilayer++)
      {
        dd4hep::DDSegmentation::CellID volumeId = 0;
        // Get VolumeID
        // m_fieldValuesSegmented: in .py systemValuesModuleTheta = [4]
        // m_activeFieldNamesSegmented: in .py activeFieldNamesModuleTheta = ["layer"]
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
        extrema[2] = std::make_pair(numCells[2], numCells[2] + (numCells[1] - 1) * moduleThetaSegmentation->mergedThetaCells(ilayer));
        // for layer N-1 of ECal barrel,  will be used for volume connecting
        // should 4 be systemValuesModuleTheta instead?
        if (ilayer == (m_activeVolumesNumbersSegmented[iSys] - 1) && m_fieldNamesSegmented[iSys] == "system" &&
            m_fieldValuesSegmented[iSys] == 4)
        {
          // eCalLastLayer = m_activeVolumesNumbersSegmented[iSys] - 1;
          extremaECalLastLayerModule = std::make_pair(0, numCells[0] - 1);
          extremaECalLastLayerTheta = std::make_pair(numCells[2], numCells[1] + numCells[2] - 1);
        }
        else if (m_fieldNamesSegmented[iSys] == "system" &&
                 m_fieldValuesSegmented[iSys] == 8 && m_readoutNamesSegmented[iSys] == "BarHCal_Readout_phitheta")
        {
          ////
          uint cellsTheta = ceil((2 * m_activeVolumesTheta[ilayer] - segmentation->gridSizeTheta()) / 2 / segmentation->gridSizeTheta()) * 2 + 1; // ceil( 2*m_activeVolumesRadii[ilayer] / segmentation->gridSizeTheta());
          uint minThetaID = int(floor((-m_activeVolumesTheta[ilayer] + 0.5 * segmentation->gridSizeTheta() - segmentation->offsetTheta()) / segmentation->gridSizeTheta()));
          numCells[1] = cellsTheta;
          numCells[2] = minThetaID;
          // for layer 0 of HCal barrel,  will be used for volume connecting
          if (ilayer == 0)
          {
            extremaHCalFirstLayerPhi = std::make_pair(0, numCells[0] - 1);
            extremaHCalFirstLayerTheta = std::make_pair(numCells[2], numCells[1] + numCells[2] - 1);
          }
        }
        debug() << "Layer: " << ilayer << endmsg;
        debug() << "Extrema[0]: " << extrema[0].first << " , " << extrema[0].second << endmsg;
        debug() << "Extrema[1]: " << extrema[1].first << " , " << extrema[1].second << endmsg;
        debug() << "Extrema[2]: " << extrema[2].first << " , " << extrema[2].second << endmsg;
        debug() << "Number of segmentation cells in (module,theta): " << numCells << endmsg;
        // Loop over segmentation cells
        for (int imodule = extrema[1].first; imodule <= extrema[1].second; imodule += moduleThetaSegmentation->mergedModules(ilayer))
        {
          for (int itheta = extrema[2].first; itheta <= extrema[2].second; itheta += moduleThetaSegmentation->mergedThetaCells(ilayer))
          {
            dd4hep::DDSegmentation::CellID cellId = volumeId;
            decoder->set(cellId, "module", imodule);
            decoder->set(cellId, "theta", itheta); // start from the minimum existing theta cell in this layer
            uint64_t id = cellId;
            map.insert(std::pair<uint64_t, std::vector<uint64_t>>(
                id,
                det::utils::neighbours_ModuleThetaMerged(
                    *moduleThetaSegmentation,
                    *decoder,
                    {m_activeFieldNamesSegmented[iSys],
                     "module", "theta"},
                    extrema,
                    id,
                    m_includeDiagonalCells)));
          }
        }
      }
    }

    if (msgLevel() <= MSG::DEBUG)
    {
      std::vector<int> counter;
      counter.assign(40, 0);
      for (const auto &item : map)
      {
        counter[item.second.size()]++;
      }
      for (uint iCount = 0; iCount < counter.size(); iCount++)
      {
        if (counter[iCount] != 0)
        {
          info() << counter[iCount] << " cells have " << iCount << " neighbours" << endmsg;
        }
      }
    }
    info() << "total number of cells:  " << map.size() << endmsg;
  }

  //////////////////////////////////////////////////
  ///      BARREL: connection ECAL + HCAL        ///
  /////////////////////////////////////////////////
  int count = 0;

  if (m_connectBarrels)
  {
    /*
    // first check if ECAL barrel (system==5) and HCal barrel (system==8) are configured
    if (decoderECalBarrel == nullptr || decoderHCalBarrel == nullptr)
    {
      error() << "ECAL barrel and/or HCal barrel are not configured correctly!" << endmsg;
      return StatusCode::FAILURE;
    }
    // print how many cells in each dimensions will be matched
    info() << "ECAL layer " << eCalLastLayer << " is a neighbour of HCAL layer 0." << endmsg;
    info() << "ECAL phi cells " << extremaECalLastLayerModule.first << " - " << extremaECalLastLayerModule.second
           << " will be matched to HCAL cells " << extremaHCalFirstLayerPhi.first
           << " - " << extremaHCalFirstLayerPhi.second << endmsg;
    info() << "ECAL theta cells " << extremaECalLastLayerTheta.first << " - " << extremaECalLastLayerTheta.second
           << " will be matched to HCAL " << extremaHCalFirstLayerTheta.first
           << " - " << extremaHCalFirstLayerTheta.second << endmsg;

    std::unordered_map<uint, std::vector<uint>> thetaNeighbours;
    std::unordered_map<uint, std::vector<uint>> phiNeighbours;
    double hCalPhiSize = 2 * M_PI / (extremaHCalFirstLayerPhi.second - extremaHCalFirstLayerPhi.first + 1);
    // loop over theta cells to match in theta
    for (int iTheta = extremaHCalFirstLayerTheta.first; iTheta < extremaHCalFirstLayerTheta.second + 1; iTheta++)
    {
      double lowTheta = hCalThetaOffset + iTheta * hCalThetaSize;
      double highTheta = hCalThetaOffset + (iTheta + 1) * hCalThetaSize;
      debug() << "HCal theta range  : " << lowTheta << " -  " << highTheta << endmsg;
      int lowId = floor((lowTheta - 0.5 * eCalThetaSize - eCalThetaOffset) / eCalThetaSize);
      int highId = floor((highTheta + 0.5 * eCalThetaSize - eCalThetaOffset) / eCalThetaSize);
      debug() << "ECal theta range  : " << lowId * eCalThetaSize + eCalThetaOffset << " -  "
              << highId * eCalThetaSize + eCalThetaOffset << endmsg;
      std::vector<uint> neighbours;
      for (int idThetaToAdd = lowId; idThetaToAdd <= highId; idThetaToAdd++)
      {
        neighbours.push_back(det::utils::cyclicNeighbour(idThetaToAdd, extremaECalLastLayerTheta));
      }
      debug() << "HCal theta id  : " << iTheta << endmsg;
      debug() << "Found ECal Neighbours in theta : " << neighbours.size() << endmsg;
      for (auto id : neighbours)
      {
        debug() << "ECal Neighbours id : " << id << endmsg;
      }
      thetaNeighbours.insert(std::pair<uint, std::vector<uint>>(iTheta, neighbours));
    }

    // loop over phi and find which phi cells to add
    for (int iPhi = 0; iPhi < extremaHCalFirstLayerPhi.second + 1; iPhi++)
    {
      double lowPhi = hCalPhiOffset + iPhi * hCalPhiSize;
      double highPhi = hCalPhiOffset + (iPhi + 1) * hCalPhiSize;
      debug() << "HCal phi range  : " << lowPhi << " -  " << highPhi << endmsg;
      int lowId = floor((lowPhi - 0.5 * eCalModuleSize - eCalPhiOffset) / eCalModuleSize);
      int highId = floor((highPhi + 0.5 * eCalModuleSize - eCalPhiOffset) / eCalModuleSize);
      debug() << "ECal phi range  : " << lowId * eCalModuleSize + eCalPhiOffset << " -  "
              << highId * eCalModuleSize + eCalPhiOffset << endmsg;
      std::vector<uint> neighbours;
      for (int idPhiToAdd = lowId; idPhiToAdd <= highId; idPhiToAdd++)
      {
        neighbours.push_back(det::utils::cyclicNeighbour(idPhiToAdd, extremaECalLastLayerModule));
      }
      debug() << "HCal phi id  : " << iPhi << endmsg;
      debug() << "Found ECal Neighbours in phi : " << neighbours.size() << endmsg;
      for (auto id : neighbours)
      {
        debug() << "ECal Neighbours id : " << id << endmsg;
      }
      phiNeighbours.insert(std::pair<uint, std::vector<uint>>(iPhi, neighbours));
    }
    // add neighbours to both ecal cell and hcal cells
    dd4hep::DDSegmentation::CellID ecalCellId = 0;
    dd4hep::DDSegmentation::CellID hcalCellId = 0;
    (*decoderECalBarrel)["system"].set(ecalCellId, 5);
    (*decoderECalBarrel)[m_activeFieldNamesSegmented[0]].set(ecalCellId, eCalLastLayer);
    (*decoderHCalBarrel)["system"].set(hcalCellId, 8);
    (*decoderHCalBarrel)[m_activeFieldNamesSegmented[1]].set(hcalCellId, 0);
    // loop over segmented hcal cells
    for (auto iThetaHCal : thetaNeighbours)
    {
      (*decoderHCalBarrel)["theta"].set(hcalCellId, iThetaHCal.first);
      for (auto iPhiHCal : phiNeighbours)
      {
        (*decoderHCalBarrel)["phi"].set(hcalCellId, iPhiHCal.first);
        for (auto iTheta : iThetaHCal.second)
        {
          (*decoderECalBarrel)["theta"].set(ecalCellId, iTheta);
          for (auto iPhi : iPhiHCal.second)
          {
            (*decoderECalBarrel)["phi"].set(ecalCellId, iPhi);
            map.find(hcalCellId)->second.push_back(ecalCellId);
            map.find(ecalCellId)->second.push_back(hcalCellId);
            count++;
          }
        }
      }
    }
    */
  }
  // end of connectBarrels

  if (msgLevel() <= MSG::DEBUG)
  {
    std::vector<int> counter;
    counter.assign(40, 0);
    for (const auto &item : map)
    {
      counter[item.second.size()]++;
    }
    for (uint iCount = 0; iCount < counter.size(); iCount++)
    {
      if (counter[iCount] != 0)
      {
        debug() << counter[iCount] << " cells have " << iCount << " neighbours" << endmsg;
      }
    }
  }
  debug() << "cells with neighbours across Calo boundaries: " << count << endmsg;

  std::unique_ptr<TFile> file(TFile::Open(m_outputFileName.c_str(), "RECREATE"));
  file->cd();
  TTree tree("neighbours", "Tree with map of neighbours");
  uint64_t saveCellId;
  std::vector<uint64_t> saveNeighbours;
  tree.Branch("cellId", &saveCellId, "cellId/l");
  tree.Branch("neighbours", &saveNeighbours);
  for (const auto &item : map)
  {
    saveCellId = item.first;
    saveNeighbours = item.second;
    tree.Fill();
  }
  file->Write();
  file->Close();

  return StatusCode::SUCCESS;
}

StatusCode CreateFCCeeCaloNeighbours::finalize() { return Service::finalize(); }
