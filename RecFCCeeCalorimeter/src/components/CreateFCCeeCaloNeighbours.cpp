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
  int ecalBarrelId = -1; // index of ecal barrel in segmentation list
  dd4hep::DDSegmentation::BitFieldCoder *decoderECalBarrel = nullptr;
  dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo *hcalBarrelSegmentation = nullptr;
  //  will be used for volume connecting
  std::pair<int, int> extremaHCalFirstLayerPhi;
  std::pair<int, int> extremaHCalFirstLayerTheta;
  int hcalBarrelId = -1; // index of hcal barrel in segmentation list
  dd4hep::DDSegmentation::BitFieldCoder *decoderHCalBarrel = nullptr;
  dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo *ecalBarrelModuleThetaSegmentation = nullptr;
  dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo *ecalBarrelPhiThetaSegmentation = nullptr;

  for (uint iSys = 0; iSys < m_readoutNamesSegmented.size(); iSys++)
  {
    // Check if readout exists
    info() << "Readout: " << m_readoutNamesSegmented[iSys] << endmsg;
    if (m_geoSvc->getDetector()->readouts().find(m_readoutNamesSegmented[iSys]) == m_geoSvc->getDetector()->readouts().end())
    {
      error() << "Readout <<" << m_readoutNamesSegmented[iSys] << ">> does not exist." << endmsg;
      return StatusCode::FAILURE;
    }

    // get theta-based segmentation
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
    if (m_connectBarrels)
    {
      if (m_fieldNamesSegmented[iSys] == "system" && m_fieldValuesSegmented[iSys] == m_ecalBarrelSysId)
      {
        decoderECalBarrel = decoder;
        ecalBarrelId = iSys;
        if (segmentationType == "FCCSWGridModuleThetaMerged_k4geo")
        {
          ecalBarrelModuleThetaSegmentation = moduleThetaSegmentation;
        }
        else if (segmentationType == "FCCSWGridPhiTheta_k4geo")
        {
          ecalBarrelPhiThetaSegmentation = phiThetaSegmentation;
        }
        else
        {
          error() << "ECAL barrel segmentation type not handled for connectBarrels option." << endmsg;
          return StatusCode::FAILURE;
        }
      }

      if (m_fieldNamesSegmented[iSys] == "system" && m_fieldValuesSegmented[iSys] == m_hcalBarrelSysId)
      {
        decoderHCalBarrel = decoder;
        hcalBarrelId = iSys;
        if (segmentationType == "FCCSWGridPhiTheta_k4geo")
        {
          hcalBarrelSegmentation = phiThetaSegmentation;
        }
        else
        {
          error() << "HCAL barrel segmentation type not handled for connectBarrels option." << endmsg;
          return StatusCode::FAILURE;
        }
      }
    }

    // Loop over all cells in the calorimeter and retrieve existing cellIDs and find neihbours
    if (segmentationType == "FCCSWGridPhiTheta_k4geo")
    {
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
        // for vlume such as HCAL for which the numberOfCells function does not work
        // (because the G4 volume does not span the whole z of the detector)
        // use theta range passed via activeVolumesTheta option
        // if (m_fieldValuesSegmented[iSys] == m_hcalBarrelSysId)
        if (m_activeVolumesTheta[iSys].size() > 0)
        {
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
        }
        else
        {
          numCells = det::utils::numberOfCells(volumeId, *phiThetaSegmentation);
        }
        // extrema 1: 0, ID of last cell in phi
        extrema[1] = std::make_pair(0, numCells[0] - 1);
        // extrema[2]: min theta ID, n theta cells
        extrema[2] = std::make_pair(numCells[2], numCells[1] + numCells[2] - 1);

        // for layer 0 of HCal barrel,  will be used for volume connecting
        if (
            m_fieldNamesSegmented[iSys] == "system" &&
            m_fieldValuesSegmented[iSys] == m_hcalBarrelSysId &&
            ilayer == 0)
        {
          extremaHCalFirstLayerPhi = std::make_pair(0, numCells[0] - 1);
          extremaHCalFirstLayerTheta = std::make_pair(numCells[2], numCells[1] + numCells[2] - 1);
        }
        // for layer N-1 of ECal barrel,  will be used for volume connecting
        if (m_fieldNamesSegmented[iSys] == "system" &&
            m_fieldValuesSegmented[iSys] == m_ecalBarrelSysId &&
            ilayer == (m_activeVolumesNumbersSegmented[iSys] - 1))
        {
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
        for (unsigned int iphi = 0; iphi < numCells[0]; iphi++)
        {
          for (unsigned int itheta = 0; itheta < numCells[1]; itheta++)
          {
            dd4hep::DDSegmentation::CellID cellId = volumeId;
            decoder->set(cellId, "phi", iphi);
            decoder->set(cellId, "theta", itheta + numCells[2]); // start from the minimum existing theta cell in this layer
            uint64_t id = cellId;
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
        if (ilayer == (m_activeVolumesNumbersSegmented[iSys] - 1) && m_fieldNamesSegmented[iSys] == "system" &&
            m_fieldValuesSegmented[iSys] == m_ecalBarrelSysId)
        {
          eCalLastLayer = m_activeVolumesNumbersSegmented[iSys] - 1;
          extremaECalLastLayerModule = std::make_pair(0, (numCells[0] - 1) * moduleThetaSegmentation->mergedModules(eCalLastLayer));
          extremaECalLastLayerTheta = std::make_pair(numCells[2], numCells[2] + (numCells[1] - 1) * moduleThetaSegmentation->mergedThetaCells(eCalLastLayer));
        }
        debug() << "Layer: " << ilayer << endmsg;
        debug() << "Extrema[0]: " << extrema[0].first << " , " << extrema[0].second << endmsg;
        debug() << "Extrema[1]: " << extrema[1].first << " , " << extrema[1].second << endmsg;
        debug() << "Extrema[2]: " << extrema[2].first << " , " << extrema[2].second << endmsg;
        debug() << "Number of segmentation cells in (module,theta): " << numCells << endmsg;
        // Loop over segmentation cells to find neighbours in ECAL
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
    // first check if ECAL barrel and HCal barrel are configured
    if (decoderECalBarrel == nullptr || decoderHCalBarrel == nullptr)
    {
      error() << "ECAL barrel and/or HCal barrel are not configured correctly, cannot link barrels!" << endmsg;
      return StatusCode::FAILURE;
    }

    // print how many cells in each dimensions will be matched
    info() << "ECAL layer " << eCalLastLayer << " is a neighbour of HCAL layer 0." << endmsg;
    if (ecalBarrelModuleThetaSegmentation)
    {
      info() << "ECAL modules " << extremaECalLastLayerModule.first << " - " << extremaECalLastLayerModule.second;
    }
    else
    {
      info() << "ECAL phi cells " << extremaECalLastLayerPhi.first << " - " << extremaECalLastLayerPhi.second;
    }
    info() << " will be matched to HCAL phi cells " << extremaHCalFirstLayerPhi.first
           << " - " << extremaHCalFirstLayerPhi.second << endmsg;
    info() << "ECAL theta cells " << extremaECalLastLayerTheta.first << " - " << extremaECalLastLayerTheta.second
           << " will be matched to HCAL theta cells " << extremaHCalFirstLayerTheta.first
           << " - " << extremaHCalFirstLayerTheta.second << endmsg;

    dd4hep::DDSegmentation::CellID volumeIdECal = 0;
    dd4hep::DDSegmentation::CellID volumeIdHCal = 0;
    (*decoderECalBarrel)["system"].set(volumeIdECal, m_ecalBarrelSysId);
    (*decoderECalBarrel)[m_activeFieldNamesSegmented[ecalBarrelId]].set(volumeIdECal, eCalLastLayer);
    (*decoderHCalBarrel)["system"].set(volumeIdHCal, m_hcalBarrelSysId);
    (*decoderHCalBarrel)[m_activeFieldNamesSegmented[hcalBarrelId]].set(volumeIdHCal, 0);

    // retrieve needed parameters
    double hCalThetaSize = hcalBarrelSegmentation->gridSizeTheta();
    double hCalThetaOffset = hcalBarrelSegmentation->offsetTheta();
    double hCalPhiOffset = hcalBarrelSegmentation->offsetPhi();
    double hCalPhiSize = hcalBarrelSegmentation->gridSizePhi();
    int hCalPhiBins = hcalBarrelSegmentation->phiBins();
    double eCalThetaOffset, eCalThetaSize;
    double eCalPhiOffset, eCalPhiSize;
    int eCalModules(-1);
    int eCalPhiBins(-1);
    if (ecalBarrelModuleThetaSegmentation)
    {
      eCalThetaOffset = ecalBarrelModuleThetaSegmentation->offsetTheta();
      eCalThetaSize = ecalBarrelModuleThetaSegmentation->gridSizeTheta();
      eCalModules = ecalBarrelModuleThetaSegmentation->nModules();
      // for ECAL barrel with module readout, need to find out phi position of module 0 of ECal barrel last layer
      // get it from volume manager
      // for a normal phi grid it would suffice to use the segmentation phi(cellId) method
      dd4hep::VolumeManager volman = m_geoSvc->getDetector()->volumeManager();
      auto detelement = volman.lookupDetElement(volumeIdECal);
      const auto &transformMatrix = detelement.nominal().worldTransformation();
      double outGlobal[3];
      double inLocal[] = {0, 0, 0};
      transformMatrix.LocalToMaster(inLocal, outGlobal);
      eCalPhiOffset = std::atan2(outGlobal[1], outGlobal[0]);
      info() << "ECAL modules : " << eCalModules << " , phi of module 0 in last ECAL layer : " << eCalPhiOffset << endmsg;
      eCalPhiSize = TMath::TwoPi() / eCalModules;
    }
    else
    {
      eCalThetaOffset = ecalBarrelPhiThetaSegmentation->offsetTheta();
      eCalThetaSize = ecalBarrelPhiThetaSegmentation->gridSizeTheta();
      eCalPhiOffset = ecalBarrelPhiThetaSegmentation->offsetPhi();
      eCalPhiSize = ecalBarrelPhiThetaSegmentation->gridSizePhi();
      eCalPhiBins = ecalBarrelPhiThetaSegmentation->phiBins();
    }

    dd4hep::DDSegmentation::CellID cellIdECal = volumeIdECal;
    dd4hep::DDSegmentation::CellID cellIdHCal = volumeIdHCal;

    // Loop over ECAL segmentation cells to find neighbours in HCAL
    if (ecalBarrelModuleThetaSegmentation)
    {
      for (int imodule = extremaECalLastLayerModule.first; imodule <= extremaECalLastLayerModule.second; imodule += ecalBarrelModuleThetaSegmentation->mergedModules(eCalLastLayer))
      {
        (*decoderECalBarrel)["module"].set(cellIdECal, imodule);
        double dphi = ecalBarrelModuleThetaSegmentation->phi(cellIdECal);
        double phiVol = eCalPhiOffset + imodule * eCalPhiSize;
        double phi = phiVol + dphi;
        double phiMin = phi - 0.5 * eCalPhiSize * ecalBarrelModuleThetaSegmentation->mergedModules(eCalLastLayer);
        double phiMax = phi + 0.5 * eCalPhiSize * ecalBarrelModuleThetaSegmentation->mergedModules(eCalLastLayer);
        int minPhiHCal = (phiMin - hCalPhiOffset) / hCalPhiSize;
        int maxPhiHCal = (phiMax - hCalPhiOffset) / hCalPhiSize;
        for (int itheta = extremaECalLastLayerTheta.first; itheta <= extremaECalLastLayerTheta.second; itheta += ecalBarrelModuleThetaSegmentation->mergedThetaCells(eCalLastLayer))
        {
          (*decoderECalBarrel)["theta"].set(cellIdECal, itheta);
          double theta = ecalBarrelModuleThetaSegmentation->theta(cellIdECal);
          double thetaMin = theta - eCalThetaSize * ecalBarrelModuleThetaSegmentation->mergedThetaCells(eCalLastLayer) / 2.;
          double thetaMax = theta + eCalThetaSize * ecalBarrelModuleThetaSegmentation->mergedThetaCells(eCalLastLayer) / 2.;
          int minThetaHCal = (thetaMin - hCalThetaOffset) / hCalThetaSize;
          int maxThetaHCal = (thetaMax - hCalThetaOffset) / hCalThetaSize;
          debug() << "ECAL cell: cellId = " << cellIdECal << " module = " << imodule << " theta bin = " << itheta << endmsg;
          debug() << "ECAL cell: phi = " << phi << " theta = " << theta << endmsg;
          debug() << "ECAL cell: thetaMin = " << thetaMin << " thetaMax = " << thetaMax << endmsg;
          debug() << "ECAL cell: phiMin = " << phiMin << " phiMax = " << phiMax << endmsg;
          debug() << "ECAL cell: theta bin of neighbours in HCAL layer 0 = " << minThetaHCal << " - " << maxThetaHCal << endmsg;
          debug() << "ECAL cell: phi bin of neighbours in HCAL layer 0 = " << minPhiHCal << " - " << maxPhiHCal << endmsg;
          bool hasNeighbours = false;
          for (int iphiHCal = minPhiHCal; iphiHCal <= maxPhiHCal; iphiHCal++)
          {
            int phiBinHCal = iphiHCal;
            // bring phi bin in 0..nbins-1 range
            if (phiBinHCal < 0)
              phiBinHCal += hCalPhiBins;
            if (phiBinHCal >= hCalPhiBins)
              phiBinHCal -= hCalPhiBins;
            if (phiBinHCal < extremaHCalFirstLayerPhi.first || phiBinHCal > extremaHCalFirstLayerPhi.second)
            {
              warning() << "phi bin " << phiBinHCal << " out of range, skipping" << endmsg;
              continue;
            }
            (*decoderHCalBarrel)["phi"].set(cellIdHCal, phiBinHCal);
            for (int ithetaHCal = minThetaHCal; ithetaHCal <= maxThetaHCal; ithetaHCal++)
            {
              if (ithetaHCal < extremaHCalFirstLayerTheta.first || ithetaHCal > extremaHCalFirstLayerTheta.second)
              {
                warning() << "theta bin " << ithetaHCal << " out of range, skipping" << endmsg;
                continue;
              }
              (*decoderHCalBarrel)["theta"].set(cellIdHCal, ithetaHCal);
              debug() << "ECAL cell: neighbour in HCAL has cellID = " << cellIdHCal << endmsg;
              map.find((uint64_t)cellIdECal)->second.push_back((uint64_t)cellIdHCal);
              hasNeighbours = true;
            }
          }
          if (hasNeighbours)
            count++;
        }
      }
    }
    else
    {
      for (int iphi = extremaECalLastLayerPhi.first; iphi <= extremaECalLastLayerPhi.second; iphi++)
      {
        (*decoderECalBarrel)["phi"].set(cellIdECal, iphi);
        double phi = ecalBarrelPhiThetaSegmentation->phi(cellIdECal);
        double phiMin = phi - 0.5 * eCalPhiSize;
        double phiMax = phi + 0.5 * eCalPhiSize;
        int minPhiHCal = (phiMin - hCalPhiOffset) / hCalPhiSize;
        int maxPhiHCal = (phiMax - hCalPhiOffset) / hCalPhiSize;
        for (int itheta = extremaECalLastLayerTheta.first; itheta <= extremaECalLastLayerTheta.second; itheta ++)
        {
          (*decoderECalBarrel)["theta"].set(cellIdECal, itheta);
          double theta = ecalBarrelPhiThetaSegmentation->theta(cellIdECal);
          double thetaMin = theta - 0.5 * eCalThetaSize;
          double thetaMax = theta + 0.5 * eCalThetaSize;
          int minThetaHCal = (thetaMin - hCalThetaOffset) / hCalThetaSize;
          int maxThetaHCal = (thetaMax - hCalThetaOffset) / hCalThetaSize;
          debug() << "ECAL cell: cellId = " << cellIdECal << " phi bin = " << iphi << " theta bin = " << itheta << endmsg;
          debug() << "ECAL cell: phi = " << phi << " theta = " << theta << endmsg;
          debug() << "ECAL cell: thetaMin = " << thetaMin << " thetaMax = " << thetaMax << endmsg;
          debug() << "ECAL cell: phiMin = " << phiMin << " phiMax = " << phiMax << endmsg;
          debug() << "ECAL cell: theta bin of neighbours in HCAL layer 0 = " << minThetaHCal << " - " << maxThetaHCal << endmsg;
          debug() << "ECAL cell: phi bin of neighbours in HCAL layer 0 = " << minPhiHCal << " - " << maxPhiHCal << endmsg;
          bool hasNeighbours = false;
          for (int iphiHCal = minPhiHCal; iphiHCal <= maxPhiHCal; iphiHCal++)
          {
            int phiBinHCal = iphiHCal;
            // bring phi bin in 0..nbins-1 range
            if (phiBinHCal < 0)
              phiBinHCal += hCalPhiBins;
            if (phiBinHCal >= hCalPhiBins)
              phiBinHCal -= hCalPhiBins;
            if (phiBinHCal < extremaHCalFirstLayerPhi.first || phiBinHCal > extremaHCalFirstLayerPhi.second)
            {
              warning() << "phi bin " << phiBinHCal << " out of range, skipping" << endmsg;
              continue;
            }
            (*decoderHCalBarrel)["phi"].set(cellIdHCal, phiBinHCal);
            for (int ithetaHCal = minThetaHCal; ithetaHCal <= maxThetaHCal; ithetaHCal++)
            {
              if (ithetaHCal < extremaHCalFirstLayerTheta.first || ithetaHCal > extremaHCalFirstLayerTheta.second)
              {
                warning() << "theta bin " << ithetaHCal << " out of range, skipping" << endmsg;
                continue;
              }
              (*decoderHCalBarrel)["theta"].set(cellIdHCal, ithetaHCal);
              debug() << "ECAL cell: neighbour in HCAL has cellID = " << cellIdHCal << endmsg;
              map.find((uint64_t)cellIdECal)->second.push_back((uint64_t)cellIdHCal);
              hasNeighbours = true;
            }
          }
          if (hasNeighbours)
            count++;
        }
      }
    }

    // Loop over HCAL segmentation cells to find neighbours in ECAL
    // Loop on phi
    for (int iphi = extremaHCalFirstLayerPhi.first; iphi <= extremaHCalFirstLayerPhi.second; iphi++)
    {
      // determine phi extent of HCAL cell
      (*decoderHCalBarrel)["phi"].set(cellIdHCal, iphi);
      double phi = hcalBarrelSegmentation->phi(cellIdHCal);
      double phiMin = phi - 0.5 * hCalPhiSize;
      double phiMax = phi + 0.5 * hCalPhiSize;

      // find ECAL barrel modules corresponding to this phi range
      int minModuleECal(-1), maxModuleECal(-1), minPhiECal(-1), maxPhiECal(-1);
      if (ecalBarrelModuleThetaSegmentation)
      {
        minModuleECal = (int)((phiMin - eCalPhiOffset) / eCalPhiSize);
        if (minModuleECal < 0)
          minModuleECal += eCalModules; // need module to be >=0 so that subtracting module%mergedModules gives the correct result
        minModuleECal -= minModuleECal % ecalBarrelModuleThetaSegmentation->mergedModules(eCalLastLayer);
        maxModuleECal = (int)((phiMax - eCalPhiOffset) / eCalPhiSize);
        if (maxModuleECal < 0)
          maxModuleECal += eCalModules;
        maxModuleECal -= maxModuleECal % ecalBarrelModuleThetaSegmentation->mergedModules(eCalLastLayer);
        // due to ciclic behaviour of modules, one could have e.g min = 1534 and max = 2 (for N=1536)
        if (maxModuleECal < minModuleECal)
          maxModuleECal += eCalModules;
      }
      else
      {
        minPhiECal = (int)((phiMin - eCalPhiOffset) / eCalPhiSize);
        maxPhiECal = (int)((phiMax - eCalPhiOffset) / eCalPhiSize);
      }

      // Loop on theta
      for (int itheta = extremaHCalFirstLayerTheta.first; itheta <= extremaHCalFirstLayerTheta.second; itheta++)
      {
        // determine theta extent of HCAL cell
        (*decoderHCalBarrel)["theta"].set(cellIdHCal, itheta);
        double theta = hcalBarrelSegmentation->theta(cellIdHCal);
        double thetaMin = theta - 0.5 * hCalThetaSize;
        double thetaMax = theta + 0.5 * hCalThetaSize;
        // find ECAL barrel theta bins corresponding to this theta range
        int minThetaECal = (thetaMin - eCalThetaOffset) / eCalThetaSize;
        int maxThetaECal = (thetaMax - eCalThetaOffset) / eCalThetaSize;
        if (ecalBarrelModuleThetaSegmentation)
        {
          minThetaECal -= minThetaECal % ecalBarrelModuleThetaSegmentation->mergedThetaCells(eCalLastLayer);
          maxThetaECal -= maxThetaECal % ecalBarrelModuleThetaSegmentation->mergedThetaCells(eCalLastLayer);
        }
        // print out some debug info
        debug() << "HCAL cell: cellId = " << cellIdHCal << " phi bin = " << iphi << " theta bin = " << itheta << endmsg;
        debug() << "HCAL cell: phi = " << phi << " theta = " << theta << endmsg;
        debug() << "HCAL cell: thetaMin = " << thetaMin << " thetaMax = " << thetaMax << endmsg;
        debug() << "HCAL cell: phiMin = " << phiMin << " phiMax = " << phiMax << endmsg;
        debug() << "HCAL cell: theta bin of neighbours in ECAL layer max = " << minThetaECal << " - " << maxThetaECal << endmsg;
        if (ecalBarrelModuleThetaSegmentation)
        {
          debug() << "HCAL cell: module of neighbours in ECAL layer max = " << minModuleECal << " - " << maxModuleECal << endmsg;
        }
        else
        {
          debug() << "HCAL cell: phi bin of neighbours in ECAL layer max = " << minPhiECal << " - " << maxPhiECal << endmsg;
        }

        // add the neighbours if within the acceptance of the layer
        bool hasNeighbours = false;
        int thetaStep = (ecalBarrelModuleThetaSegmentation) ? ecalBarrelModuleThetaSegmentation->mergedThetaCells(eCalLastLayer) : 1;
        for (int ithetaECal = minThetaECal; ithetaECal <= maxThetaECal; ithetaECal += thetaStep)
        {
          if (ithetaECal < extremaECalLastLayerTheta.first || ithetaECal > extremaECalLastLayerTheta.second)
          {
            debug() << "theta bin " << ithetaECal << " out of range, skipping" << endmsg;
            continue;
          }
          (*decoderECalBarrel)["theta"].set(cellIdECal, ithetaECal);
          if (ecalBarrelModuleThetaSegmentation)
          {
            for (int iModule = minModuleECal; iModule <= maxModuleECal; iModule += ecalBarrelModuleThetaSegmentation->mergedModules(eCalLastLayer))
            {
              int module = iModule;
              if (module < 0)
                module += eCalModules;
              if (module >= eCalModules)
                module -= eCalModules;
              if (module < extremaECalLastLayerModule.first || module > extremaECalLastLayerModule.second)
              {
                warning() << "module " << module << " out of range, skipping" << endmsg;
                continue;
              }
              (*decoderECalBarrel)["module"].set(cellIdECal, module);
              debug() << "HCAL cell: neighbour in ECAL has cellID = " << cellIdECal << endmsg;
              map.find((uint64_t)cellIdHCal)->second.push_back((uint64_t)cellIdECal);
              hasNeighbours = true;
            }
          }
          else
          {
            for (int iphiECal = minPhiECal; iphiECal <= maxPhiECal; iphiECal++)
            {
              int phiBinECal = iphiECal;
              // bring phi bin in 0..nbins-1 range
              if (phiBinECal < 0)
                phiBinECal += eCalPhiBins;
              if (phiBinECal >= eCalPhiBins)
                phiBinECal -= eCalPhiBins;
              if (phiBinECal < extremaECalLastLayerPhi.first || phiBinECal > extremaECalLastLayerPhi.second)
              {
                warning() << "phi bin " << phiBinECal << " out of range, skipping" << endmsg;
                continue;
              }
              (*decoderECalBarrel)["phi"].set(cellIdECal, phiBinECal);
              debug() << "HCAL cell: neighbour in ECAL has cellID = " << cellIdECal << endmsg;
              map.find((uint64_t)cellIdHCal)->second.push_back((uint64_t)cellIdECal);
              hasNeighbours = true;
            }
          }
          if (hasNeighbours)
            count++;
        }
      }
    }
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
