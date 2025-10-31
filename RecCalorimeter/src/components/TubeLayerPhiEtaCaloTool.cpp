#include "TubeLayerPhiEtaCaloTool.h"

// segm
#include "DD4hep/Detector.h"
#include "DD4hep/MultiSegmentation.h"
#include "detectorCommon/DetUtils_k4geo.h"

DECLARE_COMPONENT(TubeLayerPhiEtaCaloTool)


StatusCode TubeLayerPhiEtaCaloTool::collectCells(std::vector<uint64_t>& cells) const
{
  // Get the total number of active volumes in the geometry
  auto highestVol = gGeoManager->GetTopVolume();
  unsigned int numLayers;
  if (!m_activeVolumesNumber) {
    numLayers = det::utils::countPlacedVolumes(highestVol, m_activeVolumeName);
  } else {
    // used when MergeLayers tool is used. To be removed once MergeLayer gets replaced by RedoSegmentation.
    numLayers = m_activeVolumesNumber;
  }
  info() << "Number of active layers " << numLayers << endmsg;

  // get PhiEta segmentation
  const dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo* segmentation = nullptr;
  const dd4hep::DDSegmentation::MultiSegmentation* segmentationMulti = nullptr;
  segmentation = dynamic_cast<const dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo*>(
      readout().segmentation().segmentation());
  if (segmentation == nullptr) {
    segmentationMulti = dynamic_cast<const dd4hep::DDSegmentation::MultiSegmentation*>(
        readout().segmentation().segmentation());
    if (segmentationMulti == nullptr) {
      error() << "There is no phi-eta or multi- segmentation for the readout " << readoutName() << " defined."
              << endmsg;
      return StatusCode::FAILURE;
    } else {
      // check if multisegmentation contains only phi-eta sub-segmentations
      const dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo* subsegmentation = nullptr;
      for (const auto& subSegm : segmentationMulti->subSegmentations()) {
        subsegmentation = dynamic_cast<const dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo*>(subSegm.segmentation);
        if (subsegmentation == nullptr) {
          error() << "At least one of the sub-segmentations in MultiSegmentation named " << readoutName()
                  << " is not a phi-eta grid." << endmsg;
          return StatusCode::FAILURE;
        } else {
          info() << "subsegmentation for " << segmentationMulti->discriminatorName() << " from " << subSegm.key_min
                 << " to " << subSegm.key_max << endmsg;
          info() << "size in eta " << subsegmentation->gridSizeEta() << " , bins in phi " << subsegmentation->phiBins()
                 << endmsg;
          info() << "offset in eta " << subsegmentation->offsetEta() << " , offset in phi "
                 << subsegmentation->offsetPhi() << endmsg;
        }
      }
    }
  } else {
    info() << "FCCSWGridPhiEta_k4geo: size in eta " << segmentation->gridSizeEta() << " , bins in phi "
           << segmentation->phiBins() << endmsg;
    info() << "FCCSWGridPhiEta_k4geo: offset in eta " << segmentation->offsetEta() << " , offset in phi "
           << segmentation->offsetPhi() << endmsg;
  }
  // Take readout bitfield decoder from GeoSvc
  auto decoder = readout().idSpec().decoder();
  if (m_fieldNames.size() != m_fieldValues.size()) {
    error() << "Volume readout field descriptors (names and values) have different size." << endmsg;
    return StatusCode::FAILURE;
  }

  // Loop over all cells in the calorimeter and retrieve existing cellIDs
  // Loop over active layers
  for (unsigned int ilayer = 0; ilayer < numLayers; ilayer++) {
    // Get VolumeID
    dd4hep::DDSegmentation::VolumeID volumeID = 0;
    for (unsigned int it = 0; it < m_fieldNames.size(); it++) {
      decoder->set(volumeID, m_fieldNames[it], m_fieldValues[it]);
    }
    decoder->set(volumeID, m_activeFieldName, ilayer);
    decoder->set(volumeID, "eta", 0);
    decoder->set(volumeID, "phi", 0);

    if (segmentationMulti != nullptr) {
      segmentation = dynamic_cast<const dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo*>(
          &segmentationMulti->subsegmentation(volumeID));
    }

    // Get number of segmentation cells within the active volume
    auto numCells = det::utils::numberOfCells(volumeID, *segmentation);
    debug() << "Segmentation cells  (Nphi, Neta, minEta): " << numCells << endmsg;
    // Loop over segmenation cells
    for (unsigned int iphi = 0; iphi < numCells[0]; iphi++) {
      for (unsigned int ieta = 0; ieta < numCells[1]; ieta++) {
        decoder->set(volumeID, "phi", iphi);
        decoder->set(volumeID, "eta", ieta + numCells[2]); // start from the minimum existing eta cell in this layer
        dd4hep::DDSegmentation::CellID cellId = volumeID;
        cells.push_back(cellId);
      }
    }
  }
  return StatusCode::SUCCESS;
}
