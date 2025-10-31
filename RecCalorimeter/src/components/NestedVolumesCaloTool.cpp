#include "NestedVolumesCaloTool.h"

// segm
#include "DD4hep/Detector.h"
#include "detectorCommon/DetUtils_k4geo.h"

DECLARE_COMPONENT(NestedVolumesCaloTool)


StatusCode NestedVolumesCaloTool::collectCells(std::vector<uint64_t>& cells) const
{
  // Take readout bitfield decoder from GeoSvc
  auto decoder = readout().idSpec().decoder();
  if (m_fieldNames.size() != m_fieldValues.size()) {
    error() << "Volume readout field descriptors (names and values) have different size." << endmsg;
    return StatusCode::FAILURE;
  }
  // Get VolumeID
  dd4hep::DDSegmentation::CellID cID = 0;
  for (unsigned int it = 0; it < m_fieldNames.size(); it++) {
    decoder->set(cID, m_fieldNames[it], m_fieldValues[it]);
  }
  // Get the total number of given hierarchy of active volumes
  auto highestVol = gGeoManager->GetTopVolume();
  std::vector<unsigned int> numVolumes;
  numVolumes.reserve(m_activeVolumeName.size());
  for (const auto& volName : m_activeVolumeName) {
    numVolumes.push_back(det::utils::countPlacedVolumes(highestVol, volName));
    info() << "Number of active volumes named " << volName << " is " << numVolumes.back() << endmsg;
    if (numVolumes.back() == 0) {
      error() << "Volume name " << volName << " not found! Check naming in detector description." << endmsg;
      return StatusCode::FAILURE;
    }
  }
  // First sort to figure out which volume is inside which one
  std::vector<std::pair<std::string, int>> numVolumesMap;
  for (unsigned int it = 0; it < m_activeVolumeName.size(); it++) {
    // names of volumes (m_activeVolumeName) not needed anymore, only corresponding bitfield names are used
    // (m_activeFieldName)
    numVolumesMap.push_back(std::pair<std::string, int>(m_activeFieldName[it], numVolumes[it]));
  }
  std::sort(
      numVolumesMap.begin(), numVolumesMap.end(),
      [](std::pair<std::string, int> vol1, std::pair<std::string, int> vol2) { return vol1.second < vol2.second; });
  // now recompute how many volumes exist within the larger volume
  for (unsigned int it = numVolumesMap.size() - 1; it > 0; it--) {
    if (numVolumesMap[it].second % numVolumesMap[it - 1].second != 0) {
      error() << "Given volumes are not nested in each other!" << endmsg;
      return StatusCode::FAILURE;
    }
    numVolumesMap[it].second /= numVolumesMap[it - 1].second;
  }
  // Debug calculation of total number of cells
  if (msgLevel() <= MSG::DEBUG) {
    unsigned int checkTotal = 1;
    for (const auto& vol : numVolumesMap) {
      debug() << "Number of active volumes named " << vol.first << " is " << vol.second << endmsg;
      checkTotal *= vol.second;
    }
    debug() << "Total number of cells ( " << numVolumesMap.back().first << " ) is " << checkTotal << endmsg;
  }
  // Loop over all volumes in calorimeter to retrieve active cells
  std::vector<int> currentVol;
  unsigned int numVolTypes = numVolumes.size();
  currentVol.assign(numVolTypes, 0);
  unsigned int index = 0;
  do {
    index = 0;
    for (unsigned int it = 0; it < numVolTypes; it++) {
      decoder->set(cID, numVolumesMap[it].first, currentVol[it]);
    }
    cells.push_back(cID);
    if (msgLevel() <= MSG::VERBOSE) {
      verbose() << "Adding volume: " << decoder->valueString(cID) << endmsg;
    }
    // get next volume iterators
    currentVol[index]++;
    while (currentVol[index] == numVolumesMap[index].second) {
      currentVol[index] = 0;
      currentVol[++index]++;
    }
  } while (index != numVolTypes);
  return StatusCode::SUCCESS;
}
