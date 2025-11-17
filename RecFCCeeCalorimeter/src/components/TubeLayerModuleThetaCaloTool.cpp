#include "TubeLayerModuleThetaCaloTool.h"

// dd4hep
#include "DD4hep/Detector.h"
#include "DD4hep/MultiSegmentation.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"


DECLARE_COMPONENT(TubeLayerModuleThetaCaloTool)


StatusCode TubeLayerModuleThetaCaloTool::collectCells(std::vector<uint64_t>& cells) const
{
  // Get the total number of active volumes in the geometry
  unsigned int numLayers = m_activeVolumesNumber;
  info() << "Number of active layers " << numLayers << endmsg;

  // get segmentation
  const dd4hep::DDSegmentation::Segmentation* aSegmentation =
      readout().segmentation().segmentation();
  std::string segmentationType = aSegmentation->type();
  info() << "Segmentation type : " << segmentationType << endmsg;
  const dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo* moduleThetaSegmentation = nullptr;
  if (segmentationType == "FCCSWGridModuleThetaMerged_k4geo") {
    moduleThetaSegmentation = dynamic_cast<const dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo*>(aSegmentation);
    info() << "Segmentation: bins in Module " << moduleThetaSegmentation->nModules() << endmsg;
  } else {
    error() << "Unable to cast segmentation pointer!!!! Tool only applicable to FCCSWGridModuleThetaMerged_k4geo "
               "segmentation."
            << endmsg;
    return StatusCode::FAILURE;
  }

  // get decoder and extrema
  auto decoder = readout().idSpec().decoder();
  std::vector<std::pair<int, int>> extrema;
  extrema.push_back(std::make_pair(0, numLayers - 1));
  extrema.push_back(std::make_pair(0, 0));
  extrema.push_back(std::make_pair(0, 0));

  for (unsigned int ilayer = 0; ilayer < numLayers; ilayer++) {
    dd4hep::DDSegmentation::CellID volumeId = 0;
    for (unsigned int it = 0; it < m_fieldNames.size(); it++) {
      decoder->set(volumeId, m_fieldNames[it], m_fieldValues[it]);
    }
    (*decoder)["layer"].set(volumeId, ilayer);
    (*decoder)["theta"].set(volumeId, 0);
    (*decoder)["module"].set(volumeId, 0);
    // Get number of segmentation cells within the active volume, for given layer
    auto numCells = det::utils::numberOfCells(volumeId, *moduleThetaSegmentation);
    // extrema[1]: Range of module ID (0, max module ID)
    extrema[1] = std::make_pair(0, (numCells[0] - 1) * moduleThetaSegmentation->mergedModules(ilayer));
    // extrema[2]: Range of theta ID (min theta ID, max theta ID
    extrema[2] = std::make_pair(numCells[2],
                                numCells[2] + (numCells[1] - 1) * moduleThetaSegmentation->mergedThetaCells(ilayer));
    debug() << "Layer: " << ilayer << endmsg;
    debug() << "Number of segmentation cells in (module, theta): " << numCells << endmsg;
    // Loop over segmentation cells
    for (unsigned int imodule = 0; imodule < numCells[0]; imodule++) {
      for (unsigned int itheta = 0; itheta < numCells[1]; itheta++) {
        dd4hep::DDSegmentation::CellID cellId = volumeId;
        decoder->set(cellId, "module", imodule * moduleThetaSegmentation->mergedModules(ilayer));
        decoder->set(cellId, "theta",
                     numCells[2] + itheta * moduleThetaSegmentation->mergedThetaCells(
                                                ilayer)); // start from the minimum existing theta cell in this layer
        cells.push_back(cellId);
      }
    }
  }

  return StatusCode::SUCCESS;
}
