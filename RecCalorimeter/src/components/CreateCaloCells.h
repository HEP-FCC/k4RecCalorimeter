#ifndef RECCALORIMETER_CREATECALOCELLS_H
#define RECCALORIMETER_CREATECALOCELLS_H

// k4FWCore
#include "k4FWCore/DataHandle.h"

// Interfaces
#include "RecCaloCommon/ICalibrateCaloHitsTool.h"
#include "RecCaloCommon/ICaloReadCrosstalkMap.h"
#include "RecCaloCommon/ICalorimeterTool.h"
#include "RecCaloCommon/INoiseCaloCellsTool.h"

// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"

// edm4hep
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/Constants.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DD4hep/Volumes.h"
#include "TGeoManager.h"

class IGeoSvc;

/** @class CreateCaloCells
 *
 *  Algorithm for creating calorimeter cells from Geant4 hits.
 *  Tube geometry with PhiEta segmentation expected.
 *
 *  Flow of the program:
 *  1/ Merge Geant4 energy deposits with same cellID
 *  2/ Emulate cross-talk (if switched on)
 *  3/ Calibrate to electromagnetic scale (if calibration switched on)
 *  4/ Add random noise to each cell (if noise switched on)
 *  5/ Filter cells and remove those with energy below threshold (if noise +
 * filtering switched on)
 *
 *  Tools called:
 *    - CalibrateCaloHitsTool
 *    - NoiseCaloCellsTool
 *    - CaloReadCrosstalkMap
 *    - CalorimeterTool (for geometry)
 *
 *  @author Jana Faltova
 *  @author Anna Zaborowska
 *  @date   2016-09
 *
 */

class CreateCaloCells : public Gaudi::Algorithm {

public:
  CreateCaloCells(const std::string& name, ISvcLocator* svcLoc);

  virtual StatusCode initialize() override;

  virtual StatusCode execute(const EventContext&) const override;

private:
  /// Build m_caloTypes, giving calorimeter type per system ID.
  void findCaloTypes();

  /// Handle for the calorimeter cells crosstalk tool
  mutable ToolHandle<k4::recCalo::ICaloReadCrosstalkMap> m_crosstalksTool{"ReadCaloCrosstalkMap", this};
  /// Handle for tool to calibrate Geant4 energy to EM scale tool
  mutable ToolHandle<k4::recCalo::ICalibrateCaloHitsTool> m_calibTool{"CalibrateCaloHitsTool", this};
  /// Handle for the calorimeter cells noise tool
  mutable ToolHandle<k4::recCalo::INoiseCaloCellsTool> m_noiseTool{"NoiseCaloCellsFlatTool", this};
  /// Handle for the geometry tool
  ToolHandle<k4::recCalo::ICalorimeterTool> m_geoTool{"", this};

  /// Add crosstalk to cells?
  Gaudi::Property<bool> m_addCrosstalk{this, "addCrosstalk", false, "Add crosstalk effect?"};
  /// Calibrate to EM scale?
  Gaudi::Property<bool> m_doCellCalibration{this, "doCellCalibration", true, "Calibrate to EM scale?"};
  /// Add noise to cells?
  Gaudi::Property<bool> m_addCellNoise{this, "addCellNoise", true, "Add noise to cells?"};
  /// Save only cells with energy above threshold?
  Gaudi::Property<bool> m_filterCellNoise{this, "filterCellNoise", false,
                                          "Save only cells with energy above threshold?"};
  // Add position information to the cells? (based on Volumes, not cells, could be improved)
  Gaudi::Property<bool> m_addPosition{this, "addPosition", false, "Add position information to the cells?"};

  /// Handle for calo hits (input collection)
  mutable k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection> m_hits{"hits", Gaudi::DataHandle::Reader, this};
  /// Handle for calo cells (output collection)
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_cells{"cells", Gaudi::DataHandle::Writer, this};
  /// Name of the detector readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "ECalBarrelPhiEta", "Name of the detector readout"};
  /// Name of active volumes

  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;
  dd4hep::VolumeManager m_volman;
  /// Maps of cell IDs (corresponding to DD4hep IDs) on final energies to be used for clustering
  mutable std::unordered_map<uint64_t, double> m_cellsMap;
  /// Maps of cell IDs (corresponding to DD4hep IDs) on transfer of signals due to crosstalk
  mutable std::unordered_map<uint64_t, double> m_CrosstalkCellsMap;
  /// Maps of cell IDs with zero energy, for all cells in calo (needed if addCellNoise and filterCellNoise are both set)
  mutable std::unordered_map<uint64_t, double> m_emptyCellsMap;

  /// Indexed by system ID, giving the calorimeter type word.
  /// Non-calorimeter system IDs are set to 0.
  /// We record this for all system IDs since we build this during
  /// initialization, before we know which system ID we're handling.
  std::vector<int> m_caloTypes;

  /// Cell ID decoder.
  dd4hep::DDSegmentation::BitFieldCoder m_decoder;

  /// Field indices for detector ID and layer.
  unsigned m_systemIndex = -1;
  unsigned m_layerIndex = -1;
};

#endif /* RECCALORIMETER_CREATECALOCELLS_H */
