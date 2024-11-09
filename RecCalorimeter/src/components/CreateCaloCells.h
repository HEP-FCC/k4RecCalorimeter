#ifndef RECCALORIMETER_CREATECALOCELLS_H
#define RECCALORIMETER_CREATECALOCELLS_H

// k4FWCore
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/MetaDataHandle.h"
#include "k4Interface/ICalibrateCaloHitsTool.h"
#include "k4Interface/ICalorimeterTool.h"
#include "k4Interface/INoiseCaloCellsTool.h"
#include "k4Interface/ICaloReadCrosstalkMap.h"

// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"

// edm4hep
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/Constants.h"

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

  StatusCode initialize();

  StatusCode execute(const EventContext&) const;

  StatusCode finalize();

private:

  /// Handle for the calorimeter cells crosstalk tool
  mutable ToolHandle<ICaloReadCrosstalkMap> m_crosstalksTool{"ReadCaloCrosstalkMap", this};
  /// Handle for tool to calibrate Geant4 energy to EM scale tool
  mutable ToolHandle<ICalibrateCaloHitsTool> m_calibTool{"CalibrateCaloHitsTool", this};
  /// Handle for the calorimeter cells noise tool
  mutable ToolHandle<INoiseCaloCellsTool> m_noiseTool{"NoiseCaloCellsFlatTool", this};
  /// Handle for the geometry tool
  ToolHandle<ICalorimeterTool> m_geoTool{"TubeLayerPhiEtaCaloTool", this};

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
  mutable DataHandle<edm4hep::SimCalorimeterHitCollection> m_hits{"hits", Gaudi::DataHandle::Reader, this};
  /// Handle for the cellID encoding string of the input collection
  MetaDataHandle<std::string> m_hitsCellIDEncoding{m_hits, edm4hep::labels::CellIDEncoding, Gaudi::DataHandle::Reader};
  /// Handle for calo cells (output collection)
  mutable DataHandle<edm4hep::CalorimeterHitCollection> m_cells{"cells", Gaudi::DataHandle::Writer, this};
  MetaDataHandle<std::string> m_cellsCellIDEncoding{m_cells, edm4hep::labels::CellIDEncoding, Gaudi::DataHandle::Writer};
  /// Name of the detector readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "ECalBarrelPhiEta", "Name of the detector readout"};
  /// Name of active volumes
  Gaudi::Property<std::string> m_activeVolumeName{this, "activeVolumeName", "_sensitive", "Name of the active volumes"};
  /// Name of active layers for sampling calorimeter
  Gaudi::Property<std::string> m_activeFieldName{this, "activeFieldName", "active_layer",
                                                 "Name of active layers for sampling calorimeter"};
  /// Name of the bit-fields (in the readout) describing the volume
  Gaudi::Property<std::vector<std::string>> m_fieldNames{
      this, "fieldNames", {}, "Name of the bit-fields (in the readout) describing the volume"};
  /// Values of the fields that identify the volume to change segmentation (e.g.
  /// ID of the ECal)
  Gaudi::Property<std::vector<int>> m_fieldValues{this,
                                                  "fieldValues",
                                                  {},
                                                  "Value of the field that identifies the volume "
                                                  "to to change segmentation (e.g. ID of the "
                                                  "ECal)"};

  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;
  dd4hep::VolumeManager m_volman;
  /// Maps of cell IDs (corresponding to DD4hep IDs) on final energies to be used for clustering
  mutable std::unordered_map<uint64_t, double> m_cellsMap;
  /// Maps of cell IDs (corresponding to DD4hep IDs) on transfer of signals due to crosstalk
  mutable std::unordered_map<uint64_t, double> m_CrosstalkCellsMap;
  /// Maps of cell IDs with zero energy, for all cells in calo (needed if addCellNoise and filterCellNoise are both set)
  mutable std::unordered_map<uint64_t, double> m_emptyCellsMap;
};

#endif /* RECCALORIMETER_CREATECALOCELLS_H */
