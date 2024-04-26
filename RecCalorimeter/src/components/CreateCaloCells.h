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
#include "GaudiAlg/GaudiAlgorithm.h"
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
 *  2/ Calibrate to electromagnetic scale (if calibration switched on)
 *  3/ Add random noise to each cell (if noise switched on)
 *  4/ Filter cells and remove those with energy below threshold (if noise +
 * filtering switched on)
 *
 *  Tools called:
 *    - CalibrateCaloHitsTool
 *    - NoiseCaloCellsTool
 *
 *  @author Jana Faltova
 *  @author Anna Zaborowska
 *  @date   2016-09
 *
 */

class CreateCaloCells : public GaudiAlgorithm {

public:
  CreateCaloCells(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:

  /// Handle for the calorimeter cells crosstalk tool
  ToolHandle<ICaloReadCrosstalkMap> m_crosstalksTool{"ReadCaloCrosstalkMap", this};
  /// Handle for tool to calibrate Geant4 energy to EM scale tool
  ToolHandle<ICalibrateCaloHitsTool> m_calibTool{"CalibrateCaloHitsTool", this};
  /// Handle for the calorimeter cells noise tool
  ToolHandle<INoiseCaloCellsTool> m_noiseTool{"NoiseCaloCellsFlatTool", this};
  /// Handle for the geometry tool
  ToolHandle<ICalorimeterTool> m_geoTool{"TubeLayerPhiEtaCaloTool", this};

  /// Add crosstalk to cells?
  Gaudi::Property<bool> m_addCrosstalk{this, "addCrosstalk", true, "Add crosstalk effect?"};
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
  DataHandle<edm4hep::SimCalorimeterHitCollection> m_hits{"hits", Gaudi::DataHandle::Reader, this};
  /// Handle for the cellID encoding string of the input collection
  MetaDataHandle<std::string> m_hitsCellIDEncoding{m_hits, edm4hep::CellIDEncoding, Gaudi::DataHandle::Reader};
  /// Handle for calo cells (output collection)
  DataHandle<edm4hep::CalorimeterHitCollection> m_cells{"cells", Gaudi::DataHandle::Writer, this};
  MetaDataHandle<std::string> m_cellsCellIDEncoding{m_cells, edm4hep::CellIDEncoding, Gaudi::DataHandle::Writer};
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

  /** Temporary: for use with MergeLayer tool
   * MergeLayer is going to be replaced by RedoSegmentation once we can define
   * segmentation with variable cell (layer) size.
   * This property won't be needed anymore.
   */
  unsigned int m_activeVolumesNumber;
  /// Use only volume ID? If false, using PhiEtaSegmentation
  bool m_useVolumeIdOnly;

  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;
  dd4hep::VolumeManager m_volman;
  /// Maps of cell IDs (corresponding to DD4hep IDs) on (1) final energies to be used for clustering, (2) energy desposits from EM shower hits and (3) transfer of signals due to crosstalk
  std::unordered_map<uint64_t, double> m_cellsMap;
  std::unordered_map<uint64_t, double> m_HitCellsMap;
  std::unordered_map<uint64_t, double> m_CrosstalkCellsMap;
};

#endif /* RECCALORIMETER_CREATECALOCELLS_H */
