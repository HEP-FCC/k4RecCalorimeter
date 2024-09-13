#ifndef RECCALORIMETER_CREATEPOSITIONEDCALOCELLS_H
#define RECCALORIMETER_CREATEPOSITIONEDCALOCELLS_H

// k4FWCore
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/MetaDataHandle.h"
#include "k4Interface/ICellPositionsTool.h"
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

/** @class CreatePositionedCaloCells
 *
 *  Algorithm for creating positioned calorimeter cells from Geant4 hits.
 *  Based on CreateCaloCells, adding a positioning tool so that output
 *  cells from the digitisation contain the correct position in 3D space.
 * 
 *  Flow of the program:
 *  1/ Merge Geant4 energy deposits with same cellID
 *  2/ Emulate cross-talk (if switched on)
 *  3/ Calibrate to electromagnetic scale (if calibration switched on)
 *  4/ Add random noise to each cell (if noise switched on)
 *  5/ Filter cells and remove those with energy below threshold (if noise +
 *     filtering switched on)
 *  6/ Add cell positions
 *
 *  Tools called:
 *    - CalibrateCaloHitsTool
 *    - NoiseCaloCellsTool
 *    - CaloReadCrosstalkMap
 *    - CalorimeterTool (for geometry)
 *    - CellPositionsTool
 *
 *  @author Giovanni Marchiori
 *  @date   2024-07
 *
 */

class CreatePositionedCaloCells : public Gaudi::Algorithm {

public:
  CreatePositionedCaloCells(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize();

  StatusCode execute(const EventContext&) const;

  StatusCode finalize();

private:
  /// Handle for tool to get cells positions
  ToolHandle<ICellPositionsTool> m_cellPositionsTool{"CellPositionsTool", this};
  /// Handle for the calorimeter cells crosstalk tool
  mutable ToolHandle<ICaloReadCrosstalkMap> m_crosstalkTool{"ReadCaloCrosstalkMap", this};
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

  /// Handle for calo hits (input collection)
  mutable DataHandle<edm4hep::SimCalorimeterHitCollection> m_hits{"hits", Gaudi::DataHandle::Reader, this};
  /// Handle for the cellID encoding string of the input collection
  MetaDataHandle<std::string> m_hitsCellIDEncoding{m_hits, edm4hep::labels::CellIDEncoding, Gaudi::DataHandle::Reader};
  /// Handle for calo cells (output collection)
  mutable DataHandle<edm4hep::CalorimeterHitCollection> m_cells{"cells", Gaudi::DataHandle::Writer, this};
  MetaDataHandle<std::string> m_cellsCellIDEncoding{m_cells, edm4hep::labels::CellIDEncoding, Gaudi::DataHandle::Writer};
  /// Maps of cell IDs (corresponding to DD4hep IDs) on final energies to be used for clustering
  mutable std::unordered_map<uint64_t, double> m_cellsMap;
  /// Maps of cell IDs (corresponding to DD4hep IDs) on transfer of signals due to crosstalk
  mutable std::unordered_map<uint64_t, double> m_crosstalkCellsMap;
  /// Maps of cell IDs with zero energy, for all cells in calo (needed if addNoise and filterNoise are both set)
  mutable std::unordered_map<uint64_t, double> m_emptyCellsMap;
  /// Cache position vs cellID
  mutable std::unordered_map<dd4hep::DDSegmentation::CellID, edm4hep::Vector3f> m_positions_cache{};
};

#endif /* RECCALORIMETER_CREATEPOSITIONEDCALOCELLS_H */
