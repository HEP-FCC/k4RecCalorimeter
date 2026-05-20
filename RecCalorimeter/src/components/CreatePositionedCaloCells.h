#ifndef RECCALORIMETER_CREATEPOSITIONEDCALOCELLS_H
#define RECCALORIMETER_CREATEPOSITIONEDCALOCELLS_H

// k4FWCore
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/Transformer.h"

#include <tuple>

// Interfaces
#include "RecCaloCommon/ICalibrateCaloHitsTool.h"
#include "RecCaloCommon/ICaloReadCrosstalkMap.h"
#include "RecCaloCommon/ICalorimeterTool.h"
#include "RecCaloCommon/ICellPositionsTool.h"
#include "RecCaloCommon/INoiseCaloCellsTool.h"

// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"

// edm4hep
#include "edm4hep/CaloHitSimCaloHitLinkCollection.h"
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

struct CreatePositionedCaloCells final
    : k4FWCore::MultiTransformer<
          std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>(
              const edm4hep::SimCalorimeterHitCollection&)> {

  CreatePositionedCaloCells(const std::string& name, ISvcLocator* svcLoc)
      : MultiTransformer(name, svcLoc, {KeyValue("hits", "")}, /// Hits from which to create cells (input)
                         {
                             KeyValue("cells", ""), /// The created calorimeter cells (output)
                             KeyValue("links", ""), /// The links between hits and cells (output)
                         }) {
    declareProperty("positionsTool", m_cellPositionsTool, "Handle for cell positions tool");
    declareProperty("crosstalkTool", m_crosstalkTool, "Handle for the cell crosstalk tool");
    declareProperty("calibTool", m_calibTool, "Handle for tool to calibrate Geant4 energy to EM scale tool");
    declareProperty("noiseTool", m_noiseTool, "Handle for the calorimeter cells noise tool");
    declareProperty("geometryTool", m_geoTool, "Handle for the geometry tool");
    m_decoder = nullptr;
  }

  virtual ~CreatePositionedCaloCells();

  StatusCode initialize() override;
  std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>
  operator()(const edm4hep::SimCalorimeterHitCollection& hits) const override;

private:
  /// Handle for tool to get cells positions
  ToolHandle<k4::recCalo::ICellPositionsTool> m_cellPositionsTool{"CellPositionsTool", this};
  /// Handle for the calorimeter cells crosstalk tool
  mutable ToolHandle<k4::recCalo::ICaloReadCrosstalkMap> m_crosstalkTool{"ReadCaloCrosstalkMap", this};
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

  /// Name for calo hits (input collection)
  mutable std::string m_hits;
  /// Name for calo cells (output collection)
  mutable std::string m_cells;

  /// Maps of cell IDs (corresponding to DD4hep IDs) vs digitised cell energies
  mutable std::unordered_map<uint64_t, double> m_cellsMap;
  /// Maps of cell IDs (corresponding to DD4hep IDs) on transfer of signals due to crosstalk
  mutable std::unordered_map<uint64_t, double> m_crosstalkCellsMap;
  /// Maps of cell IDs with zero energy, for all cells in calo (needed if addNoise and filterNoise are both set)
  mutable std::unordered_map<uint64_t, double> m_emptyCellsMap;
  /// Cache position vs cellID
  mutable std::unordered_map<dd4hep::DDSegmentation::CellID, edm4hep::Vector3f> m_positions_cache{};

  /// For cell type - for PandoraPFA
  mutable int m_calotype;
  mutable int m_caloid;
  mutable int m_layout;
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
};

#endif /* RECCALORIMETER_CREATEPOSITIONEDCALOCELLS_H */
