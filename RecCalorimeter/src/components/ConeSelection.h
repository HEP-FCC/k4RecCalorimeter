#ifndef RECCALORIMETER_CONESELECTION_H
#define RECCALORIMETER_CONESELECTION_H

// Gaudi
#include "Gaudi/Algorithm.h"

// Key4HEP
#include "k4FWCore/DataHandle.h"
#include "k4Interface/ICellPositionsTool.h"

// EDM4HEP
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/MCParticleCollection.h"

class IGeoSvc;

/** @class ConeSelection
 *
 *  Algorithm select cells within a cone around the generated particles.
 *
 *  @author Coralie Neubueser
 *  @date   2018-11
 *
 */

class ConeSelection : public Gaudi::Algorithm {

public:
  ConeSelection(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize();

  StatusCode execute(const EventContext&) const;

  StatusCode finalize();

private:
  /// Handle for tool to get cell positions
  ToolHandle<ICellPositionsTool> m_cellPositionsTool{"CellPositionsTool", this};

  /// Handle for calo hits (input collection)
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_cells{"cells", Gaudi::DataHandle::Reader, this};
  /// Handle for calo hits (input collection)
  mutable k4FWCore::DataHandle<edm4hep::MCParticleCollection> m_particles{"particles", Gaudi::DataHandle::Reader, this};
  /// Handle for calo cells (output collection)
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_selCells{"selCells", Gaudi::DataHandle::Writer,
                                                                             this};
  /// Map of cell IDs (corresponding to DD4hep IDs) and energy
  mutable std::unordered_map<uint64_t, double> m_cellsMap;

  Gaudi::Property<double> m_r{this, "radius", 0.4, "radius of selection cone"};
};

#endif /* RECCALORIMETER_CONESELECTION_H */
