#ifndef RECCALORIMETER_CREATEEMPTYCALOCELLSCOLLECTION_H
#define RECCALORIMETER_CREATEEMPTYCALOCELLSCOLLECTION_H

// k4FWCore
#include "k4FWCore/DataHandle.h"

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"

// edm4hep
#include "edm4hep/CalorimeterHitCollection.h"

class IGeoSvc;

/** @class CreateEmptyCaloCellsCollection
 *
 *  Algorithm for creating empty calorimeter cells collection
 *  Use-case:
 *  Input for clustering algorithm are collections from ALL calorimeter systems 
 *  (ECAL barrel, HCAL barrel+extended barrel, ECAL + HCAL endcaps, ECAL + HCAL forward calorimeter). 
 *  If not all collections are available (e.g. running simulations with just a subset of calorimeters), 
 *  the user is supposed to create dummy cell collections for the other calorimeters.
 *
 *  @author Jana Faltova
 *  @date   2017-08
 *
 */

class CreateEmptyCaloCellsCollection : public GaudiAlgorithm {

public:
  CreateEmptyCaloCellsCollection(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  /// Handle for the calo cells (output collection)
  DataHandle<edm4hep::CalorimeterHitCollection> m_caloCells{"cells", Gaudi::DataHandle::Writer, this};
};

#endif /* RECCALORIMETER_CREATEEMPTYCALOCELLSCOLLECTION_H */
