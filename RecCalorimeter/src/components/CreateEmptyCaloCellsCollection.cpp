#include "CreateEmptyCaloCellsCollection.h"

// datamodel
#include "edm4hep/CalorimeterHitCollection.h"

DECLARE_COMPONENT(CreateEmptyCaloCellsCollection)

CreateEmptyCaloCellsCollection::CreateEmptyCaloCellsCollection(const std::string& name, ISvcLocator* svcLoc) : Gaudi::Algorithm(name, svcLoc) {
  declareProperty("cells", m_caloCells, "Empty calorimeter cells output)");
}

StatusCode CreateEmptyCaloCellsCollection::initialize() {
  StatusCode sc = Gaudi::Algorithm::initialize();
  if (sc.isFailure()) return sc;
  info() << "Create dummy cells collection initialized" << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode CreateEmptyCaloCellsCollection::execute(const EventContext&) const { 
  //create a new empty calo cells collection
  m_caloCells.createAndPut();
  return StatusCode::SUCCESS; 
}

StatusCode CreateEmptyCaloCellsCollection::finalize() { return Gaudi::Algorithm::finalize(); }
