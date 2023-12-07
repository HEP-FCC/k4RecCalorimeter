#include "Merge2CaloHitCollections.h"
#include <GaudiKernel/MsgStream.h>


DECLARE_COMPONENT(Merge2CaloHitCollections)


Merge2CaloHitCollections::Merge2CaloHitCollections(
    const std::string& name,
    ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
  declareProperty("inHitsA", m_inCollA, "First calo hit collection (input)");
  declareProperty("inHitsB", m_inCollB, "Second calo hit collection (input)");
  declareProperty("outHits", m_outColl, "Resulting hit collection (output)");
}


StatusCode Merge2CaloHitCollections::initialize() {
  {
    StatusCode sc = GaudiAlgorithm::initialize();
    if (sc.isFailure()) {
      return sc;
    }
  }

  debug() << "Merge 2 calo hit collections initialized!" << endmsg;

  return StatusCode::SUCCESS;
}


StatusCode Merge2CaloHitCollections::execute() {
  // Create a new empty calo hits collection
  edm4hep::CalorimeterHitCollection* outColl = m_outColl.createAndPut();

  // Clone hits from the first collection
  const edm4hep::CalorimeterHitCollection* inCollA = m_inCollA.get();
  debug() << "Cloning " << inCollA->size() << " from the \""
          << m_inCollA.fullKey().key() << "\" collection." << endmsg;
  for (const auto& hit : *inCollA) {
    outColl->push_back(hit.clone());
  }

  // Clone hits from the Second collection
  const edm4hep::CalorimeterHitCollection* inCollB = m_inCollB.get();
  debug() << "Cloning " << inCollB->size() << " from the \""
          << m_inCollB.fullKey().key() << "\" collection." << endmsg;
  for (const auto& hit : *inCollB) {
    outColl->push_back(hit.clone());
  }

  debug() << "Output collection \"" << m_outColl.fullKey().key()
          << "\" contains " << outColl->size() << " calo hits." << endmsg;
  return StatusCode::SUCCESS;
}


StatusCode Merge2CaloHitCollections::finalize() {
  return GaudiAlgorithm::finalize();
}
