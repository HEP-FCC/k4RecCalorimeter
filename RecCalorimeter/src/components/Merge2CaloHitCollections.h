#ifndef RECCALORIMETER_MERGE_2_CALO_HIT_COLLECTIONS_H
#define RECCALORIMETER_MERGE_2_CALO_HIT_COLLECTIONS_H

// k4FWCore
#include "k4FWCore/DataHandle.h"

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"

// Datamodel
#include "edm4hep/CalorimeterHitCollection.h"

class IGeoSvc;

/** @class Merge2CaloHitCollections
 *
 *  Algorithm for merging two collections together
 *  Use-case:
 *  Merging multiple calorimeter collections together.
 *  (ECAL barrel, HCAL barrel+extended barrel, ECAL + HCAL endcaps,
 *  ECAL + HCAL forward calorimeter). 
 *
 *  @author Juraj Smiesko
 *  @date   2023-11-06
 *
 */

class Merge2CaloHitCollections : public GaudiAlgorithm {

public:
  Merge2CaloHitCollections(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  /// Handle for the first hit collection (input)
  DataHandle<edm4hep::CalorimeterHitCollection> m_inCollA{
    "inHitsA", Gaudi::DataHandle::Reader, this};

  /// Handle for the second hit collection (input)
  DataHandle<edm4hep::CalorimeterHitCollection> m_inCollB{
    "inHitsB", Gaudi::DataHandle::Reader, this};

  /// Handle for the resulting hit collection (output)
  DataHandle<edm4hep::CalorimeterHitCollection> m_outColl{
    "outHits", Gaudi::DataHandle::Writer, this};
};

#endif /* RECCALORIMETER_MERGE_2_CALO_HIT_COLLECTIONS_H */
