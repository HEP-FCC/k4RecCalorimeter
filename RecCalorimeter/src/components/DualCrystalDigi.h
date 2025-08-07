#ifndef DUAL_TEST_BEAM_DUALCRYSTALDIGI_H
#define DUAL_TEST_BEAM_DUALCRYSTALDIGI_H

#include "Gaudi/Algorithm.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "k4Interface/IGeoSvc.h"
#include "k4FWCore/DataHandle.h"

namespace DUAL {

class DualCrystalDigi : public Gaudi::Algorithm {
public:
  DualCrystalDigi(const std::string& name, ISvcLocator* svcLoc);
  StatusCode initialize() override;
  StatusCode execute(const EventContext& ctx) const override;
  StatusCode finalize() override;

private:
  mutable k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection> m_simHitsA{
    "SimCalorimeterHitsA", Gaudi::DataHandle::Reader, this};
  mutable k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection> m_simHitsB{
    "SimHitsB", Gaudi::DataHandle::Reader, this};
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_outputHits{
    "OutputHits", Gaudi::DataHandle::Writer, this};

  Gaudi::Property<double> m_gainA{"gainA",1000.0};
  Gaudi::Property<double> m_gainB{"gainB",1000.0};
  Gaudi::Property<double> m_pedestal{"pedestal",0.0};
  Gaudi::Property<double> m_noise{"noise",1.0};
  Gaudi::Property<uint32_t> m_adcBits{"adcBits",12};

  SmartIF<IGeoSvc> m_geoSvc;

  float digitizeEnergy(double energy, double gain) const;
};

} // namespace DUAL

#endif // DUAL_TEST_BEAM_DUALCRYSTALDIGI_H
