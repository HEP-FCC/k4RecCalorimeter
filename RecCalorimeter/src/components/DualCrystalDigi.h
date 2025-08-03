#ifndef DUAL_TEST_BEAM_DualCrysCalDigi_H
#define DUAL_TEST_BEAM_DualCrysCalDigi_H

#include "FWCore/Algorithm.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "DetInterface/IGeoSvc.h"

namespace DUAL {

class DualCrysCalDigi : public GaudiAlgorithm {
public:
  DualCrysCalDigi(const std::string& name, ISvcLocator* svcLoc);
  virtual StatusCode initialize() override;
  virtual StatusCode execute() override;
  virtual StatusCode finalize() override;

private:
  Gaudi::Property<std::string> m_simHitsA{"simHitsA","SimHitsCrystalA","Sim hits for crystal A"};
  Gaudi::Property<std::string> m_simHitsB{"simHitsB","SimHitsCrystalB","Sim hits for crystal B"};
  Gaudi::Property<std::string> m_outputHits{"outputHits","DualCrysCalDigi","Digitized output"};

  Gaudi::Property<double> m_gainA{"gainA",1000.0};
  Gaudi::Property<double> m_gainB{"gainB",1000.0};
  Gaudi::Property<double> m_pedestal{"pedestal",0.0};
  Gaudi::Property<double> m_noise{"noise",1.0};
  Gaudi::Property<uint32_t> m_adcBits{"adcBits",12};

  SmartIF<IGeoSvc> m_geoSvc;

  float digitizeEnergy(double energy, double gain) const;
};

} // namespace DUAL

#endif
