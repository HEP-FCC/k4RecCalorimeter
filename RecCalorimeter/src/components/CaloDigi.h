#ifndef REC_CALORIMETER_CaloDIGI_H
#define REC_CALORIMETER_CaloDIGI_H

#include "FWCore/Algorithm.h"
#include "DetInterface/IGeoSvc.h"
#include "GaudiKernel/ToolHandle.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "DDSegmentation/Segmentation.h"

namespace RECAL {

class CaloDIGI : public GaudiAlgorithm {
public:
  CaloDIGI(const std::string& name, ISvcLocator* svcLoc);
  virtual StatusCode initialize() override;
  virtual StatusCode execute() override;
  virtual StatusCode finalize() override;

private:
  Gaudi::Property<std::string> m_inputSimHits{"inputSimHits", "SimCalorimeterHits", "Input simulated hits"};
  Gaudi::Property<std::string> m_outputHits{"outputHits", "DigitizedCalorimeterHits", "Output digitized hits"};

  Gaudi::Property<double> m_adcGain{"adcGain", 1000.0, "ADC conversion gain (ADC/GeV)"};
  Gaudi::Property<double> m_pedestal{"pedestal", 0.0, "ADC pedestal"};
  Gaudi::Property<double> m_noiseRMS{"noiseRMS", 0.0, "Gaussian ADC noise (ADC units)"};
  Gaudi::Property<double> m_timeSmear{"timeSmear", 0.0, "Gaussian time smearing (ns)"};
  Gaudi::Property<double> m_energyThreshold{"energyThreshold", 0.0, "Minimum energy to keep hit (MeV)"};
  Gaudi::Property<uint32_t> m_adcBits{"adcBits", 12, "ADC resolution in bits"};
  Gaudi::Property<bool> m_verbose{"verbose", false, "Print detailed digitization info"};

  SmartIF<IGeoSvc> m_geoSvc;

  double digitizeToMeV(double energyGeV) const;
  double smearTime(double time) const;
  double calculateTOF(double radius) const;
};

} // namespace RECAL

#endif
