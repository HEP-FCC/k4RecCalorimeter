#ifndef REC_CALORIMETER_CaloDIGI_H
#define REC_CALORIMETER_CaloDIGI_H

#include "Gaudi/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartIF.h"
#include "GaudiKernel/DataHandle.h"
#include "k4FWCore/DataHandle.h"

#include "k4Interface/IGeoSvc.h"
#include "DDSegmentation/Segmentation.h"

#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"

namespace RECAL {

class CaloDIGI : public Gaudi::Algorithm {
public:
  CaloDIGI(const std::string& name, ISvcLocator* svcLoc);
  StatusCode initialize();
  StatusCode execute(const EventContext&) const;
  StatusCode finalize();


private:
  // Input simulated calorimeter hits handle
  mutable k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection> m_inputSimHits{
    "SimCalorimeterHits", Gaudi::DataHandle::Reader, this};

  // Output digitized calorimeter hits handle
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_outputHits{
    "DigitizedCalorimeterHits", Gaudi::DataHandle::Writer, this};

  // Digitization parameters
  Gaudi::Property<double> m_adcGain{
      this, "adcGain", 1000.0, "ADC conversion gain (ADC/GeV)"};

  Gaudi::Property<double> m_pedestal{
      this, "pedestal", 0.0, "ADC pedestal"};

  Gaudi::Property<double> m_noiseRMS{
      this, "noiseRMS", 0.0, "Gaussian ADC noise (ADC units)"};

  Gaudi::Property<double> m_timeSmear{
      this, "timeSmear", 0.0, "Gaussian time smearing (ns)"};

  Gaudi::Property<double> m_energyThreshold{
      this, "energyThreshold", 0.0, "Minimum energy to keep hit (MeV)"};

  Gaudi::Property<uint32_t> m_adcBits{
      this, "adcBits", 12, "ADC resolution in bits"};

  Gaudi::Property<bool> m_verbose{
      this, "verbose", false, "Print detailed digitization info"};

  SmartIF<IGeoSvc> m_geoSvc;

  // Internal helper functions
  double digitizeToMeV(double energyGeV) const;
  double smearTime(double time) const;
  double calculateTOF(double radius) const;
};

} // namespace RECAL

#endif

