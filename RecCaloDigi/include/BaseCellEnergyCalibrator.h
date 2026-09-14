#ifndef BASECELLENERGYCALIBRATOR_H
#define BASECELLENERGYCALIBRATOR_H 1

#include <k4FWCore/Transformer.h>

#include <edm4hep/CaloHitSimCaloHitLinkCollection.h>
#include <edm4hep/CalorimeterHit.h>
#include <edm4hep/CalorimeterHitCollection.h>

#include "k4Interface/IGeoSvc.h"

#include <DDSegmentation/BitFieldCoder.h>

#include <string>
#include <vector>

/** === BaseCellEnergyCalibrator algorithm === <br>
    calibration of digitised calorimeter hits
    e.g. apply sampling fraction correction
    abstract base class, technology independent
    D.Jeans 02/2016.

    24 March 2016: removed gap corrections - to be put into separate processor
    changed relations: now keep relation between reconstructed and simulated hits.
*/

struct BaseCellEnergyCalibrator
    : k4FWCore::MultiTransformer<
          std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>(
              const edm4hep::CaloHitSimCaloHitLinkCollection&)> {

public:
  BaseCellEnergyCalibrator(const std::string& name, ISvcLocator* svcLoc);
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  StatusCode initialize();

  /** Called for every event.
   */
  std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>
  operator()(const edm4hep::CaloHitSimCaloHitLinkCollection& inputLinks) const;

protected:
  float getLayerCalib(int ilayer) const;
  virtual float reconstructEnergy(const edm4hep::CalorimeterHit& hit,
                                  int layer) const = 0; // to be overloaded, technology-specific

  // parameters
  // Grouping of calo layers
  Gaudi::Property<std::vector<int>> m_calLayers{this, "calibration_layergroups", {}, "Grouping of calo layers"};
  // Calibration coefficients for layers groups
  Gaudi::Property<std::vector<float>> m_calibrCoeff{
      this, "calibration_factorsMipGev", {}, "Calibration coefficients (MIP->shower GeV) of layers groups"};
  // name of the DD4hep constant holding the calorimeter cellID encoding
  Gaudi::Property<std::string> m_encodingStringVariable{
      this, "EncodingStringParameterName", "GlobalCalorimeterReadoutID",
      "The name of the DD4hep constant that contains the encoding string for calorimeters"};

  SmartIF<IGeoSvc> m_geoSvc;
  // built once in initialize() from the geometry encoding string
  dd4hep::DDSegmentation::BitFieldCoder m_bitFieldCoder{};
  // calibration coefficient per layer, expanded once in initialize() from the layer groups
  // so that the per-hit lookup is a bounds check and an index
  std::vector<float> m_layerCalib{};
};

#endif
