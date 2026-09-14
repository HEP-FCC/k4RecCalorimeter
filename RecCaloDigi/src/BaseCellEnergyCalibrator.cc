#include "BaseCellEnergyCalibrator.h"

#include <edm4hep/CaloHitSimCaloHitLinkCollection.h>
#include <edm4hep/MutableCalorimeterHit.h>
#include <edm4hep/SimCalorimeterHit.h>

#include <exception>
#include <string>

using namespace std;

BaseCellEnergyCalibrator::BaseCellEnergyCalibrator(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc, {KeyValue("inputLinkCollection", "CaloHitLinks")},
                       {KeyValue("outputHitCollection", "CalorimeterHitsRec"),
                        KeyValue("outputRelationCollection", "CaloHitLinksRec")}) {}

StatusCode BaseCellEnergyCalibrator::initialize() {
  m_geoSvc = serviceLocator()->service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to retrieve the GeoSvc" << endmsg;
    return StatusCode::FAILURE;
  }

  // Validate the calibration configuration. These were asserts, which compile away in the
  // release builds this package actually ships in; as initialize() checks they also get to
  // tell the user what is wrong.
  if (m_calibrCoeff.empty()) {
    error() << "No calibration coefficients given: calibration_factorsMipGev is empty" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_calibrCoeff.size() != m_calLayers.size()) {
    error() << "calibration_factorsMipGev has " << m_calibrCoeff.size() << " entries but calibration_layergroups has "
            << m_calLayers.size() << ": they must describe the same groups" << endmsg;
    return StatusCode::FAILURE;
  }
  for (std::size_t k = 0; k < m_calibrCoeff.size(); ++k) {
    if (m_calibrCoeff[k] <= 0.f) {
      error() << "calibration_factorsMipGev[" << k << "] is " << m_calibrCoeff[k] << ", must be > 0" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_calLayers[k] <= 0) {
      error() << "calibration_layergroups[" << k << "] is " << m_calLayers[k]
              << ", must be > 0: a group spanning no layers has an unusable calibration constant" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // Expand the layer groups into one coefficient per layer. m_calLayers holds group sizes, so
  // group k covers the next m_calLayers[k] layers after the groups before it.
  m_layerCalib.clear();
  for (std::size_t k = 0; k < m_calLayers.size(); ++k) {
    m_layerCalib.insert(m_layerCalib.end(), static_cast<std::size_t>(m_calLayers[k]), m_calibrCoeff[k]);
  }
  debug() << "Calibration covers " << m_layerCalib.size() << " layers in " << m_calLayers.size() << " groups" << endmsg;

  // the cellID encoding does not change during the job, so decode it once here
  m_bitFieldCoder = dd4hep::DDSegmentation::BitFieldCoder(m_geoSvc->constantAsString(m_encodingStringVariable.value()));

  // both algorithms decode the layer number out of the cellID, so the encoding must provide it.
  // Catching this here turns a per-event exception into a clear configuration error.
  try {
    m_bitFieldCoder.index("layer");
  } catch (const std::exception&) {
    error() << "The encoding string '" << m_encodingStringVariable.value() << "' ("
            << m_bitFieldCoder.fieldDescription() << ") has no \"layer\" field" << endmsg;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

//-----------------------------------------------------------------------------------------------

std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>
BaseCellEnergyCalibrator::operator()(const edm4hep::CaloHitSimCaloHitLinkCollection& inputLinks) const {
  // * Reading Collections of digitised calorimeter Hits *
  edm4hep::CalorimeterHitCollection newcol;
  edm4hep::CaloHitSimCaloHitLinkCollection relcol;
  debug() << " number of elements = " << inputLinks.size() << endmsg;

  for (const auto& link : inputLinks) {
    const edm4hep::CalorimeterHit hit = link.get<edm4hep::CalorimeterHit>();
    edm4hep::MutableCalorimeterHit calhit = newcol.create(); // make new hit

    const auto cellID = hit.getCellID();
    float energy =
        reconstructEnergy(hit, m_bitFieldCoder.get(cellID, "layer")); // overloaded method, technology dependent

    calhit.setCellID(cellID);
    calhit.setEnergy(energy);
    calhit.setTime(hit.getTime());
    calhit.setPosition(hit.getPosition());
    calhit.setType(hit.getType());

    // the setters stay positional: podio's templated set() rejects the mutable handles that
    // create() returns, and the argument type makes the direction unambiguous here anyway
    auto newLink = relcol.create();
    newLink.setFrom(calhit);
    newLink.setTo(link.get<edm4hep::SimCalorimeterHit>());
    newLink.setWeight(1.0);
  }

  return std::make_tuple(std::move(newcol), std::move(relcol));
}

float BaseCellEnergyCalibrator::getLayerCalib(int ilayer) const {
  // retrieve calibration constants
  if (ilayer < 0 || static_cast<std::size_t>(ilayer) >= m_layerCalib.size()) {
    warning() << "Layer " << ilayer << " lies outside the calibrated range [0, " << m_layerCalib.size()
              << "): calibrating this hit to zero energy" << endmsg;
    return 0.f;
  }
  return m_layerCalib[ilayer];
}
