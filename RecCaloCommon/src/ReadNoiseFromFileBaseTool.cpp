/**
 * @file ReadNoiseFromFileBaseTool.cpp
 * @author scott snyder <snyder@bnl.gov>, Jana Faltova, Coralie Neubueser, Giovanni Marchior
 * @date Jul, 2026
 * @brief Base class for reading noise data taken from ROOT histograms
 *        read from a file.
 */

#include "RecCaloCommon/ReadNoiseFromFileBaseTool.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"

// k4FWCore
#include "k4FWCore/GaudiChecks.h"
#include "k4Interface/IGeoSvc.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DD4hep/Readout.h"
#include "DDSegmentation/Segmentation.h"

// ROOT
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TSystem.h"

#include <format>

StatusCode ReadNoiseFromFileBaseTool::initialize() {
  K4_GAUDI_CHECK(AlgTool::initialize());

  K4_GAUDI_CHECK(m_geoSvc.retrieve());
  K4_GAUDI_CHECK(m_indexerSvc.retrieve());
  K4_GAUDI_CHECK(m_constantsSvc.retrieve());

  if (m_numRadialLayers > 0) {
    m_numHistograms = m_numRadialLayers;
    warning() << "numRadialLayers property is deprecated; use numHistograms instead." << endmsg;
  }

  // Take readout bitfield decoder from GeoSvc
  dd4hep::Readout readout = m_geoSvc->getDetector()->readout(m_readoutName);
  m_decoder = readout.idSpec().decoder();

  m_index_activeField = m_decoder->index(m_activeFieldName);

  int detID = readout.segmentation().detector()->id;
  m_indexer = m_indexerSvc->indexer(detID);
  K4_GAUDI_CHECK(m_indexer != nullptr);

  // get noise constants.
  m_data = m_constantsSvc->getObj<NoiseData>(m_noiseFileName);
  if (!m_data) {
    NoiseData noise;
    K4_GAUDI_CHECK(ReadNoiseFromFileBaseTool::initNoiseFromFile(noise));
    K4_GAUDI_CHECK(m_constantsSvc->putObj(m_noiseFileName, std::move(noise)));
    m_data = m_constantsSvc->getObj<NoiseData>(m_noiseFileName);
    K4_GAUDI_CHECK(m_data != nullptr);
  }

  return StatusCode::SUCCESS;
}

/// Return noise RMS for a given CellID.
double ReadNoiseFromFileBaseTool::getNoiseRMSPerCell(CellID aCellId) const {
  // Get cell coordinates: bin and radial layer
  unsigned ndx = m_indexer->index(aCellId);
  int ibin = m_data->m_bins.at(ndx).first;

  unsigned cellLayer = m_decoder->get(aCellId, m_index_activeField);
  return getNoiseRMSPerCell(ibin, cellLayer);
}

/// Return noise offset for a given CellID.
double ReadNoiseFromFileBaseTool::getNoiseOffsetPerCell(CellID aCellId) const {

  if (!m_setNoiseOffset)
    return 0.;

  // Get cell coordinates: bin and radial layer
  unsigned ndx = m_indexer->index(aCellId);
  int ibin = m_data->m_bins.at(ndx).second;

  unsigned cellLayer = m_decoder->get(aCellId, m_index_activeField);
  return getNoiseOffsetPerCell(ibin, cellLayer);
}

/// Return both noise RMS (first) and offset (second)  for a given CellID.
std::pair<double, double> ReadNoiseFromFileBaseTool::getNoisePerCell(CellID aCellId) const {
  // Get cell coordinates: bin and radial layer
  unsigned ndx = m_indexer->index(aCellId);
  const std::pair<unsigned, unsigned>& bins = m_data->m_bins.at(ndx);

  unsigned cellLayer = m_decoder->get(aCellId, m_index_activeField);

  double rms = getNoiseRMSPerCell(bins.first, cellLayer);
  double offset = m_setNoiseOffset ? getNoiseOffsetPerCell(bins.second, cellLayer) : 0;
  return std::make_pair(rms, offset);
}

/// Return noise RMS for a given bin/layer.
double ReadNoiseFromFileBaseTool::getNoiseRMSPerCell(int ibin, unsigned cellLayer) const {
  // Check that there are not more layers than the constants are provided for
  if (cellLayer >= m_data->m_histoElecNoiseRMS.size()) {
    error() << "More radial layers than we have noise for!!!!" << endmsg;
    return 0;
  }

  double elecNoiseRMS = m_data->m_histoElecNoiseRMS.at(cellLayer)->GetBinContent(ibin);

  // Total noise: electronics noise + pileup
  double totalNoiseRMS = 0;
  if (m_addPileup) {
    double pileupNoiseRMS = m_data->m_histoPileupNoiseRMS.at(cellLayer)->GetBinContent(ibin);
    totalNoiseRMS = std::hypot(elecNoiseRMS, pileupNoiseRMS) * m_scaleFactor;
  } else { // avoid useless math operations if no pileup
    totalNoiseRMS = std::abs(elecNoiseRMS) * m_scaleFactor;
  }

  if (totalNoiseRMS < 1e-6) {
    warning() << std::format("Zero noise: cell bin {} layer {} noise {}", ibin, cellLayer, totalNoiseRMS) << endmsg;
  }

  return totalNoiseRMS;
}

/// Return noise offset for a given bin/layer.
double ReadNoiseFromFileBaseTool::getNoiseOffsetPerCell(int ibin, unsigned cellLayer) const {
  // Check that there are not more layers than the constants are provided for
  if (cellLayer >= m_data->m_histoElecNoiseOffset.size()) {
    error() << "More radial layers than we have noise for!!!!" << endmsg;
    return 0;
  }

  double elecNoiseOffset = m_data->m_histoElecNoiseOffset.at(cellLayer)->GetBinContent(ibin);

  // Total noise: electronics noise + pileup
  double pileupNoiseOffset = m_addPileup ? m_data->m_histoPileupNoiseOffset.at(cellLayer)->GetBinContent(ibin) : 0;
  double totalNoiseOffset = (elecNoiseOffset + pileupNoiseOffset) * m_scaleFactor;

  // No warning is printed if offset is zero because that is the usual scenario

  return totalNoiseOffset;
}

/// Open file and read noise histograms into a NoiseData structure.
StatusCode ReadNoiseFromFileBaseTool::initNoiseFromFile(NoiseData& data) const {
  // Check if file exists
  if (m_noiseFileName.empty()) {
    error() << "Name of the file with the noise values not provided!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (gSystem->AccessPathName(m_noiseFileName.value().c_str())) {
    error() << "Provided file with the noise values not found!" << endmsg;
    error() << "File path: " << m_noiseFileName.value() << endmsg;
    return StatusCode::FAILURE;
  }
  std::unique_ptr<TFile> noiseFile(TFile::Open(m_noiseFileName.value().c_str(), "READ"));
  if (noiseFile->IsZombie()) {
    error() << "Unable to open the file with the noise values!" << endmsg;
    error() << "File path: " << m_noiseFileName.value() << endmsg;
    return StatusCode::FAILURE;
  } else {
    info() << "Using the following file with the noise values: " << m_noiseFileName.value() << endmsg;
  }

  auto getHist = [&](const std::string& name) {
    debug() << "Getting histogram with a name " << name << endmsg;
    TH1* h = dynamic_cast<TH1*>(noiseFile->Get(name.c_str()));
    if (!h) {
      error() << "Cannot read noise histogram " << name << " from file " << m_noiseFileName << endmsg;
    } else if (h->GetNbinsX() < 1) {
      error() << std::format("Histogram {} has 0 bins! check the file with noise and the name of the histogram!", name)
              << endmsg;
      h = nullptr;
    } else {
      h->SetDirectory(nullptr);
    }

    return std::unique_ptr<TH1>(h);
  };

  auto getHists = [&](const std::string& nameRoot, std::vector<std::unique_ptr<TH1>>& hists) {
    for (unsigned i = 0; i < m_numHistograms; i++) {
      std::string name = std::format("{}{}", nameRoot, i + 1);
      hists.push_back(getHist(name));
      if (!hists.back())
        return false;
    }
    return true;
  };

  if (!getHists(m_elecNoiseHistoName, data.m_histoElecNoiseRMS))
    return StatusCode::FAILURE;
  if (m_setNoiseOffset) {
    if (!getHists(m_elecNoiseOffsetHistoName, data.m_histoElecNoiseOffset))
      return StatusCode::FAILURE;
  }
  if (m_addPileup) {
    if (!getHists(m_pileupHistoName, data.m_histoPileupNoiseRMS))
      return StatusCode::FAILURE;
    if (m_setNoiseOffset) {
      if (!getHists(m_pileupOffsetHistoName, data.m_histoPileupNoiseOffset))
        return StatusCode::FAILURE;
    }
  }

  K4_GAUDI_CHECK(initBinning(data, *m_indexer));

  return StatusCode::SUCCESS;
}

// Need a copy ctor to put these in std::any.
ReadNoiseFromFileBaseTool::NoiseData::NoiseData(const NoiseData& other) : m_bins(other.m_bins) {
  auto copyHists = [](std::vector<std::unique_ptr<TH1>>& dst, const std::vector<std::unique_ptr<TH1>>& src) {
    dst.reserve(src.size());
    for (size_t i = 0; i < src.size(); i++) {
      TH1* hcopy = dynamic_cast<TH1*>(src[i]->Clone());
      hcopy->SetDirectory(nullptr);
      dst.emplace_back(hcopy);
    }
  };

  copyHists(m_histoPileupNoiseRMS, other.m_histoPileupNoiseRMS);
  copyHists(m_histoElecNoiseRMS, other.m_histoElecNoiseRMS);
  copyHists(m_histoPileupNoiseOffset, other.m_histoPileupNoiseOffset);
  copyHists(m_histoElecNoiseOffset, other.m_histoElecNoiseOffset);
}
