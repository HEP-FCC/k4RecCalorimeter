#include "FilterDoubleLayerHits.h"
#include <iostream>
#include <cstdlib>
#include <string.h>

#include <edm4hep/TrackerHitPlane.h>
#include <GaudiKernel/GaudiException.h>

#include <DDSegmentation/BitFieldCoder.h>

#include "DD4hep/DD4hepUnits.h"

#include <TH1F.h>
#include <TH2F.h>

DECLARE_COMPONENT(FilterDoubleLayerHits)

FilterDoubleLayerHits::FilterDoubleLayerHits(const std::string& name, ISvcLocator* svcLoc) : Transformer(name, svcLoc,
    KeyValues("InputCollection", {"VXDTrackerHitPlanes"}),
    KeyValues("OutputCollection", {"VXDTrackerHitPlanes_DLFiltered"})) {}


dd4hep::rec::Vector2D FilterDoubleLayerHits::globalToLocal(long int cellID, const dd4hep::rec::Vector3D& posGlobal, dd4hep::rec::ISurface** surfptr=nullptr) const{
  dd4hep::rec::ISurface* surf;
  // Using directly the provided surface object if available
  if (surfptr && *surfptr) {
    surf = *surfptr;
  } else {
    // Finding the surface corresponding to the cellID
    dd4hep::rec::SurfaceMap::const_iterator surfIt = m_map->find( cellID );
    if( surfIt == m_map->end() ){
      throw GaudiException(" FilterDoubleLayerHits::processEvent(): no surface found for cellID: " + cellID, "Fail", StatusCode::FAILURE);
    }
    surf = surfIt->second;
    // Saving the surface object outside the function to be reused for the same cellID
    if (surfptr) *surfptr = surf;
  }
  // Converting global position to local in [cm]
  dd4hep::rec::Vector2D posLocal = surf->globalToLocal(  dd4hep::mm * posGlobal );

  return dd4hep::rec::Vector2D( posLocal.u() / dd4hep::mm, posLocal.v() / dd4hep::mm );
}


StatusCode FilterDoubleLayerHits::initialize() {
  m_geoSvc = serviceLocator()->service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to retrieve the GeoSvc" << endmsg;
    return StatusCode::FAILURE;
  }
  // Get Histogram and Data Services
	SmartIF<ITHistSvc> histSvc;
	histSvc = serviceLocator()->service("THistSvc");

  debug() << "   init called  " << endmsg;

  // Extracting double-layer cut configurations
  m_dlCuts.resize( m_dlCutConfigs.size() / 4 ) ;

  unsigned i=0,index=0 ;
  while( i < m_dlCutConfigs.size() ){
    m_dlCuts[index].layer0      = std::atoi( m_dlCutConfigs[ i++ ].c_str() ) ;
    m_dlCuts[index].layer1      = std::atoi( m_dlCutConfigs[ i++ ].c_str() ) ;
    m_dlCuts[index].dPhi_max    = std::atof( m_dlCutConfigs[ i++ ].c_str() ) / 1e3 ;  // converting mrad -> rad
    m_dlCuts[index].dTheta_max  = std::atof( m_dlCutConfigs[ i++ ].c_str() ) / 1e3 ;  // converting mrad -> rad
    ++index ;
  }

  // Get the surface map from the SurfaceManager
  dd4hep::Detector& theDetector = *(m_geoSvc->getDetector());
  dd4hep::rec::SurfaceManager& surfMan = *theDetector.extension<dd4hep::rec::SurfaceManager>();
  dd4hep::DetElement det = theDetector.detector( m_subDetName );

  m_map = surfMan.map( det.name() );

  if( !m_map ) {
    throw GaudiException( " Could not find surface map for detector: " + m_subDetName + " in SurfaceManager", "Fail", StatusCode::FAILURE);
  }

  // Booking diagnostic histograms for each configured cut
  char hname[100];
  for(size_t iCut=0, nCuts=m_dlCuts.size() ; iCut<nCuts ; ++iCut){
    const DoubleLayerCut& cut = m_dlCuts.at(iCut);
    if (m_fillHistos) {
      // Properties of the closest hits across 2 sublayers
      sprintf(hname, "h_dU_layers_%d_%d", cut.layer0, cut.layer1);
      m_histos[ std::string(hname) ] = new TH1F( hname , ";#DeltaU [mm]; Hit pairs", 2000, -5, 5 );
      sprintf(hname, "h2_dU_dPhi_layers_%d_%d", cut.layer0, cut.layer1);
      m_histos[ std::string(hname) ] = new TH2F( hname , ";#DeltaU [mm]; #Delta#phi [mrad]", 500, -5, 5, 500, 0, 50);
      sprintf(hname, "h_dTheta_layers_%d_%d", cut.layer0, cut.layer1);
      m_histos[ std::string(hname) ] = new TH1F( hname , ";#Delta#Theta [mrad]; Hit pairs", 3000, -30, 30 );
      sprintf(hname, "h_dPhi_layers_%d_%d", cut.layer0, cut.layer1);
      m_histos[ std::string(hname) ] = new TH1F( hname , ";|#Delta#phi| [mrad]; Hit pairs", 1000, 0, 50 );
      sprintf(hname, "h_dt_layers_%d_%d", cut.layer0, cut.layer1);
      m_histos[ std::string(hname) ] = new TH1F( hname , ";|#Deltat| [ns]; Hit pairs", 2000, -10, 10 );
      // Hit properties for each individual layer
      std::vector<unsigned int> layers{cut.layer0, cut.layer1};
      for (auto layer : layers) {
        sprintf(hname, "h2_posUV_rejected_layer_%d", layer);
        m_histos[ std::string(hname) ] = new TH2F( hname , ";U [mm]; V [mm]", 500, -100, 100, 1000, -200, 200 );
      }
    }

    for (const auto& [name, histo] : m_histos) {
        (void)m_histSvc->regHist("/histos/" + name, histo);
    }

    // Printing the configured cut
    debug() <<  iCut << ". layers: " << cut.layer0 << " >> " << cut.layer1 << ";  dPhi: "
    << cut.dPhi_max << " rad;  dTheta: " << cut.dTheta_max << " rad" << endmsg;
  }

  return StatusCode::SUCCESS;

}


edm4hep::TrackerHitPlaneCollection FilterDoubleLayerHits::operator()(const edm4hep::TrackerHitPlaneCollection& inputTrackerHitCollection) const{
  std::string initString;
  initString = m_geoSvc->constantAsString(m_encodingStringVariable.value());
  dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(initString);  // check!

  //---- create the output collection
  edm4hep::TrackerHitPlaneCollection outCol;
  outCol.setSubsetCollection();

  // Set acceptance flags for all hits to FALSE
  const size_t nHit = inputTrackerHitCollection.size(); 
  bool hitAccepted[NHITS_MAX]; ////Array of flags for hits to be accepted
  memset(&hitAccepted, false, nHit);

  ////Map of vectors of hits grouped by position in the detector
  std::map<SensorPosition, std::vector<size_t> > hitsGrouped{};

  // Splitting hits by sensor ids for faster association
  for (size_t iHit = 0; iHit < nHit ; iHit++) {

    edm4hep::TrackerHitPlane h = inputTrackerHitCollection.at( iHit );

    unsigned int layerID  = bitFieldCoder.get(h.getCellID(), "layer");
    unsigned int sideID   = bitFieldCoder.get(h.getCellID(), "side");
    unsigned int ladderID = bitFieldCoder.get(h.getCellID(), "module");
    unsigned int moduleID = bitFieldCoder.get(h.getCellID(), "sensor");

    SensorPosition sensPos = {layerID, sideID, ladderID, moduleID};
    if (hitsGrouped.find(sensPos) == hitsGrouped.end()) {
      hitsGrouped[sensPos] = std::vector<size_t>();
      hitsGrouped[sensPos].reserve(nHit);
    }
    hitsGrouped[sensPos].push_back(iHit);
  }

  //---- loop over hits
  for (size_t iHit = 0; iHit < nHit ; iHit++) {

    // Skipping hits that are already accepted
    if (hitAccepted[iHit]) continue;

    edm4hep::TrackerHitPlane h = inputTrackerHitCollection.at( iHit );

    unsigned int layerID  = bitFieldCoder.get(h.getCellID(), "layer");
    unsigned int sideID   = bitFieldCoder.get(h.getCellID(), "side");
    unsigned int ladderID = bitFieldCoder.get(h.getCellID(), "module");
    unsigned int moduleID = bitFieldCoder.get(h.getCellID(), "sensor");
    debug() << " Checking 1st hit " << iHit << " / " << nHit << " at layer: " << layerID << "  ladder: " << ladderID << "  module: " << moduleID <<  endmsg ;

    const SensorPosition sensPos = {layerID, sideID, ladderID, moduleID};

    // Checking if the hit is at the inner double layer to be filtered
    const DoubleLayerCut* dlCut(0);
    for (int iCut=0, nCuts=m_dlCuts.size(); iCut<nCuts; ++iCut) {
      const DoubleLayerCut& cut = m_dlCuts.at(iCut);
      if( ( layerID != cut.layer0 ) && ( layerID != cut.layer1 ) ) continue;
      dlCut = &cut;
      break;
    }

    // Accepting hit immediately if it's not affected by any double-layer cut
    if (!dlCut) {
      hitAccepted[iHit] = true;
      continue;
    }

    // Skipping the hit from the first pass if it belongs to the outer sublayer of the cut
    if (layerID == dlCut->layer1) continue;

    // Getting local and global hit positions
    dd4hep::rec::Vector3D posGlobal( h.getPosition().x, h.getPosition().y, h.getPosition().z );
    dd4hep::rec::Vector2D posLocal = globalToLocal( h.getCellID(), posGlobal );

    // Setting the values for closest hits
    double dR_min(999.0);
    double dU_closest(0.0);
    double dTheta_closest(0.0);
    double dPhi_closest(0.0);
    double dt_closest(0.0);

    // Looking for the compliment hits in the 2nd sublayer
    size_t nCompatibleHits(0);
    SensorPosition sensPos2 = sensPos;
    sensPos2.layer = dlCut->layer1;
    dd4hep::rec::ISurface* surf=nullptr;
    // Checking if there are any hits in the corresponding sensor at the other sublayer
    if (hitsGrouped.find(sensPos2) == hitsGrouped.end()) continue;
    for (size_t iHit2 : hitsGrouped.at(sensPos2)) {
      edm4hep::TrackerHitPlane h2 = inputTrackerHitCollection.at( iHit2 );
      unsigned int layerID2 = bitFieldCoder.get(h2.getCellID(), "layer");

      // Checking whether hit is in the time acceptance window
      double dt = h2.getTime() - h.getTime();
      if (m_dtMax >= 0.0 && std::fabs(dt) > m_dtMax) continue;

      // Getting the local and global hit positions
      dd4hep::rec::Vector3D posGlobal2( h2.getPosition().x, h2.getPosition().y, h2.getPosition().z );
      dd4hep::rec::Vector2D posLocal2 = globalToLocal( h2.getCellID(), posGlobal2, &surf );

      // Checking whether hit is close enough to the 1st one
      double dU = posLocal2.u() - posLocal.u();
      double dTheta = posGlobal2.theta() - posGlobal.theta();
      double dPhi = std::fabs(posGlobal2.phi() - posGlobal.phi());
      if (dPhi > dd4hep::pi) dPhi = dd4hep::twopi - dPhi;
      double dR = sqrt(dPhi*dPhi + dTheta*dTheta);
      debug() << " Checking 2nd hit at layer: " << layerID2 << ";  dPhi: " << dPhi << ";  dTheta: " <<  dTheta << endmsg;

      // Updating the minimal values
      if (dR < dR_min) {
        dR_min = dR;
        dU_closest = dU;
        dTheta_closest = dTheta;
        dPhi_closest = dPhi;
        dt_closest = dt;
      }

      // Skipping if the hit is outside the cut window
      if (std::fabs(dPhi) > dlCut->dPhi_max) continue;
      if (std::fabs(dTheta) > dlCut->dTheta_max) continue;

      nCompatibleHits++;
      hitAccepted[iHit2] = true;
      debug() << " Accepted 2nd hit at layer: " << layerID2 << ";  dPhi: " << dPhi << ";  dTheta: " <<  dTheta << endmsg;
    }
    // Filling diagnostic histograms
    if (m_fillHistos && dR_min < 998) {
        char hname[100];
        sprintf(hname, "h_dU_layers_%d_%d", dlCut->layer0, dlCut->layer1);
        m_histos.find(hname)->second->Fill(dU_closest);
        sprintf(hname, "h2_dU_dPhi_layers_%d_%d", dlCut->layer0, dlCut->layer1);
        m_histos.find(hname)->second->Fill(dU_closest, dPhi_closest*1e3);
        sprintf(hname, "h_dTheta_layers_%d_%d", dlCut->layer0, dlCut->layer1);
        m_histos.find(hname)->second->Fill(dTheta_closest*1e3);
        sprintf(hname, "h_dPhi_layers_%d_%d", dlCut->layer0, dlCut->layer1);
        m_histos.find(hname)->second->Fill(dPhi_closest*1e3);
        sprintf(hname, "h_dt_layers_%d_%d", dlCut->layer0, dlCut->layer1);
        m_histos.find(hname)->second->Fill(dt_closest);
      }

    // Accepting the first hit if it has at least one compatible pair
    if (nCompatibleHits > 0) {
      hitAccepted[iHit] = true;
      debug() << " Accepted 1st hit at layer: " << layerID << endmsg;
    }
  }

  // Adding accepted hits to the output collection
  size_t nHitsAccepted(0);
  for (size_t iHit = 0; iHit < nHit; iHit++) {
    if (!hitAccepted[iHit]) {
      // Filling the positions of rejected hits
      if (m_fillHistos) {
        edm4hep::TrackerHitPlane h = inputTrackerHitCollection.at( iHit );
        unsigned int layerID = bitFieldCoder.get(h.getCellID(), "layer");

        // Getting local hit position
        dd4hep::rec::Vector3D posGlobal( h.getPosition().x, h.getPosition().y, h.getPosition().z );
        dd4hep::rec::Vector2D posLocal = globalToLocal( h.getCellID(), posGlobal );

        char hname[100];
        sprintf(hname, "h2_posUV_rejected_layer_%d", layerID);
        m_histos.find(hname)->second->Fill(posLocal.u(), posLocal.v());
      }
      continue;
    }
    outCol.push_back( inputTrackerHitCollection.at( iHit ) );
    nHitsAccepted++;
  }
	info() << " " << nHitsAccepted << " hits added to collection." << endmsg;

    return outCol;
}

StatusCode FilterDoubleLayerHits::finalize(){
  return StatusCode::SUCCESS;
}