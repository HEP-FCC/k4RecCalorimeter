#ifndef FilterDoubleLayerHits_h
#define FilterDoubleLayerHits_h 1

#include <k4FWCore/Transformer.h>
#include "k4Interface/IGeoSvc.h"
#include <GaudiKernel/ITHistSvc.h>

#include "DDRec/SurfaceManager.h"
#include <edm4hep/TrackerHitPlaneCollection.h>
#include <edm4hep/TrackerHitPlane.h>

#include <string>

#include <TH1.h>



/** Utility processor that removes tracker hits in double layers if they don't have
 *  a corresponding close-by hit in the other sublayer.
 *  Pairs of considered layers are configurable and extracted from the cellID word.
 *  Works for all four lcio hit classes.
 *
 *  @parameter InputCollection name of the hit collection with (Sim)TrackerHits/(Sim)CalorimeterHits
 *  @parameter OutputCollections ( ColName  StartLayer EndLayer )
 *
 * @author N. Bartosik, INFN Torino, S. Ferraro
 * @date  17 June 2020
 * @version $Id: $
 */

struct FilterDoubleLayerHits : public k4FWCore::Transformer<edm4hep::TrackerHitPlaneCollection(
  const edm4hep::TrackerHitPlaneCollection &)> {

protected:

  static const size_t NHITS_MAX = 10000000;

  struct SensorPosition{
    unsigned int layer;
    unsigned int side;
    unsigned int ladder;
    unsigned int module;

    bool operator<(const SensorPosition& rhs) const {
      return std::tie(layer, side, ladder, module) < std::tie(rhs.layer, rhs.side, rhs.ladder, rhs.module);
    }
  };

  /// Double layer cut struct
  struct DoubleLayerCut{
    unsigned int layer0 ;
    unsigned int layer1 ;
    double dPhi_max ;
    double dTheta_max ;
  };


 public:
  FilterDoubleLayerHits(const std::string& name, ISvcLocator* svcLoc) ;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual StatusCode initialize() ;

  /** Called for every event - the working horse.
   */
  edm4hep::TrackerHitPlaneCollection operator()(const edm4hep::TrackerHitPlaneCollection& inputTrackerHitCollection) const override;

  /** Called after data processing for clean up.
   */
  virtual StatusCode finalize();


 protected:
  dd4hep::rec::Vector2D globalToLocal(long int cellID, const dd4hep::rec::Vector3D& posGlobal, dd4hep::rec::ISurface** surf) const;

  Gaudi::Property<std::string> m_subDetName{this, "SubDetectorName", "Vertex", "Name of sub detector"};
  Gaudi::Property<bool> m_fillHistos{this, "FillHistograms", false, "Whether to fill diagnostic histograms"};
  Gaudi::Property<double> m_dtMax{this, "DeltaTimeMax", -1.0, "Maximum time difference between hits in a doublet [ns]"};
  Gaudi::Property<std::string> m_encodingStringVariable{this, "EncodingStringParameterName", "GlobalTrackerReadoutID", "The name of the DD4hep constant that contains the Encoding string for tracking detectors"};
  Gaudi::Property<std::vector<std::string>> m_dlCutConfigs{this, "DoubleLayerCuts" , {"0", "1", "0.5", "0.05"}, "Layer IDs and angular cuts [mrad] to be applied: <layer 0> <layer 1> <dPhi> <dTheta>"};

  ////Double layer cuts configuration
  std::vector<DoubleLayerCut> m_dlCuts {};
  ////Surface map for getting local hit positions at sensor surface
  const dd4hep::rec::SurfaceMap* m_map {nullptr};
  ////Monitoring histograms
  std::map<std::string, TH1*> m_histos {};
  // GeoSvc
  SmartIF<IGeoSvc>   m_geoSvc;
  SmartIF<ITHistSvc> m_histSvc;
} ;

#endif