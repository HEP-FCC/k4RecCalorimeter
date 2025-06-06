#ifndef RECCALORIMETER_CREATECALOCLUSTERS_H
#define RECCALORIMETER_CREATECALOCLUSTERS_H

// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/ITHistSvc.h"
#include "GaudiKernel/ToolHandle.h"

// Key4HEP
#include "k4FWCore/DataHandle.h"
#include "k4Interface/ICalorimeterTool.h"
#include "k4Interface/ICellPositionsTool.h"
#include "k4Interface/INoiseConstTool.h"

// DD4hep
#include "DDSegmentation/Segmentation.h"

// EDM4HEP
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/Cluster.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/MCParticleCollection.h"

class IGeoSvc;
namespace DD4hep {
namespace DDSegmentation {
  class Segmentation;
}
} // namespace DD4hep

class TH2F;
class TH1F;

/** @class CreateCaloClusters
 *
 * Applies hadronic calibration to cluster.
 * By default only the cluster that consist of E and HCal cells or only HCal cells get calibrated to hadronic scale.
 * By specifying m_doCryoCorrection=True the benchmark method is used, based on the ATLAS calibration scheme for
 * LAr+Tile testbeams. In this case all clusters get calibrated.
 *
 *  Tools called:
 *    - ICellPositionsTool
 *
 *  @author Coralie Neubueser
 *  @date   2019-03
 *
 */

class CreateCaloClusters : public Gaudi::Algorithm {

public:
  CreateCaloClusters(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize();

  StatusCode execute(const EventContext&) const;

  StatusCode finalize();

private:
  /// Pointer to the interface of histogram service
  SmartIF<ITHistSvc> m_histSvc;
  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;

  /// Handle for calo clusters (input collection)
  mutable k4FWCore::DataHandle<edm4hep::ClusterCollection> m_clusters{"calo/clusters", Gaudi::DataHandle::Reader, this};
  /// Handle for calo clusters (input collection)
  mutable k4FWCore::DataHandle<edm4hep::MCParticleCollection> m_genParticles{"calo/genParticles",
                                                                             Gaudi::DataHandle::Reader, this};
  /// Handle for calo clusters (output collection)
  mutable k4FWCore::DataHandle<edm4hep::ClusterCollection> m_newClusters{"calo/calibClusters",
                                                                         Gaudi::DataHandle::Writer, this};
  // Handle for calo cells (output collection)
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_newCells{"calo/calibClusterCells",
                                                                             Gaudi::DataHandle::Writer, this};

  /// Handle for tool to get positions in ECal Barrel
  ToolHandle<ICellPositionsTool> m_cellPositionsECalTool{"CellPositionsECalBarrelTool", this};
  /// Handle for tool to get positions in HCal Barrel and Ext Barrel, no Segmentation
  ToolHandle<ICellPositionsTool> m_cellPositionsHCalTool{"CellPositionsHCalBarrelTool", this};
  /// Handle for tool to get positions in HCal Barrel and Ext Barrel, no Segmentation
  ToolHandle<ICellPositionsTool> m_cellPositionsHCalNoSegTool{"CellPositionsHCalBarrelNoSegTool", this};

  mutable TH1F* m_energyScale;
  mutable TH1F* m_benchmark;
  mutable TH1F* m_fractionEMcluster;
  mutable TH2F* m_energyScaleVsClusterEnergy;
  mutable TH1F* m_totEnergy;
  mutable TH1F* m_totCalibEnergy;
  mutable TH1F* m_totBenchmarkEnergy;
  mutable TH1F* m_clusterEnergy;
  mutable TH1F* m_sharedClusterEnergy;
  mutable TH1F* m_clusterEnergyCalibrated;
  mutable TH1F* m_clusterEnergyBenchmark;
  mutable TH1F* m_nCluster;
  mutable TH1F* m_nCluster_1GeV;
  mutable TH1F* m_nCluster_halfTrueEnergy;
  mutable TH1F* m_energyCalibCluster_1GeV;
  mutable TH1F* m_energyCalibCluster_halfTrueEnergy;

  /// bool if calibration is applied
  bool m_doCalibration = true;

  /// e/h of ECal
  double m_ehECal;
  /// e/h of HCal
  double m_ehHCal;
  /// bool if energy loss needs correction is applied
  bool m_doCryoCorrection = false;

  dd4hep::DDSegmentation::BitFieldCoder* m_decoder = new dd4hep::DDSegmentation::BitFieldCoder("system:4");
  dd4hep::DDSegmentation::BitFieldCoder* m_decoderECal;
  dd4hep::DDSegmentation::BitFieldCoder* m_decoderHCal;

  /// System id by default Barrel, EC(6,7), Fwd(10,11)
  /// with .. x^b2 .. : -5.53466e-07,4.73147e-11,-1.73903e-05,1515.84,0.823583,-4.87235,150252,9.8425e+09,0.326512
  /// from chi2 minimisation, no C: 0.975799,-2.54738e-06,0.822663,-0.140975,-2.18657e-05,-0.0193682
  Gaudi::Property<float> m_a1{this, "a1", 0.957, "scaling of ECal energy"};                  // no Bfield: 0.9867
  Gaudi::Property<float> m_a2{this, "a2", 3.772e-08, "scaling of ECal energy"};              // no Bfield: 0.9867
  Gaudi::Property<float> m_a3{this, "a3", 1.028, "scaling of ECal energy"};                  // no Bfield: 0.9867
  Gaudi::Property<float> m_b1{this, "b1", 0.243, "scaling of energy loss in cryostat"};      // no Bfield: 0.432
  Gaudi::Property<float> m_b2{this, "b2", 0.097, "scaling of energy loss in cryostat"};      // no Bfield: 0.432
  Gaudi::Property<float> m_b3{this, "b3", -0.217, "scaling of energy loss in cryostat"};     // no Bfield: 0.432
  Gaudi::Property<float> m_c1{this, "c1", -1.546e-10, "scaling of energy loss in cryostat"}; // no Bfield: -5.567E-6
  Gaudi::Property<float> m_c2{this, "c2", -4.519e+06, "scaling of energy loss in cryostat"}; // no Bfield: -5.567E-6
  Gaudi::Property<float> m_c3{this, "c3", 7.026, "scaling of energy loss in cryostat"};      // no Bfield: -5.567E-6
  /// no segmentation used in HCal
  Gaudi::Property<bool> m_noSegmentationHCal{this, "noSegmentationHCal", true,
                                             "HCal readout w/o eta-phi segementation?"};
  Gaudi::Property<int> m_lastECalLayer{this, "lastECalLayer", 7, "Layer id of last ECal layer"};
  Gaudi::Property<int> m_firstHCalLayer{this, "firstHCalLayer", 0, "Layer id of first HCal layer"};

  Gaudi::Property<uint> m_systemIdECal{this, "systemECal", 5, "System id of ECal"};
  Gaudi::Property<uint> m_systemIdHCal{this, "systemHCal", 8, "System id of HCal"};
  Gaudi::Property<std::string> m_readoutECal{this, "readoutECal", "Readout of ECal"};
  Gaudi::Property<std::string> m_readoutHCal{this, "readoutHCal", "Readout of HCal"};

  Gaudi::Property<double> m_fractionECal{this, "fractionECal", 0.7,
                                         "Fraction of clsuter energy in ECal to be flagged as EM"};
};

#endif /* RECCALORIMETER_CREATECALOCLUSTERS_H */
