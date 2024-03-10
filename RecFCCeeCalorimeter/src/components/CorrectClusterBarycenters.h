#ifndef RECFCCEECALORIMETER_CORRECTCLUSTERBARYCENTERS_H
#define RECFCCEECALORIMETER_CORRECTCLUSTERBARYCENTERS_H

// Key4HEP
#include "k4FWCore/DataHandle.h"
class IGeoSvc;


// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/RndmGenerators.h"

//EDM4HEP
namespace edm4hep {
    class ClusterCollection;
}

namespace dd4hep{
    namespace DDSegmentation {
        class FCCSWGridPhiTheta_k4geo;
        class MultiSegmentation;
        class BitFieldCoder;
    }
}

class CorrectClusterBarycenters : public GaudiAlgorithm{

    public :
    CorrectClusterBarycenters(const std::string& name, ISvcLocator* svcLoc);
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();

    private:

    unsigned int phiNeighbour(int aIPhi, int aMaxPhi) const;

    /// Handle for clusters (input collection)
    DataHandle<edm4hep::ClusterCollection> m_inClusters{"clusters", Gaudi::DataHandle::Reader, this};
    /// Handle for corrected clusters (output collection)
    DataHandle<edm4hep::ClusterCollection> m_correctedClusters{"correctedClusters", Gaudi::DataHandle::Writer, this};
    /// Pointer to the geometry service
    ServiceHandle<IGeoSvc> m_geoSvc;
    /// Number of layers
    Gaudi::Property<uint> m_numLayers{this, "numLayers", 12, "Number of layers in the theta direction"};
    /// Weights for each detector layer for theta position log-weighting
    Gaudi::Property<std::vector<double>> m_thetaRecalcLayerWeights{
        this,
        "ThetaRecalcWeights",
            {0, 3.0, 3.0, 3.0, 4.25, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0},
        "Weights for each detector layer for theta position log-weighting"};
    Gaudi::Property<std::vector<double>> m_thetaLayerResolutionSampling{
        this,
        "ThetaLayerResolutionSampling",
            {0, 1.405, 2.051, 3.195, 5.791, 10.28, 17.54, 30.5},
        "sampling term of Theta resolution"};
    Gaudi::Property<std::vector<double>> m_thetaLayerResolutionConst{
        this,
        "ThetaLayerResolutionConst",
            {2.623, 0.1083, 0.1152, 0, 0, 0, 0, 0},
        "const term of Theta resolution"};
    /// General decoder to encode the calorimeter sub-system to determine which positions tool to use
    dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
    /// IDs of the detectors
    Gaudi::Property<uint> m_systemId{this, "systemId", 4, "IDs of systems"};
    /// Names of the detector readout, corresponding to system IDs in m_systemId
    Gaudi::Property<std::string> m_readoutName{
        this, "readoutName", "ECalBarrelModuleThetaMerged", "Names of the detector readout, corresponding to systemId"};
    /// map of system Id to segmentation, created based on m_readoutName and m_systemId
    dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo* m_segmentationPhiTheta;
    // dd4hep::DDSegmentation::MultiSegmentation* m_segmentationMulti;
    /// Name of the layer/cell field
    Gaudi::Property<std::string> m_layerFieldName{this, "layerFieldName", "layer", "Identifier of layers"};
    /// Id of the first layer
    Gaudi::Property<uint> m_firstLayerId{this, "firstLayerId", 0, "ID of first layer"};
    /// Size of the window in phi for the final cluster building, optimised for each layer  (in units of cell size)
    /// If empty use same size for each layer, as in *nPhiFinal*
    Gaudi::Property<std::vector<int>> m_nPhiFinal{this, "nPhiOptimFinal", {}};
    // Recalculate to half size N (window size = 2*N+1)
    std::vector<int> m_halfPhiFin;
    /// Size of the window in Theta for the final cluster building, optimised for each layer  (in units of cell size)
    /// If empty use same size for each layer, as in *nThetaFinal*
    Gaudi::Property<std::vector<int>> m_nThetaFinal{this, "nThetaOptimFinal", {}};
    // Recalculate to half size N (window size = 2*N+1)
    std::vector<int> m_halfThetaFin;
};

#endif  /* RECFCCEECALORIMETER_CORRECTCLUSTERBARYCENTERS_H */