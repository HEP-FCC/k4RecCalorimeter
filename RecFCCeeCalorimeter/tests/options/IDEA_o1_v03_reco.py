from Gaudi.Configuration import *
from Configurables import ApplicationMgr

from Configurables import k4DataSvc
dataservice = k4DataSvc("EventDataSvc", input="testIDEA_o1_v03_digi.root")

# detector geometry
# if K4GEO is empty, this should use relative path to working directory
from Configurables import GeoSvc
import os
geoservice = GeoSvc("GeoSvc")
path_to_detector = "/afs/cern.ch/work/s/sako/private/kfc-dream/k4geo/" # os.environ.get("K4GEO", "")
detectors_to_use = [
    'FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml'
]

geoservice.detectors = [
    os.path.join(path_to_detector, _det) for _det in detectors_to_use
]

geoservice.OutputLevel = INFO

from Configurables import PodioInput
podioinput = PodioInput("PodioInput",
    collections = [
        "DRcaloSiPMreadout_scint",
        "DRcaloSiPMreadout_scintContributions",
        "DRcaloSiPMreadoutSimHit",
        "DRcaloSiPMreadoutDigiHit_scint",
        "DRcaloSiPMreadoutDigiHit"
    ],
    OutputLevel = DEBUG
)

# use const noise tool for the moment
from Configurables import ConstNoiseTool
constNoiseTool = ConstNoiseTool("ConstNoiseTool",
    detectors = ["FiberDRCalo"],
    detectorsNoiseRMS = [0.001], # ad-hoc small value
    detectorsNoiseOffset = [0.],
    OutputLevel = INFO
)

from Configurables import CaloTopoClusterFCCee
topoClusterCheren = CaloTopoClusterFCCee("topoClusterCheren",
    cells = ["DRcaloSiPMreadoutDigiHit"],
    clusters = "TopoClusterCheren",
    clusterCells = "TopoClusterCherenCells",
    useNeighborMap = False,
    readoutName = "DRcaloSiPMreadout",
    neigboursTool = None,
    noiseTool = constNoiseTool,
    seedSigma = 4,
    neighbourSigma = 2,
    lastNeighbourSigma = 0,
    OutputLevel = INFO
)

topoClusterScint = CaloTopoClusterFCCee("topoClusterScint",
    cells = ["DRcaloSiPMreadoutDigiHit_scint"],
    clusters = "TopoClusterScint",
    clusterCells = "TopoClusterScintCells",
    useNeighborMap = False,
    readoutName = "DRcaloSiPMreadout",
    neigboursTool = None,
    noiseTool = constNoiseTool,
    seedSigma = 4,
    neighbourSigma = 2,
    lastNeighbourSigma = 0,
    OutputLevel = INFO
)

from Configurables import PodioOutput
podiooutput = PodioOutput("PodioOutput", filename = "testIDEA_o1_v03_reco.root", OutputLevel = DEBUG)
podiooutput.outputCommands = ["keep *"]

ApplicationMgr(
    TopAlg = [
        podioinput,
        topoClusterCheren,
        topoClusterScint,
        podiooutput
    ],
    EvtSel = 'NONE',
    EvtMax = 10,
    ExtSvc = [dataservice,geoservice]
)
