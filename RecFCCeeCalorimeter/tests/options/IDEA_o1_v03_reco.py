from Gaudi.Configuration import *
from Configurables import EventDataSvc
from k4FWCore import ApplicationMgr, IOSvc

iosvc = IOSvc()
iosvc.Input = "testIDEA_o1_v03_digi.root"

dataservice = EventDataSvc("EventDataSvc")

# detector geometry
# if K4GEO is empty, this should use relative path to working directory
from Configurables import GeoSvc
import os

geoservice = GeoSvc("GeoSvc")
path_to_detector = os.environ.get("K4GEO", "")
detectors_to_use = ["FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml"]

geoservice.detectors = [
    os.path.join(path_to_detector, _det) for _det in detectors_to_use
]

# retrieve subdetector IDs
import xml.etree.ElementTree as ET

tree = ET.parse(
    path_to_detector + "/FCCee/IDEA/compact/IDEA_o1_v03/DectDimensions_IDEA_o1_v03.xml"
)
root = tree.getroot()
IDs = {}
for constant in root.find("define").findall("constant"):
    if constant.get("name") == "DetID_FiberDRCalo":
        IDs[constant.get("name")[6:]] = int(constant.get("value"))
# debug
print("Subdetector IDs:")
print(IDs)

geoservice.OutputLevel = INFO

iosvc.CollectionNames = [
    "DRcaloSiPMreadout_scint",
    "DRcaloSiPMreadout_scintContributions",
    "DRcaloSiPMreadoutSimHit",
    "DRcaloSiPMreadoutDigiHit_scint",
    "DRcaloSiPMreadoutDigiHit",
]

# use const noise tool for the moment
from Configurables import ConstNoiseTool

constNoiseTool = ConstNoiseTool(
    "ConstNoiseTool",
    detectors=["FiberDRCalo"],
    systemEncoding="system:5",
    detectorsNoiseRMS=[0.001],  # ad-hoc small value
    detectorsNoiseOffset=[0.0],
    OutputLevel=INFO,
)

from Configurables import CaloTopoClusterFCCee

caloIDs = [IDs["FiberDRCalo"]]

topoClusterAll = CaloTopoClusterFCCee(
    "topoClusterAll",
    cells=["DRcaloSiPMreadoutDigiHit", "DRcaloSiPMreadoutDigiHit_scint"],
    clusters="TopoClusterAll",
    clusterCells="TopoClusterAllCells",
    useNeighborMap=False,
    readoutName="DRcaloSiPMreadout",
    neigboursTool=None,
    noiseTool=constNoiseTool,
    systemEncoding="system:5",
    seedSigma=4,
    neighbourSigma=2,
    lastNeighbourSigma=0,
    calorimeterIDs=caloIDs,
    createClusterCellCollection=True,
    OutputLevel=INFO,
)

iosvc.Output = "testIDEA_o1_v03_reco.root"
iosvc.outputCommands = ["keep *"]

ApplicationMgr(
    TopAlg=[topoClusterAll],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[dataservice, geoservice],
)
