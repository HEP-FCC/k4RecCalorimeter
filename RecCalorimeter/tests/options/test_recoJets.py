#
# IMPORTS
#

# Logger
from Gaudi.Configuration import INFO
# units and physical constants
from GaudiKernel.PhysicalConstants import pi

#
# SETTINGS
#
inputfile = "recoJets_sim.root"   # input file produced with ddsim
outputfile = "recoJets_rec.root"  # output file to be produced

#
# ALGORITHMS AND SERVICES SETUP
#
TopAlg = []  # alg sequence
ExtSvc = []  # list of external services


# CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
ExtSvc += [audsvc]


# Detector geometry
# prefix all xmls with path_to_detector
# if K4GEO is empty, this should use relative path to working directory
import os
from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc")
path_to_detector = os.environ.get("K4GEO", "") + '/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/'
detectors_to_use = [
    'ALLEGRO_o1_v03.xml'
]
# prefix all xmls with path_to_detector
geoservice.detectors = [
    os.path.join(path_to_detector, _det) for _det in detectors_to_use
]
geoservice.OutputLevel = INFO
ExtSvc += [geoservice]

# retrieve subdetector IDs
import xml.etree.ElementTree as ET
tree = ET.parse(path_to_detector + 'DectDimensions.xml')
root = tree.getroot()
IDs = {}
for constant in root.find('define').findall('constant'):
    if (constant.get('name') == 'DetID_VXD_Barrel' or
        constant.get('name') == 'DetID_VXD_Disks' or
        constant.get('name') == 'DetID_DCH' or
        constant.get('name') == 'DetID_SiWr_Barrel' or
        constant.get('name') == 'DetID_SiWr_Disks' or
        constant.get('name') == 'DetID_ECAL_Barrel' or
        constant.get('name') == 'DetID_ECAL_Endcap' or
        constant.get('name') == 'DetID_HCAL_Barrel' or
        constant.get('name') == 'DetID_HCAL_Endcap' or
        constant.get('name') == 'DetID_Muon_Barrel'):
        IDs[constant.get("name")[6:]] = int(constant.get('value'))
    if (constant.get('name') == 'DetID_Muon_Endcap_1'):
        IDs[constant.get("name")[6:-2]] = int(constant.get('value'))

# Input/Output handling
from k4FWCore import IOSvc
from Configurables import EventDataSvc
io_svc = IOSvc("IOSvc")
io_svc.Input = inputfile
io_svc.Output = outputfile
ExtSvc += [EventDataSvc("EventDataSvc")]


# Calorimeter digitisation
# Detector readouts
# ECAL
ecalBarrelReadoutName = "ECalBarrelModuleThetaMerged"
ecalEndcapReadoutName = "ECalEndcapTurbine"

# EM scale calibration (sampling fraction)
from Configurables import CalibrateInLayersTool
calibEcalBarrel = CalibrateInLayersTool("CalibrateECalBarrel",
                                        samplingFraction=[0.3864252122990472] * 1 + [0.13597644835735828] * 1 + [0.14520427829645913] * 1 + [0.1510076084632846] * 1 + [0.1552347580991012] * 1 + [0.159694330729184] * 1 + [0.1632954482794191] * 1 + [0.16720711037339814] * 1 + [0.17047749048884808] * 1 + [0.17461698117974286] * 1 + [0.1798984163980135] * 1 + [0.17920355117405806] * 1,
                                        readoutName=ecalBarrelReadoutName,
                                        layerFieldName="layer")


calibEcalEndcap = CalibrateInLayersTool("CalibrateECalEndcap",
                                        samplingFraction = [0.16419] * 1 + [0.192898] * 1 + [0.18783] * 1 + [0.193203] * 1 + [0.193928] * 1 + [0.192286] * 1 + [0.199959] * 1 + [0.200153] * 1 + [0.212635] * 1 + [0.180345] * 1 + [0.18488] * 1 + [0.194762] * 1 + [0.197775] * 1 + [0.200504] * 1 + [0.205555] * 1 + [0.203601] * 1 + [0.210877] * 1 + [0.208376] * 1 + [0.216345] * 1 + [0.201452] * 1 + [0.202134] * 1 + [0.207566] * 1 + [0.208152] * 1 + [0.209889] * 1 + [0.211743] * 1 + [0.213188] * 1 + [0.215864] * 1 + [0.22972] * 1 + [0.192515] * 1 + [0.0103233] * 1,
                                        readoutName=ecalEndcapReadoutName,
                                        layerFieldName="layer")

# Step 1: merge hits into cells according to initial segmentation
ecalBarrelCellsName = "ECalBarrelCells"
from Configurables import CreateCaloCells
createEcalBarrelCells = CreateCaloCells("CreateECalBarrelCells",
                                        doCellCalibration=True,
                                        calibTool=calibEcalBarrel,
                                        crosstalksTool=None,
                                        addCrosstalk=False,
                                        addCellNoise=False,
                                        filterCellNoise=False,
                                        addPosition=True,
                                        OutputLevel=INFO,
                                        hits=ecalBarrelReadoutName,
                                        cells=ecalBarrelCellsName)
TopAlg += [createEcalBarrelCells]


# Add to Ecal barrel cells the position information
# (good for physics, all coordinates set properly)
from Configurables import CellPositionsECalBarrelModuleThetaSegTool
cellPositionEcalBarrelTool = CellPositionsECalBarrelModuleThetaSegTool(
    "CellPositionsECalBarrel",
    readoutName=ecalBarrelReadoutName,
    OutputLevel=INFO
)
ecalBarrelPositionedCellsName = "ECalBarrelPositionedCells"
from Configurables import CreateCaloCellPositionsFCCee
createEcalBarrelPositionedCells = CreateCaloCellPositionsFCCee(
    "CreateECalBarrelPositionedCells",
    OutputLevel=INFO
)
TopAlg += [createEcalBarrelPositionedCells]
createEcalBarrelPositionedCells.positionsTool = cellPositionEcalBarrelTool
createEcalBarrelPositionedCells.hits.Path = ecalBarrelCellsName
createEcalBarrelPositionedCells.positionedHits.Path = ecalBarrelPositionedCellsName

# Create cells in ECal endcap
ecalEndcapCellsName = "ECalEndcapCells"
createEcalEndcapCells = CreateCaloCells("CreateEcalEndcapCaloCells",
                                        doCellCalibration=True,
                                        calibTool=calibEcalEndcap,
                                        addCellNoise=False,
                                        filterCellNoise=False,
                                        OutputLevel=INFO,
                                        hits=ecalEndcapReadoutName,
                                        cells=ecalEndcapCellsName)
TopAlg += [createEcalEndcapCells]

# Add to Ecal endcap cells the position information
# (good for physics, all coordinates set properly)
from Configurables import CellPositionsECalEndcapTurbineSegTool
cellPositionEcalEndcapTool = CellPositionsECalEndcapTurbineSegTool(
    "CellPositionsECalEndcap",
    readoutName=ecalEndcapReadoutName,
     OutputLevel=INFO
)
ecalEndcapPositionedCellsName = "ECalEndcapPositionedCells"
createEcalEndcapPositionedCells = CreateCaloCellPositionsFCCee(
    "CreateECalEndcapPositionedCells",
    OutputLevel=INFO
)
TopAlg += [createEcalEndcapPositionedCells]
createEcalEndcapPositionedCells.positionsTool = cellPositionEcalEndcapTool
createEcalEndcapPositionedCells.hits.Path = ecalEndcapCellsName
createEcalEndcapPositionedCells.positionedHits.Path = ecalEndcapPositionedCellsName

#Empty cells for parts of calorimeter not implemented yet
from Configurables import CreateEmptyCaloCellsCollection
createemptycells = CreateEmptyCaloCellsCollection("CreateEmptyCaloCells")
createemptycells.cells.Path = "emptyCaloCells"
TopAlg += [createemptycells]

# Produce sliding window clusters
from Configurables import CaloTowerToolFCCee
towers = CaloTowerToolFCCee("towers",
                            #deltaThetaTower=4 * 0.009817477 / 4, deltaPhiTower=2 * 2 * pi / 1536.,
                            deltaThetaTower=0.022180,
                            deltaPhiTower=2 * pi / 256.,
                            cells=[ecalBarrelPositionedCellsName,
                                   ecalEndcapPositionedCellsName],
                            calorimeterIDs=[
                                IDs["ECAL_Barrel"],
                                IDs["ECAL_Endcap"],
                                ],
                            OutputLevel=INFO)

# Cluster variables (not optimized)
windT = 18
windP = 34
posT = 10
posP = 22
dupT = 14
dupP = 26
finT = 18
finP = 34
# Minimal energy to create a cluster in GeV (FCC-ee detectors have to reconstruct low energy particles)
threshold = 0.5

from Configurables import CreateCaloClustersSlidingWindowFCCee
createClusters = CreateCaloClustersSlidingWindowFCCee("CreateClusters",
                                                      towerTool=towers,
                                                      nThetaWindow=windT, nPhiWindow=windP,
                                                      nThetaPosition=posT, nPhiPosition=posP,
                                                      nThetaDuplicates=dupT, nPhiDuplicates=dupP,
                                                      nThetaFinal=finT, nPhiFinal=finP,
                                                      energyThreshold=threshold,
                                                      energySharingCorrection=False,
                                                      createClusterCellCollection=True,
                                                      OutputLevel=INFO
                                                      )
createClusters.clusters.Path = "CaloClusters"
createClusters.clusterCells.Path = "CaloClusterCells"
TopAlg += [createClusters]


from Configurables import CreateCaloJet
createJets = CreateCaloJet(
    "createJets",
    InputClusterCollection=[createClusters.clusters.Path],
    OutputJetCollection=["Jets"],
    # JetAlg="antikt",
    # JetRadius=0.4,
    # OutputLevel=INFO,
    MinPt=5
)
TopAlg += [createJets]

# Configure output
io_svc.outputCommands = ["keep *",
                         "drop emptyCaloCells", 
                         "drop ECalBarrelHits", "drop HCal*", "drop ECalBarrelCellsStep*", "drop ECalBarrelPositionedHits", "drop ECalEndcapPositionedHits", "drop emptyCaloCells", "drop CaloClusterCells"]


# configure the application
print(TopAlg)
print(ExtSvc)
from k4FWCore import ApplicationMgr
applicationMgr = ApplicationMgr(
    TopAlg=TopAlg,
    EvtSel='NONE',
    EvtMax=-1,
    ExtSvc=ExtSvc,
    StopOnSignal=True,
)

for algo in applicationMgr.TopAlg:
    algo.AuditExecute = True
    # for debug
    # algo.OutputLevel = DEBUG
