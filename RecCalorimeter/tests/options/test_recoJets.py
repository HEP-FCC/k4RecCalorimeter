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

ecalBarrelLayers = 11
ecalBarrelSamplingFraction = [0.3800493723322256] * 1 + [0.13494147915064658] * 1 + [0.142866851721152] * 1 + [0.14839315921940666] * 1 + [0.15298362570665006] * 1 + [0.15709704561942747] * 1 + [0.16063717490147533] * 1 + [0.1641723795419055] * 1 + [0.16845490287689746] * 1 + [0.17111520115997653] * 1 + [0.1730605163148862] * 1
ecalEndcapLayers = 98
ecalEndcapSamplingFraction = [0.0897818] * 1+ [0.221318] * 1+ [0.0820002] * 1+ [0.994281] * 1+ [0.0414437] * 1+ [0.1148] * 1+ [0.178831] * 1+ [0.142449] * 1+ [0.181206] * 1+ [0.342843] * 1+ [0.137479] * 1+ [0.176479] * 1+ [0.153273] * 1+ [0.195836] * 1+ [0.0780405] * 1+ [0.150202] * 1+ [0.17846] * 1+ [0.164886] * 1+ [0.175758] * 1+ [0.10836] * 1+ [0.160243] * 1+ [0.183373] * 1+ [0.171818] * 1+ [0.194848] * 1+ [0.111899] * 1+ [0.170704] * 1+ [0.188455] * 1+ [0.178164] * 1+ [0.209113] * 1+ [0.105241] * 1+ [0.180637] * 1+ [0.192206] * 1+ [0.186096] * 1+ [0.211962] * 1+ [0.112019] * 1+ [0.180344] * 1+ [0.195684] * 1+ [0.190778] * 1+ [0.218259] * 1+ [0.118516] * 1+ [0.207786] * 1+ [0.204474] * 1+ [0.207048] * 1+ [0.225913] * 1+ [0.111325] * 1+ [0.147875] * 1+ [0.195625] * 1+ [0.173326] * 1+ [0.175449] * 1+ [0.104087] * 1+ [0.153645] * 1+ [0.161263] * 1+ [0.165499] * 1+ [0.171758] * 1+ [0.175789] * 1+ [0.180657] * 1+ [0.184563] * 1+ [0.187876] * 1+ [0.191762] * 1+ [0.19426] * 1+ [0.197959] * 1+ [0.199021] * 1+ [0.204428] * 1+ [0.195709] * 1+ [0.151751] * 1+ [0.171477] * 1+ [0.165509] * 1+ [0.172565] * 1+ [0.172961] * 1+ [0.175534] * 1+ [0.177989] * 1+ [0.18026] * 1+ [0.181898] * 1+ [0.183912] * 1+ [0.185654] * 1+ [0.187515] * 1+ [0.190408] * 1+ [0.188794] * 1+ [0.193699] * 1+ [0.192287] * 1+ [0.19755] * 1+ [0.190943] * 1+ [0.218553] * 1+ [0.161085] * 1+ [0.373086] * 1+ [0.122495] * 1+ [0.21103] * 1+ [1] * 1+ [0.138686] * 1+ [0.0545171] * 1+ [1] * 1+ [1] * 1+ [0.227945] * 1+ [0.0122872] * 1+ [0.00437334] * 1+ [0.00363533] * 1+ [1] * 1+ [1] * 1


# EM scale calibration (sampling fraction)
from Configurables import CalibrateInLayersTool
calibEcalBarrel = CalibrateInLayersTool("CalibrateECalBarrel",
                                        samplingFraction=ecalBarrelSamplingFraction,
                                        readoutName=ecalBarrelReadoutName,
                                        layerFieldName="layer")


calibEcalEndcap = CalibrateInLayersTool("CalibrateECalEndcap",
                                        samplingFraction = ecalEndcapSamplingFraction,
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
    MinPt=2
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
