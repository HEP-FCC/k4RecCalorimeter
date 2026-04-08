# Number of events
Nevts = 3
# input/output files
inputfile = "output_fullCalo_SimDigi_e_50GeV_"+str(Nevts)+"events.root"
outputfile = "digi_cellPositions_50GeVelectrons.root"

#
# ALGORITHMS AND SERVICES SETUP
#
TopAlg = []  # alg sequence
ExtSvc = []  # list of external services

from Gaudi.Configuration import INFO, WARNING, DEBUG

# Event counter
from Configurables import EventCounter
eventCounter = EventCounter("EventCounter",
                            OutputLevel=INFO,
                            Frequency=1)
TopAlg += [eventCounter]

# CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
ExtSvc += [audsvc]

# Detector geometry
# prefix all xmls with path_to_detector
# if K4GEO is empty, this should use relative path to working directory
from Configurables import GeoSvc
import os
geoservice = GeoSvc("GeoSvc",
                    OutputLevel=WARNING
                    )

path_to_detector = os.environ.get("K4GEO", "") + "/FCChh/compact/FCChhBaseline/"
detectors_to_use = [
    'FCChh_DectMaster.xml'
]
geoservice.detectors = [
    os.path.join(path_to_detector, _det) for _det in detectors_to_use
]

#detectors_to_use=['file:Detector/DetFCChhBaseline1/compact/FCChh_DectEmptyMaster.xml',
#                  # Tracker disabled to save cpu time
#                  #'file:Detector/DetFCChhTrackerTkLayout/compact/Tracker.xml',
#                  'file:Detector/DetFCChhECalInclined/compact/FCChh_ECalBarrel_withCryostat.xml',
#                  'file:Detector/DetFCChhHCalTile/compact/FCChh_HCalBarrel_TileCal.xml',
#                  'file:Detector/DetFCChhHCalTile/compact/FCChh_HCalExtendedBarrel_TileCal.xml',
#                  'file:Detector/DetFCChhCalDiscs/compact/Endcaps_coneCryo.xml',
#                  'file:Detector/DetFCChhCalDiscs/compact/Forward_coneCryo.xml',
#                  ]
# geoservice.detector = detectors_to_use
ExtSvc += [geoservice]

# Input/Output handling
from k4FWCore import IOSvc
from Configurables import EventDataSvc
io_svc = IOSvc("IOSvc")
io_svc.Input = inputfile
io_svc.Output = outputfile
ExtSvc += [EventDataSvc("EventDataSvc")]

# Configure tools for calo reconstruction

# ECAL readouts
ecalBarrelReadoutName = "ECalBarrelEta"
ecalBarrelReadoutNamePhiEta = "ECalBarrelPhiEta"
ecalEndcapReadoutName = "EMECPhiEtaReco"
ecalFwdReadoutName = "EMFwdPhiEtaReco"
# HCAL readouts
hcalBarrelReadoutName = "HCalBarrelReadout"
hcalExtBarrelReadoutName = "HCalExtBarrelReadout"
hcalEndcapReadoutName = "HECPhiEtaReco"
hcalFwdReadoutName = "HFwdPhiEtaReco"
# Number of events
num_events = 3

# Configure tools for calo cell positions
from Configurables import CellPositionsECalBarrelTool, CellPositionsHCalBarrelNoSegTool, CellPositionsCaloDiscsTool, CellPositionsTailCatcherTool 
ECalBcells = CellPositionsECalBarrelTool("CellPositionsECalBarrel", 
                                         readoutName = ecalBarrelReadoutNamePhiEta, 
                                         OutputLevel = INFO)
EMECcells = CellPositionsCaloDiscsTool("CellPositionsEMEC", 
                                       readoutName = ecalEndcapReadoutName, 
                                       OutputLevel = INFO)
ECalFwdcells = CellPositionsCaloDiscsTool("CellPositionsECalFwd", 
                                          readoutName = ecalFwdReadoutName, 
                                          OutputLevel = INFO)
HCalBcells = CellPositionsHCalBarrelNoSegTool("CellPositionsHCalBarrel", 
                                              readoutName = hcalBarrelReadoutName, 
                                              OutputLevel = INFO)
HCalExtBcells = CellPositionsHCalBarrelNoSegTool("CellPositionsHCalExtBarrel", 
                                                 readoutName = hcalExtBarrelReadoutName, 
                                                 OutputLevel = INFO)
HECcells = CellPositionsCaloDiscsTool("CellPositionsHEC", 
                                      readoutName = hcalEndcapReadoutName, 
                                      OutputLevel = INFO)
HCalFwdcells = CellPositionsCaloDiscsTool("CellPositionsHCalFwd", 
                                          readoutName = hcalFwdReadoutName, 
                                          OutputLevel = INFO)

# cell positions
from Configurables import CreateCellPositions
positionsEcalBarrel = CreateCellPositions("positionsEcalBarrel", 
                                          positionsTool=ECalBcells, 
                                          hits = "ECalBarrelCells", 
                                          positionedHits = "ECalBarrelCellPositions", 
                                          OutputLevel = INFO)
positionsHcalBarrel = CreateCellPositions("positionsHcalBarrel", 
                                          positionsTool=HCalBcells, 
                                          hits = "HCalBarrelCells", 
                                          positionedHits = "HCalBarrelCellPositions", 
                                          OutputLevel = INFO)
positionsHcalExtBarrel = CreateCellPositions("positionsHcalExtBarrel", 
                                             positionsTool=HCalExtBcells, 
                                             hits = "HCalExtBarrelCells", 
                                             positionedHits = "HCalExtBarrelCellPositions", 
                                             OutputLevel = INFO)
positionsEcalEndcap = CreateCellPositions("positionsEcalEndcap", 
                                          positionsTool=EMECcells, 
                                          hits = "ECalEndcapCells", 
                                          positionedHits = "ECalEndcapCellPositions", 
                                          OutputLevel = INFO)
positionsHcalEndcap = CreateCellPositions("positionsHcalEndcap", 
                                          positionsTool=HECcells, 
                                          hits = "HCalEndcapCells", 
                                          positionedHits = "HCalEndcapCellPositions", 
                                          OutputLevel = INFO)
positionsEcalFwd = CreateCellPositions("positionsEcalFwd", 
                                       positionsTool=ECalFwdcells, 
                                       hits = "ECalFwdCells", 
                                       positionedHits = "ECalFwdCellPositions", 
                                       OutputLevel = INFO)
positionsHcalFwd = CreateCellPositions("positionsHcalFwd", 
                                       positionsTool=HCalFwdcells, 
                                       hits = "HCalFwdCells", 
                                       positionedHits = "HCalFwdCellPositions", 
                                       OutputLevel = INFO)

TopAlg += [    
    positionsEcalBarrel,
    positionsEcalEndcap,
    positionsEcalFwd, 
    positionsHcalBarrel, 
    positionsHcalExtBarrel,
    positionsHcalEndcap, 
    positionsHcalFwd,
]

io_svc.outputCommands = ["keep *","drop ECalBarrelCells","drop ECalEndcapCells","drop ECalFwdCells","drop HCalBarrelCells", "drop HCalExtBarrelCells", "drop HCalEndcapCells", "drop HCalFwdCells"]


# configure the application
print(TopAlg)
print(ExtSvc)
from k4FWCore import ApplicationMgr
applicationMgr = ApplicationMgr(
    TopAlg=TopAlg,
    EvtSel='NONE',
    EvtMax=Nevts,
    ExtSvc=ExtSvc,
    StopOnSignal=True,
)

for algo in applicationMgr.TopAlg:
    algo.AuditExecute = True
