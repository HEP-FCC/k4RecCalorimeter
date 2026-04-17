# simulations setup
Nevts = 3
inputfile = "FCChh_sim_e_50GeV.root"
outputfile = "output_fullCalo_SimAndDigi_e_50GeV_"+str(Nevts)+"events.root"

# ECAL readouts
ecalBarrelReadoutName = "ECalBarrelEta"
ecalBarrelReadoutNamePhiEta = "ECalBarrelPhiEta"
ecalEndcapReadoutName = "EMECPhiEta"
ecalFwdReadoutName = "EMFwdPhiEta"
# HCAL readouts
hcalReadoutName = "HCalBarrelReadout"
extHcalReadoutName = "HCalExtBarrelReadout"
hcalEndcapReadoutName = "HECPhiEta"
hcalFwdReadoutName = "HFwdPhiEta"
# layers to be merged in endcaps & forward calo
ecalEndcapNumberOfLayersToMerge = [26]*5+[27]
ecalFwdNumberOfLayersToMerge = [7]*5+[8]
hcalEndcapNumberOfLayersToMerge = [13]+[14]*5
hcalFwdNumberOfLayersToMerge = [8]+[9]*5
identifierName = "layer"
volumeName = "layer"

#
# ALGORITHMS AND SERVICES SETUP
#
TopAlg = []  # alg sequence
ExtSvc = []  # list of external services

from Gaudi.Configuration import INFO, DEBUG

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
                    OutputLevel=INFO
                    # OutputLevel=DEBUG  # set to DEBUG to print dd4hep::DEBUG messages in k4geo C++ drivers
                    )

path_to_detector = os.environ.get("K4GEO", "") + "/FCChh/compact/FCChhBaseline/"
detectors_to_use = [
    'FCChh_DectMaster.xml'
]
geoservice.detectors = [
    os.path.join(path_to_detector, _det) for _det in detectors_to_use
]
ExtSvc += [geoservice]

# Input/Output handling
from k4FWCore import IOSvc
from Configurables import EventDataSvc
io_svc = IOSvc("IOSvc")
io_svc.Input = inputfile
io_svc.Output = outputfile
ExtSvc += [EventDataSvc("EventDataSvc")]

# Configure tools for calo reconstruction
# EM scale calibration
from Configurables import CalibrateCaloHitsTool
calibHcells = CalibrateCaloHitsTool("CalibrateHCal", invSamplingFraction="41.66")

from Configurables import CalibrateInLayersTool
calibEcalBarrel = CalibrateInLayersTool("CalibrateECalBarrel",
                                   # sampling fraction obtained using SamplingFractionInLayers from DetStudies package
                                   samplingFraction = [0.36571381189697705] * 1 + [0.09779064189677973] * 1 + [0.12564152224404024] * 1 + [0.14350599973146283] * 1 + [0.1557126972314961] * 1 + [0.16444759076233928] * 1 + [0.17097165096847836] * 1 + [0.17684775359805122] * 1 + [0.18181154293837265] * 1 + [0.18544247938196395] * 1 + [0.18922747431624687] * 1 + [0.21187001375505543] * 1,
                                   readoutName = ecalBarrelReadoutName,
                                   layerFieldName = "layer")

calibEcalEndcap = CalibrateCaloHitsTool("CalibrateECalEndcap", invSamplingFraction="13.89")
calibEcalFwd = CalibrateCaloHitsTool("CalibrateECalFwd", invSamplingFraction="303.03")
calibHcalEndcap = CalibrateCaloHitsTool("CalibrateHCalEndcap", invSamplingFraction="33.62")
calibHcalFwd = CalibrateCaloHitsTool("CalibrateHCalFwd", invSamplingFraction="1207.7")

# Create cells in ECal barrel
# 1. step - merge hits into cells with default Eta segmentation
# 2. step - rewrite the cellId using the Phi-Eta segmentation
from Configurables import CreateCaloCells
createEcalBarrelCellsStep1 = CreateCaloCells("CreateECalBarrelCellsStep1",
                                             doCellCalibration=True,
                                             calibTool = calibEcalBarrel,
                                             addCellNoise=False,
                                             filterCellNoise=False,
                                             OutputLevel=INFO,
                                             hits=ecalBarrelReadoutName,
                                             cells="ECalBarrelCellsStep1",
                                             addPosition=True)
TopAlg += [createEcalBarrelCellsStep1]

# Ecal barrel cell positions
# does not work - should not needed since only x-y position is needed, which is correctly
# calculated by CreateCaloCells with addPosition=True
# from Configurables import CreateVolumeCaloPositions
# positionsEcalBarrel = CreateVolumeCaloPositions("positionsBarrelEcal", OutputLevel = INFO)
# positionsEcalBarrel.hits.Path = "ECalBarrelCellsStep1"
# positionsEcalBarrel.positionedHits.Path = "ECalBarrelPositions"
# TopAlg += [positionsEcalBarrel]

# Use Phi-Eta segmentation in ECal barrel
from Configurables import RedoSegmentation
resegmentEcalBarrel = RedoSegmentation("ReSegmentationEcal",
                                       # old bitfield (readout)
                                       oldReadoutName = ecalBarrelReadoutName,
                                       # specify which fields are going to be altered (deleted/rewritten)
                                       oldSegmentationIds = ["module"],
                                       # new bitfield (readout), with new segmentation
                                       newReadoutName = ecalBarrelReadoutNamePhiEta,
                                       OutputLevel = INFO,
                                       # inhits = "ECalBarrelPositions",
                                       inhits = "ECalBarrelCellsStep1",
                                       outhits = "ECalBarrelCellsStep2")
TopAlg += [resegmentEcalBarrel]

createEcalBarrelCells = CreateCaloCells("CreateECalBarrelCells",
                                        doCellCalibration=False,
                                        addCellNoise=False,
                                        filterCellNoise=False,
                                        OutputLevel=INFO,
                                        hits="ECalBarrelCellsStep2",
                                        cells="ECalBarrelCells")
TopAlg += [createEcalBarrelCells]

# Create Ecal cells in endcaps
# 1. step - merge layer IDs
# 2. step - create cells
from Configurables import MergeLayers
mergelayersEcalEndcap = MergeLayers("MergeLayersEcalEndcap",
                                    # take the bitfield description from the geometry service
                                    readout = ecalEndcapReadoutName,
                                    # cells in which field should be merged
                                    identifier = identifierName,
                                    volumeName = volumeName,
                                    # how many cells to merge
                                    merge = ecalEndcapNumberOfLayersToMerge,
                                    OutputLevel = INFO)
mergelayersEcalEndcap.inhits.Path = ecalEndcapReadoutName
mergelayersEcalEndcap.outhits.Path = "mergedECalEndcapHits"
TopAlg += [mergelayersEcalEndcap]
              
createEcalEndcapCells = CreateCaloCells("CreateEcalEndcapCaloCells",
                                        doCellCalibration=True,
                                        calibTool=calibEcalEndcap,
                                        addCellNoise=False,
                                        filterCellNoise=False,
                                        OutputLevel=INFO)
createEcalEndcapCells.hits.Path="mergedECalEndcapHits"
createEcalEndcapCells.cells.Path="ECalEndcapCells"
TopAlg += [createEcalEndcapCells]

# Create Ecal cells in forward
mergelayersEcalFwd = MergeLayers("MergeLayersEcalFwd",
                                 # take the bitfield description from the geometry service
                                 readout = ecalFwdReadoutName,
                                 # cells in which field should be merged
                                 identifier = identifierName,
                                 volumeName = volumeName,
                                 # how many cells to merge
                                 merge = ecalFwdNumberOfLayersToMerge,
                                 OutputLevel = INFO)
mergelayersEcalFwd.inhits.Path = ecalFwdReadoutName
mergelayersEcalFwd.outhits.Path = "mergedECalFwdHits"
TopAlg += [mergelayersEcalFwd]

createEcalFwdCells = CreateCaloCells("CreateEcalFwdCaloCells",
                                     doCellCalibration=True,
                                     calibTool=calibEcalFwd,
                                     addCellNoise=False,
                                     filterCellNoise=False,
                                     OutputLevel=INFO)
createEcalFwdCells.hits.Path="mergedECalFwdHits"
createEcalFwdCells.cells.Path="ECalFwdCells"
TopAlg += [createEcalFwdCells]

# Create cells in HCal
# 1. step - merge hits into cells with the default readout
createHcalCells = CreateCaloCells("CreateHCaloCells",
                                  doCellCalibration=True,
                                  calibTool=calibHcells,
                                  addCellNoise = False,
                                  filterCellNoise = False,
                                  OutputLevel = INFO,
                                  hits=hcalReadoutName,
                                  cells="HCalBarrelCells")
TopAlg += [createHcalCells]

# Hcal extended barrel cells
createExtHcalCells = CreateCaloCells("CreateExtHcalCaloCells",
                                     doCellCalibration=True,
                                     calibTool=calibHcells,
                                     addCellNoise = False,
                                     filterCellNoise = False,
                                     OutputLevel = INFO,
                                     hits=extHcalReadoutName,
                                     cells="HCalExtBarrelCells")
TopAlg += [createExtHcalCells]

# Create Hcal cells in endcaps
mergelayersHcalEndcap = MergeLayers("MergeLayersHcalEndcap",
                                    # take the bitfield description from the geometry service
                                    readout = hcalEndcapReadoutName,
                                    # cells in which field should be merged
                                    identifier = identifierName,
                                    volumeName = volumeName,
                                    # how many cells to merge
                                    merge = hcalEndcapNumberOfLayersToMerge,
                                    OutputLevel = INFO)
mergelayersHcalEndcap.inhits.Path = hcalEndcapReadoutName
mergelayersHcalEndcap.outhits.Path = "mergedHCalEndcapHits"
TopAlg += [mergelayersHcalEndcap]

createHcalEndcapCells = CreateCaloCells("CreateHcalEndcapCaloCells",
                                        doCellCalibration=True,
                                        calibTool=calibHcalEndcap,
                                        addCellNoise=False,
                                        filterCellNoise=False,
                                        OutputLevel=INFO)
createHcalEndcapCells.hits.Path="mergedHCalEndcapHits"
createHcalEndcapCells.cells.Path="HCalEndcapCells"
TopAlg += [createHcalEndcapCells]

# Create Hcal cells in forward
mergelayersHcalFwd = MergeLayers("MergeLayersHcalFwd",
                                 # take the bitfield description from the geometry service
                                 readout = hcalFwdReadoutName,
                                 # cells in which field should be merged
                                 identifier = identifierName,
                                 volumeName = volumeName,
                                 # how many cells to merge
                                 merge = hcalFwdNumberOfLayersToMerge,
                                 OutputLevel = INFO)
mergelayersHcalFwd.inhits.Path = hcalFwdReadoutName
mergelayersHcalFwd.outhits.Path = "mergedHCalFwdHits"
TopAlg += [mergelayersHcalFwd]

createHcalFwdCells = CreateCaloCells("CreateHcalFwdCaloCells",
                                     doCellCalibration=True,
                                     calibTool=calibHcalFwd,
                                     addCellNoise=False, filterCellNoise=False,
                                     OutputLevel=INFO)
createHcalFwdCells.hits.Path="mergedHCalFwdHits"
createHcalFwdCells.cells.Path="HCalFwdCells"
TopAlg += [createHcalFwdCells]

# Configure output
io_svc.outputCommands = ["drop *", "keep ECalBarrelCells", "keep ECalEndcapCells", "keep ECalFwdCells", "keep HCalBarrelCells", "keep HCalExtBarrelCells", "keep HCalEndcapCells", "keep HCalFwdCells", "keep GenParticles","keep GenVertices"]
io_svc.Output = outputfile

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
