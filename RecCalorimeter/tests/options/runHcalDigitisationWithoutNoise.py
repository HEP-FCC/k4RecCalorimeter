from Gaudi.Configuration import *

from Configurables import EventDataSvc
from k4FWCore import ApplicationMgr, IOSvc

iosvc = IOSvc()
iosvc.Input = "output_hcalSim_e50GeV_eta036_10events.root"
iosvc.CollectionNames = ["HCalHits", "ExtHCalHits"]

podioevent = EventDataSvc("EventDataSvc")

from Configurables import GeoSvc

geoservice = GeoSvc(
    "GeoSvc",
    detectors=[
        "file:Detector/DetFCChhBaseline1/compact/FCChh_DectEmptyMaster.xml",
        "file:Detector/DetFCChhHCalTile/compact/FCChh_HCalBarrel_TileCal.xml",
        "file:Detector/DetFCChhHCalTile/compact/FCChh_HCalExtendedBarrel_TileCal.xml",
    ],
    OutputLevel=INFO,
)

# Configure tools for calo reconstruction
from Configurables import (
    CalibrateCaloHitsTool,
    CreateVolumeCaloPositions,
    RedoSegmentation,
    RewriteBitfield,
)

calibHcells = CalibrateCaloHitsTool("CalibrateHCal", invSamplingFraction="41.7 ")

# Use Phi-Eta segmentation in Hcal barrel
rewriteHCal = RewriteBitfield(
    "RewriteHCal",
    # old bitfield (readout)
    oldReadoutName="HCalBarrelReadout",
    # specify which fields are going to be deleted
    removeIds=["row"],
    # new bitfield (readout), with new segmentation
    newReadoutName="BarHCal_Readout_phieta",
    debugPrint=10,
    OutputLevel=INFO,
)
# clusters are needed, with deposit position and cellID in bits
rewriteHCal.inhits.Path = "HCalHits"
rewriteHCal.outhits.Path = "HCalBarrelCellsStep2"

rewriteExtHCal = RewriteBitfield(
    "RewriteExtHCal",
    # old bitfield (readout)
    oldReadoutName="HCalExtBarrelReadout",
    # specify which fields are going to be deleted
    removeIds=["row"],
    # new bitfield (readout), with new segmentation
    newReadoutName="ExtBarHCal_Readout_phieta",
    debugPrint=10,
    OutputLevel=INFO,
)
# clusters are needed, with deposit position and cellID in bits
rewriteExtHCal.inhits.Path = "ExtHCalHits"
rewriteExtHCal.outhits.Path = "HCalExtBarrelCellsStep2"

# Use Phi segmentation in Hcal barrel
resegmentHcalBarrel = RewriteBitfield(
    "ReSegmentationHcal",
    # old bitfield (readout)
    oldReadoutName="HCalBarrelReadout",
    # specify which fields are going to be altered (deleted/rewritten)
    removeIds=["eta"],
    # new bitfield (readout), with new segmentation
    newReadoutName="BarHCal_Readout_phi",
    inhits="HCalHits",
    outhits="HCalBarrelCellsStep1",
)
resegmentHcalExtBarrel = RewriteBitfield(
    "ReSegmentationExtHcal",
    # old bitfield (readout)
    oldReadoutName="HCalExtBarrelReadout",
    # specify which fields are going to be altered (deleted/rewritten)
    removeIds=["eta"],
    # new bitfield (readout), with new segmentation
    newReadoutName="ExtBarHCal_Readout_phi",
    inhits="ExtHCalHits",
    outhits="HCalExtBarrelCellsStep1",
)

from Configurables import CreateCaloCells

createHcalBarrelTiles = CreateCaloCells(
    "CreateHCalBarrelTiles",
    calibTool=calibHcells,
    doCellCalibration=True,
    addCellNoise=False,
    filterCellNoise=False,
    OutputLevel=INFO,
    hits="HCalBarrelCellsStep1",
    cells="HCalAllTiles",
)
createHcalExtBarrelTiles = CreateCaloCells(
    "CreateHCalExtBarrelTiles",
    calibTool=calibHcells,
    doCellCalibration=True,
    addCellNoise=False,
    filterCellNoise=False,
    OutputLevel=INFO,
    hits="HCalExtBarrelCellsStep1",
    cells="ExtHCalAllTiles",
)

createHcalBarrelCells = CreateCaloCells(
    "CreateHCalBarrelCells",
    calibTool=calibHcells,
    doCellCalibration=True,
    addCellNoise=False,
    filterCellNoise=False,
    OutputLevel=INFO,
    hits="HCalBarrelCellsStep2",
    cells="HCalCells",
)
createHcalExtBarrelCells = CreateCaloCells(
    "CreateHCalExtBarrelCells",
    calibTool=calibHcells,
    doCellCalibration=True,
    addCellNoise=False,
    filterCellNoise=False,
    OutputLevel=INFO,
    hits="HCalExtBarrelCellsStep2",
    cells="ExtHCalCells",
)

iosvc.Output = "output_HCalCells_digitisation_noNoise.root"
iosvc.outputCommands = ["keep *"]

# CPU information
from Configurables import AuditorSvc, ChronoAuditor

chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
createHcalBarrelCells.AuditExecute = True

ApplicationMgr(
    TopAlg=[
        resegmentHcalBarrel,
        resegmentHcalExtBarrel,
        rewriteHCal,
        rewriteExtHCal,
        createHcalBarrelTiles,
        createHcalExtBarrelTiles,
        createHcalBarrelCells,
        createHcalExtBarrelCells,
    ],
    EvtSel="NONE",
    EvtMax=1,
    ExtSvc=[podioevent, geoservice],
)
