#
# Do cell reconstruction for FCChh.
#

from Gaudi.Configuration import *

inputfile = 'FCChh_sim_e.root'
outputfile = 'FCChh_e_cells.root'

ExtSvc = []
TopAlg = []

# ECAL readouts
ecalBarrelReadoutName = "ECalBarrelEta"
ecalBarrelReadoutNamePhiEta = "ECalBarrelPhiEta"
ecalEndcapReadoutName = "EMECPhiEta"
ecalFwdReadoutName = "EMFwdPhiEta"
# HCAL readouts
hcalBarrelReadoutName = "HCalBarrelReadout"
hcalExtBarrelReadoutName = "HCalExtBarrelReadout"
hcalEndcapReadoutName = "HECPhiEta"
hcalFwdReadoutName = "HFwdPhiEta"

# Input/Output handling
from k4FWCore import IOSvc
from Configurables import EventDataSvc
io_svc = IOSvc("IOSvc")
io_svc.Input = inputfile
io_svc.Output = outputfile
ExtSvc += [EventDataSvc("EventDataSvc")]


from Configurables import GeoSvc
detectors_to_use=['file:$K4GEO/FCChh/compact/FCChhBaseline/FCChh_DectEmptyMaster.xml',
                  # Tracker disabled to save cpu time
                  #'file:$K4GEO/FCChh/compact/FCChhBaseline/TrackerTkLayout/Tracker.xml',
                  'file:$K4GEO/FCChh/compact/FCChhBaseline/ECalInclined/FCChh_ECalBarrel_withCryostat.xml',
                  'file:$K4GEO/FCChh/compact/FCChhBaseline/HCalTile/FCChh_HCalBarrel_TileCal.xml',
                  'file:$K4GEO/FCChh/compact/FCChhBaseline/HCalTile/FCChh_HCalExtendedBarrel_TileCal.xml',
                  'file:$K4GEO/FCChh/compact/FCChhBaseline/CalDiscs/Endcaps_coneCryo.xml',
                  'file:$K4GEO/FCChh/compact/FCChhBaseline/CalDiscs/Forward_coneCryo.xml',
                  ]

geoservice = GeoSvc("GeoSvc", detectors = detectors_to_use, OutputLevel = WARNING)
ExtSvc.append (geoservice)


def makeCells (readoutName, postoolcls):
    from Configurables import CreatePositionedCaloCells
    postool = postoolcls (readoutName + 'Positions',
                          readoutName = readoutName)
    alg = CreatePositionedCaloCells (readoutName + 'Cells',
                                     addCellNoise = False,
                                     doCellCalibration = False,
                                     positionsTool = postool,
                                     hits = readoutName,
                                     cells = readoutName + 'Cells',
                                     links = readoutName + 'CellsSimCaloHitLink')
    TopAlg.append (alg)
    return


from Configurables import CellPositionsECalBarrelTool, CellPositionsHCalBarrelNoSegTool, CellPositionsCaloDiscsTool, CellPositionsTailCatcherTool 
makeCells (ecalBarrelReadoutNamePhiEta, CellPositionsECalBarrelTool)
makeCells (ecalEndcapReadoutName, CellPositionsCaloDiscsTool)
makeCells (ecalFwdReadoutName, CellPositionsCaloDiscsTool)
makeCells (hcalBarrelReadoutName, CellPositionsHCalBarrelNoSegTool)
makeCells (hcalExtBarrelReadoutName, CellPositionsHCalBarrelNoSegTool)
makeCells (hcalEndcapReadoutName, CellPositionsCaloDiscsTool)
makeCells (hcalFwdReadoutName, CellPositionsCaloDiscsTool)

# Configure output
io_svc.outputCommands = ["keep *","drop ECalBarrelCells","drop ECalEndcapCells","drop ECalFwdCells","drop HCalBarrelCells", "drop HCalExtBarrelCells", "drop HCalEndcapCells", "drop HCalFwdCells"]

#CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
ExtSvc.append (audsvc)

from k4FWCore import ApplicationMgr
ApplicationMgr(
    TopAlg = TopAlg,
    EvtSel = 'NONE',
    EvtMax = -1,
    ExtSvc = ExtSvc,
)
