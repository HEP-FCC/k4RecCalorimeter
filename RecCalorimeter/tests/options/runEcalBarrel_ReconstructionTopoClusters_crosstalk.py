from Configurables import ApplicationMgr
from Configurables import EventCounter
from Configurables import AuditorSvc, ChronoAuditor
from Configurables import PodioOutput
from Configurables import CreateEmptyCaloCellsCollection
from Configurables import CreateCaloCellPositionsFCCee
from Configurables import CellPositionsECalBarrelModuleThetaSegTool
from Configurables import CreateCaloCells
from Configurables import CalibrateCaloHitsTool
from Configurables import CalibrateInLayersTool
from Configurables import SimG4Alg
from Configurables import SimG4PrimariesFromEdmTool
from Configurables import SimG4SaveCalHits
from Configurables import SimG4ConstantMagneticFieldTool
from Configurables import SimG4Svc
from Configurables import SimG4FullSimActions
from Configurables import SimG4SaveParticleHistory
from Configurables import GeoSvc
from Configurables import HepMCToEDMConverter
from Configurables import GenAlg
from Configurables import FCCDataSvc
from Configurables import CaloTopoClusterInputTool
from Configurables import TopoCaloNeighbours
from Configurables import TopoCaloNoisyCells
from Configurables import CaloTopoClusterFCCee
from Configurables import RewriteBitfield
from Configurables import ReadCaloCrosstalkMap
from Gaudi.Configuration import INFO
# , VERBOSE, DEBUG
# from Gaudi.Configuration import *

import os

from GaudiKernel.SystemOfUnits import GeV, tesla, mm
from GaudiKernel.PhysicalConstants import pi, halfpi, twopi
from math import cos, sin, tan

# for big productions, save significant space removing hits and cells
# however, hits and cluster cells might be wanted for small productions for detailed event displays
# also, cluster cells are needed for the MVA training
saveHits = False
saveCells = False
saveClusterCells = True
doCrosstalk = True # switch on/off the crosstalk

# Input for simulations (momentum is expected in GeV!)
# Parameters for the particle gun simulations, dummy if use_pythia is set
# to True
# theta from 80 to 100 degrees corresponds to -0.17 < eta < 0.17
# reminder: cell granularity in theta = 0.5625 degrees
# (in strips: 0.5625/4=0.14)

Nevts = 5

# particle momentum and direction
momentum = 50  # in GeV
thetaMin = 50  # degrees
thetaMax = 130  # degrees
phiMin = 0
phiMax = twopi

# particle type: 11 electron, 13 muon, 22 photon, 111 pi0, 211 pi+
pdgCode = 11

# Set to true if history from Geant4 decays is needed (e.g. to get the
# photons from pi0)
saveG4Hist = False
if (pdgCode == 111):
    saveG4Hist = True

magneticField = False


podioevent = FCCDataSvc("EventDataSvc")

# Particle gun setup

genAlg = GenAlg()
from Configurables import MomentumRangeParticleGun
pgun = MomentumRangeParticleGun("ParticleGun")
pgun.PdgCodes = [pdgCode]
pgun.MomentumMin = momentum * GeV
pgun.MomentumMax = momentum * GeV
pgun.PhiMin = phiMin
pgun.PhiMax = phiMax
pgun.ThetaMin = thetaMin * pi / 180.
pgun.ThetaMax = thetaMax * pi / 180.
genAlg.SignalProvider = pgun

genAlg.hepmc.Path = "hepmc"

# hepMC -> EDM converter
hepmc_converter = HepMCToEDMConverter()
hepmc_converter.hepmc.Path = "hepmc"
genParticlesOutputName = "genParticles"
hepmc_converter.GenParticles.Path = genParticlesOutputName
hepmc_converter.hepmcStatusList = []
hepmc_converter.OutputLevel = INFO

# Simulation setup
# Detector geometry
geoservice = GeoSvc("GeoSvc")
# if K4GEO is empty, this should use relative path to working directory
path_to_detector = os.environ.get("K4GEO", "")
print(path_to_detector)
detectors_to_use = [
    'FCCee/ALLEGRO/compact/ALLEGRO_o1_v02/ALLEGRO_o1_v02.xml'
]
# prefix all xmls with path_to_detector
geoservice.detectors = [
    os.path.join(path_to_detector, _det) for _det in detectors_to_use
]
geoservice.OutputLevel = INFO

# Geant4 service
# Configures the Geant simulation: geometry, physics list and user actions
actions = SimG4FullSimActions()

if saveG4Hist:
    actions.enableHistory = True
    actions.energyCut = 1.0 * GeV
    saveHistTool = SimG4SaveParticleHistory("saveHistory")

geantservice = SimG4Svc(
    "SimG4Svc",
    detector='SimG4DD4hepDetector',
    physicslist="SimG4FtfpBert",
    actions=actions
)

# Fixed seed to have reproducible results, change it for each job if you
# split one production into several jobs
# Mind that if you leave Gaudi handle random seed and some job start within
# the same second (very likely) you will have duplicates
geantservice.randomNumbersFromGaudi = False
geantservice.seedValue = 4242

# Range cut
geantservice.g4PreInitCommands += ["/run/setCut 0.1 mm"]

# Magnetic field
if magneticField == 1:
    field = SimG4ConstantMagneticFieldTool(
        "SimG4ConstantMagneticFieldTool",
        FieldComponentZ=-2 * tesla,
        FieldOn=True,
        IntegratorStepper="ClassicalRK4"
    )
else:
    field = SimG4ConstantMagneticFieldTool(
        "SimG4ConstantMagneticFieldTool",
        FieldOn=False
    )

# Geant4 algorithm
# Translates EDM to G4Event, passes the event to G4, writes out outputs
# via tools and a tool that saves the calorimeter hits

# Detector readouts
# ECAL
ecalBarrelReadoutName = "ECalBarrelModuleThetaMerged"

# Configure saving of calorimeter hits
ecalBarrelHitsName = "ECalBarrelPositionedHits"
saveECalBarrelTool = SimG4SaveCalHits(
    "saveECalBarrelHits",
    readoutName=ecalBarrelReadoutName,
    OutputLevel=INFO
)
saveECalBarrelTool.CaloHits.Path = ecalBarrelHitsName

# next, create the G4 algorithm, giving the list of names of tools ("XX/YY")
particle_converter = SimG4PrimariesFromEdmTool("EdmConverter")
particle_converter.GenParticles.Path = genParticlesOutputName

outputTools = [
    saveECalBarrelTool
]

if saveG4Hist:
    outputTools += [saveHistTool]

geantsim = SimG4Alg("SimG4Alg",
                    outputs=outputTools,
                    eventProvider=particle_converter,
                    OutputLevel=INFO)

# Digitization (Merging hits into cells, EM scale calibration)
# EM scale calibration (sampling fraction)
calibEcalBarrel = CalibrateInLayersTool("CalibrateECalBarrel",
                                        samplingFraction=[0.3864252122990472] * 1 + [0.13597644835735828] * 1 + [0.14520427829645913] * 1 + [0.1510076084632846] * 1 + [0.1552347580991012] * 1 + [0.159694330729184] * 1 + [0.1632954482794191] * 1 + [0.16720711037339814] * 1 + [0.17047749048884808] * 1 + [0.17461698117974286] * 1 + [0.1798984163980135] * 1 + [0.17920355117405806] * 1,
                                        readoutName=ecalBarrelReadoutName,
                                        layerFieldName="layer")

# read the crosstalk map
readCrosstalkMap = ReadCaloCrosstalkMap("ReadCrosstalkMap",
                                       fileName="https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v02/neighbours_map_barrel_thetamodulemerged.root",
                                       OutputLevel=INFO)

# merge hits into cells according to the detector segmentation
ecalBarrelCellsName = "ECalBarrelCells"
createEcalBarrelCells = CreateCaloCells("CreateECalBarrelCells",
                                        doCellCalibration=True,
                                        calibTool=calibEcalBarrel,
					crosstalksTool=readCrosstalkMap,
					addCrosstalk=doCrosstalk,
                                        addCellNoise=False,
                                        filterCellNoise=False,
                                        addPosition=True,
                                        OutputLevel=INFO,
                                        hits=ecalBarrelHitsName,
                                        cells=ecalBarrelCellsName)

# Add to Ecal barrel cells the position information
# (good for physics, all coordinates set properly)

cellPositionEcalBarrelTool = CellPositionsECalBarrelModuleThetaSegTool(
    "CellPositionsECalBarrel",
    readoutName=ecalBarrelReadoutName,
    OutputLevel=INFO
)
ecalBarrelPositionedCellsName = "ECalBarrelPositionedCells"
createEcalBarrelPositionedCells = CreateCaloCellPositionsFCCee(
    "CreateECalBarrelPositionedCells",
    OutputLevel=INFO
)
createEcalBarrelPositionedCells.positionsTool = cellPositionEcalBarrelTool
createEcalBarrelPositionedCells.hits.Path = ecalBarrelCellsName
createEcalBarrelPositionedCells.positionedHits.Path = ecalBarrelPositionedCellsName

# Empty cells for parts of calorimeter not implemented yet
createemptycells = CreateEmptyCaloCellsCollection("CreateEmptyCaloCells")
createemptycells.cells.Path = "emptyCaloCells"

# TOPO CLUSTERS PRODUCTION
createTopoInput = CaloTopoClusterInputTool("CreateTopoInput",
                                           ecalBarrelReadoutName=ecalBarrelReadoutName,
                                           ecalEndcapReadoutName="",
                                           ecalFwdReadoutName="",
                                           hcalBarrelReadoutName="",
                                           hcalExtBarrelReadoutName="",
                                           hcalEndcapReadoutName="",
                                           hcalFwdReadoutName="",
                                           OutputLevel=INFO)

createTopoInput.ecalBarrelCells.Path = ecalBarrelPositionedCellsName
createTopoInput.ecalEndcapCells.Path = "emptyCaloCells"
createTopoInput.ecalFwdCells.Path = "emptyCaloCells"
createTopoInput.hcalBarrelCells.Path = "emptyCaloCells"
createTopoInput.hcalExtBarrelCells.Path = "emptyCaloCells"
createTopoInput.hcalEndcapCells.Path = "emptyCaloCells"
createTopoInput.hcalFwdCells.Path = "emptyCaloCells"
cellPositionHcalBarrelNoSegTool = None
cellPositionHcalExtBarrelTool = None

neighboursMap = "/LAr_scripts/data/neighbours_map_barrel_thetamodulemerged.root"
noiseMap = "/LAr_scripts/data/cellNoise_map_electronicsNoiseLevel_thetamodulemerged.root"

readNeighboursMap = TopoCaloNeighbours("ReadNeighboursMap",
                                       fileName=os.environ['FCCBASEDIR'] + neighboursMap,
                                       OutputLevel=INFO)

# Noise levels per cell
readNoisyCellsMap = TopoCaloNoisyCells("ReadNoisyCellsMap",
                                       fileName=os.environ['FCCBASEDIR'] + noiseMap,
                                       OutputLevel=INFO)

createTopoClusters = CaloTopoClusterFCCee("CreateTopoClusters",
                                          TopoClusterInput=createTopoInput,
                                          # expects neighbours map from cellid->vec < neighbourIds >
                                          neigboursTool=readNeighboursMap,
                                          # tool to get noise level per cellid
                                          noiseTool=readNoisyCellsMap,
                                          # cell positions tools for all sub - systems
                                          positionsECalBarrelTool=cellPositionEcalBarrelTool,
                                          positionsHCalBarrelTool=None,
                                          # positionsHCalBarrelNoSegTool=cellPositionHcalBarrelNoSegTool,
                                          # positionsHCalExtBarrelTool=cellPositionHcalExtBarrelTool,
                                          # positionsHCalExtBarrelTool = HCalExtBcells,
                                          # positionsEMECTool = EMECcells,
                                          # positionsHECTool = HECcells,
                                          # positionsEMFwdTool = ECalFwdcells,
                                          # positionsHFwdTool = HCalFwdcells,
                                          noSegmentationHCal=False,
                                          seedSigma=4,
                                          neighbourSigma=2,
                                          lastNeighbourSigma=0,
                                          OutputLevel=INFO)
createTopoClusters.clusters.Path = "CaloTopoClusters"
createTopoClusters.clusterCells.Path = "CaloTopoClusterCells"

# Output
out = PodioOutput("out",
                  OutputLevel=INFO)

out.outputCommands = ["keep *"]
#out.outputCommands = ["keep *", "drop HCal*", "drop emptyCaloCells"]

if not saveCells:
    out.outputCommands.append("drop ECal*Cells*")
if not saveClusterCells:
    out.outputCommands.append("drop *ClusterCells*")
if not saveHits:
    out.outputCommands.append("drop ECal*Hits*")

out.filename = "./output_test_EcalBarrel_crosstalk.root"

# CPU information
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
genAlg.AuditExecute = True
hepmc_converter.AuditExecute = True
geantsim.AuditExecute = True
createEcalBarrelCells.AuditExecute = True
createEcalBarrelPositionedCells.AuditExecute = True
createTopoClusters.AuditExecute = True
out.AuditExecute = True

event_counter = EventCounter('event_counter')
event_counter.Frequency = 10

ExtSvc = [geoservice, podioevent, geantservice, audsvc]

TopAlg = [
    event_counter,
    genAlg,
    hepmc_converter,
    geantsim,
    createEcalBarrelCells,
    createEcalBarrelPositionedCells,
]

TopAlg += [
    createemptycells,
    createTopoClusters
]

TopAlg += [
    out
]

ApplicationMgr(
    TopAlg=TopAlg,
    EvtSel='NONE',
    EvtMax=Nevts,
    ExtSvc=ExtSvc,
    StopOnSignal=True,
)
