import os
import copy

from GaudiKernel.SystemOfUnits import MeV, GeV, tesla

use_pythia = False

## script to obtain energy calibration parameters using CalibrateBenchmarkMethod - to be used when shooting pions into ECal+HCa;

# Input for simulations (momentum is expected in GeV!)
# Parameters for the particle gun simulations, dummy if use_pythia is set to True
# theta from 80 to 100 degrees corresponds to -0.17 < eta < 0.17 
momentum = 50 # in GeV
#thetaMin = 90.25 # degrees
#thetaMax = 90.25 # degrees
thetaMin = 69.805 # degrees corresponds to eta = 0.36
thetaMax = 69.805 # degrees
#thetaMin = 50 # degrees
#thetaMax = 130 # degrees
pdgCode = 211 # 11 electron, 13 muon, 22 photon, 111 pi0, 211 pi+
magneticField = False

from Gaudi.Configuration import *

from Configurables import FCCDataSvc
podioevent  = FCCDataSvc("EventDataSvc")

################## Particle gun setup
_pi = 3.14159

from Configurables import GenAlg
genAlg = GenAlg()

from Configurables import  MomentumRangeParticleGun
pgun = MomentumRangeParticleGun("ParticleGun_Electron")
pgun.PdgCodes = [pdgCode]
pgun.MomentumMin = momentum * GeV
pgun.MomentumMax = momentum * GeV
pgun.PhiMin = 0
#pgun.PhiMax = 0
pgun.PhiMax = 2 * _pi
pgun.ThetaMin = thetaMin * _pi / 180.
pgun.ThetaMax = thetaMax * _pi / 180.
genAlg.SignalProvider = pgun

genAlg.hepmc.Path = "hepmc"

from Configurables import HepMCToEDMConverter
hepmc_converter = HepMCToEDMConverter()
hepmc_converter.hepmc.Path="hepmc"
genParticlesOutputName = "genParticles"
hepmc_converter.GenParticles.Path = genParticlesOutputName
hepmc_converter.hepmcStatusList = []

################## Simulation setup
# Detector geometry
from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc")
# if FCC_DETECTORS is empty, this should use relative path to working directory
path_to_detector = os.environ.get("FCCDETECTORS", "")
print(path_to_detector)
detectors_to_use=[
                    'Detector/DetFCCeeIDEA-LAr/compact/FCCee_DectMaster.xml',
                  ]
# prefix all xmls with path_to_detector
geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]
geoservice.OutputLevel = INFO

# Geant4 service
# Configures the Geant simulation: geometry, physics list and user actions

from Configurables import SimG4FullSimActions, SimG4Alg, SimG4PrimariesFromEdmTool, SimG4SaveParticleHistory
actions = SimG4FullSimActions()
# Uncomment if history from Geant4 decays is needed (e.g. to get the photons from pi0) and set actions=actions in SimG4Svc + uncomment saveHistTool in SimG4Alg
#actions.enableHistory=True
#actions.energyCut = 0.2 * GeV
#saveHistTool = SimG4SaveParticleHistory("saveHistory")

from Configurables import SimG4Svc
geantservice = SimG4Svc("SimG4Svc", detector='SimG4DD4hepDetector', physicslist="SimG4FtfpBert", actions=actions)

# Fixed seed to have reproducible results, change it for each job if you split one production into several jobs
# Mind that if you leave Gaudi handle random seed and some job start within the same second (very likely) you will have duplicates
geantservice.randomNumbersFromGaudi = False
geantservice.seedValue = 4242

# Range cut
geantservice.g4PreInitCommands += ["/run/setCut 0.1 mm"]

# Magnetic field
from Configurables import SimG4ConstantMagneticFieldTool
if magneticField == 1:
    field = SimG4ConstantMagneticFieldTool("SimG4ConstantMagneticFieldTool", FieldComponentZ=-2*tesla, FieldOn=True,IntegratorStepper="ClassicalRK4")
else:
    field = SimG4ConstantMagneticFieldTool("SimG4ConstantMagneticFieldTool",FieldOn=False)

# Geant4 algorithm
# Translates EDM to G4Event, passes the event to G4, writes out outputs via tools
# and a tool that saves the calorimeter hits

# Detector readouts
# ECAL
ecalBarrelReadoutName = "ECalBarrelEta"
ecalBarrelReadoutNamePhiEta = "ECalBarrelPhiEta"
# HCAL
hcalBarrelReadoutName = "HCalBarrelReadout"
extHcalReadoutName = "HCalExtBarrelReadout"

# Configure saving of calorimeter hits
ecalBarrelHitsName = "ECalBarrelPositionedHits"
from Configurables import SimG4SaveCalHits
saveECalBarrelTool = SimG4SaveCalHits("saveECalBarrelHits", readoutNames = [ecalBarrelReadoutName])
saveECalBarrelTool.CaloHits.Path = ecalBarrelHitsName

hcalBarrelHitsName = "HCalBarrelPositionedHits"
saveHCalTool = SimG4SaveCalHits("saveHCalBarrelHits", readoutNames = [hcalBarrelReadoutName])
saveHCalTool.CaloHits.Path = hcalBarrelHitsName

# next, create the G4 algorithm, giving the list of names of tools ("XX/YY")
from Configurables import SimG4PrimariesFromEdmTool
particle_converter = SimG4PrimariesFromEdmTool("EdmConverter")
particle_converter.GenParticles.Path = genParticlesOutputName

from Configurables import SimG4Alg
geantsim = SimG4Alg("SimG4Alg",
                       outputs= [saveECalBarrelTool,
                                 saveHCalTool,
                                 #saveHistTool
                       ],
                       eventProvider=particle_converter,
                       OutputLevel=INFO)

############## Digitization (Merging hits into cells, EM scale calibration)
# EM scale calibration (sampling fraction)
from Configurables import CalibrateInLayersTool
calibEcalBarrel = CalibrateInLayersTool("CalibrateECalBarrel",
                                   samplingFraction = [0.3632447480841956] * 1 + [0.13187261040190248] * 1 + [0.14349714292943705] * 1 + [0.150266118277841] * 1 + [0.15502683375826457] * 1 + [0.15954408786354762] * 1 + [0.16375302347299436] * 1 + [0.16840384714588075] * 1 + [0.17219540619311383] * 1 + [0.1755068643940401] * 1 + [0.17816980262822366] * 1 + [0.18131266048670405] * 1,
                                   readoutName = ecalBarrelReadoutName,
                                   layerFieldName = "layer")

from Configurables import CalibrateCaloHitsTool
calibHcells = CalibrateCaloHitsTool("CalibrateHCal", invSamplingFraction="35.5388") 
# HCal at EM scale 30.6681
# HCal at HAD scale 35.5388

# Create cells in ECal barrel
# 1. step - merge hits into cells with Eta and module segmentation (phi module is a 'physical' cell i.e. lead + LAr + PCB + LAr +lead)
# 2. step - rewrite the cellId using the Eta-Phi segmentation (merging several modules into one phi readout cell). Add noise at this step if you derived the noise already assuming merged cells
from Configurables import CreateCaloCells
createEcalBarrelCellsStep1 = CreateCaloCells("CreateECalBarrelCellsStep1",
                               doCellCalibration=True,
                               calibTool = calibEcalBarrel,
                               addCellNoise=False, filterCellNoise=False,
                               addPosition=True,
                               OutputLevel=INFO,
                               hits=ecalBarrelHitsName,
                               cells="ECalBarrelCellsStep1")

## Use Phi-Theta segmentation in ECal barrel
from Configurables import RedoSegmentation
resegmentEcalBarrel = RedoSegmentation("ReSegmentationEcal",
                             # old bitfield (readout)
                             oldReadoutName = ecalBarrelReadoutName,
                             # specify which fields are going to be altered (deleted/rewritten)
                             oldSegmentationIds = ["module"],
                             # new bitfield (readout), with new segmentation
                             newReadoutName = ecalBarrelReadoutNamePhiEta,
                             OutputLevel = INFO,
                             inhits = "ECalBarrelCellsStep1",
                             outhits = "ECalBarrelCellsStep2")

# cells without noise
EcalBarrelCellsName = "ECalBarrelCells"
createEcalBarrelCells = CreateCaloCells("CreateECalBarrelCells",
                               doCellCalibration=False,
                               addCellNoise=False, filterCellNoise=False,
                               OutputLevel=INFO,
                               hits="ECalBarrelCellsStep2",
                               cells=EcalBarrelCellsName)
cell_creator_to_use = createEcalBarrelCells


# Create cells in HCal
# 1. step - merge hits into cells with the default readout
HcalBarrelCellsName = "HCalBarrelCells"
from Configurables import CreateCaloCells
createHcalBarrelCells = CreateCaloCells("CreateHCalCells",
                               doCellCalibration=True,
                               calibTool=calibHcells,
                               addCellNoise = False, 
                               filterCellNoise = False,
                               addPosition=True,
                               OutputLevel = INFO,
                               hits=hcalBarrelHitsName,
                               cells=HcalBarrelCellsName)

from Configurables import CalibrateBenchmarkMethod
benchmark_calib = CalibrateBenchmarkMethod("CalibrateBenchmarkMethod",
                                      readoutNames=[ecalBarrelReadoutName, hcalBarrelReadoutName],
                                      energy=momentum,
                                      ECalSystemID=4,
                                      HCalSystemID=8,
                                      numLayersECal=12,
                                      firstLayerHCal=0,
                                      OutputLevel=INFO)
benchmark_calib.ecalBarrelCells.Path = EcalBarrelCellsName
benchmark_calib.hcalBarrelCells.Path = HcalBarrelCellsName 

THistSvc().Output = ["rec DATAFILE='benchmark_calibration_output_pMin_"+str(momentum*1000)+"_thetaMin_"+str(thetaMin)+".root' TYP='ROOT' OPT='RECREATE'"]
THistSvc().PrintAll=True
THistSvc().AutoSave=True
THistSvc().AutoFlush=False
THistSvc().OutputLevel=INFO

#CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
geantsim.AuditExecute = True
benchmark_calib.AuditExecute = True

import uuid
from Configurables import PodioOutput
### PODIO algorithm
out = PodioOutput("out", OutputLevel=WARNING)
out.outputCommands = ["keep *", "drop ECalBarrelHits", "drop ECalBarrelCellsStep*", "drop ECalBarrelPositionedHits", "drop HCalBarrelHits", "drop HCalBarrelPositionedHits"]
#out.filename = "fccee_caloBenchmarkCalib_%ideg_%igev_%s.root" % (thetaMin, momentum, uuid.uuid4().hex[0:16])
out.filename = "fccsw_output_pdgID_211_pMin_%i_thetaMin_%i.root" % (momentum*1000, thetaMin)

from Configurables import EventCounter
event_counter = EventCounter('event_counter')
event_counter.Frequency = 10

from Configurables import ApplicationMgr
ApplicationMgr(
    TopAlg = [
              event_counter,
              genAlg,
              hepmc_converter,
              geantsim,
              createEcalBarrelCellsStep1,
              resegmentEcalBarrel,
              #createEcalBarrelCells,
              cell_creator_to_use,
              createHcalBarrelCells,
              benchmark_calib, 
              out
              ],
                EvtSel = 'NONE',
                EvtMax = 10,
                # order is important, as GeoSvc is needed by G4SimSvc
                ExtSvc = [podioevent, geoservice, geantservice, audsvc],
                OutputLevel = WARNING
)
