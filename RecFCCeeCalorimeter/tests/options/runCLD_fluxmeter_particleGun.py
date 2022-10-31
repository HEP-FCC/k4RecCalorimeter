import os
from Gaudi.Configuration import INFO, DEBUG
from GaudiKernel.PhysicalConstants import pi
from GaudiKernel.SystemOfUnits import GeV, tesla

from Configurables import ApplicationMgr
ApplicationMgr().EvtSel = 'NONE'
ApplicationMgr().EvtMax = 2
ApplicationMgr().OutputLevel = INFO
ApplicationMgr().StopOnSignal = True
ApplicationMgr().ExtSvc += ['RndmGenSvc']

from Configurables import FCCDataSvc
## Data service
podioevent = FCCDataSvc("EventDataSvc")
podioevent.input = "p8_ee_Zqq_ecm91.root"
from Configurables import PodioInput
podioinput = PodioInput("PodioReader",
                        collections=["GenParticles"])
ApplicationMgr().ExtSvc += [podioevent]
ApplicationMgr().TopAlg += [podioinput]

# from Configurables import MomentumRangeParticleGun
# guntool = MomentumRangeParticleGun()
# guntool.ThetaMin = 0.
# guntool.ThetaMax = pi
# guntool.PhiMin = 0.
# guntool.PhiMax = 2. * pi
# guntool.MomentumMin = 2. * GeV
# guntool.MomentumMax = 2. * GeV
# guntool.PdgCodes = [2112]  # 11 electron, 13 muon, 22 photon, 111 pi0, 211 pi+
#
# from Configurables import GenAlg
# gen = GenAlg()
# gen.SignalProvider = guntool
# gen.hepmc.Path = "hepmc"
# # gen.OutputLevel = DEBUG
# ApplicationMgr().TopAlg += [gen]
#
# from Configurables import HepMCToEDMConverter
# # Reads an HepMC::GenEvent from the data service and writes a collection of
# # EDM Particles
# hepmc_converter = HepMCToEDMConverter("Converter")
# hepmc_converter.hepmc.Path = "hepmc"
# hepmc_converter.GenParticles.Path = "GenParticles"
# ApplicationMgr().TopAlg += [hepmc_converter]


################## Simulation setup
# Detector geometry
from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc")
# if FCC_DETECTORS is empty, this should use relative path to working directory
path_to_detector = os.environ.get("FCCDETECTORS", "")
detectors_to_use = [
    # 'Detector/DetFCCeeIDEA-LAr/compact/FCCee_DectMaster.xml',
    'Detector/DetFCCeeCLD/compact/FCCee_o2_v02/FCCee_o2_v02.xml',
]
# prefix all xmls with path_to_detector
geoservice.detectors = [os.path.join(path_to_detector, _det) for _det
                        in detectors_to_use]
geoservice.OutputLevel = INFO
ApplicationMgr().ExtSvc += [geoservice]


# Geant4 service
# Configures the Geant simulation: geometry, physics list and user actions
from Configurables import SimG4Svc
# giving the names of tools will initialize the tools of that type
geantservice = SimG4Svc("SimG4Svc")
geantservice.detector = "SimG4DD4hepDetector"
geantservice.physicslist = "SimG4FtfpBert"
geantservice.actions = "SimG4FullSimActions"
# Fixed seed to have reproducible results, change it for each job if you split
# one production into several jobs
# Mind that if you leave Gaudi handle random seed and some job start within the
# same second (very likely) you will have duplicates
geantservice.randomNumbersFromGaudi = False
geantservice.seedValue = 42424
# Range cut
geantservice.g4PreInitCommands += ["/run/setCut 0.1 mm"]
ApplicationMgr().ExtSvc += [geantservice]


# Magnetic field
from Configurables import SimG4ConstantMagneticFieldTool
field = SimG4ConstantMagneticFieldTool("SimG4ConstantMagneticFieldTool")
field.FieldComponentZ = -2 * tesla
field.FieldOn = True
field.IntegratorStepper = "ClassicalRK4"

# Geant4 algorithm
# Translates EDM to G4Event, passes the event to G4, writes out outputs via
# tools and a tool that saves the calorimeter hits


# Detector readouts
from Configurables import SimG4Alg
from Configurables import SimG4SaveFluxHits
saveFluxTool = SimG4SaveFluxHits("saveFluxHits")
saveFluxTool.readoutNames = ["FluxPhiTheta"]
# saveFluxTool.FluxHits.Path = "FluxPositionedHits"
SimG4Alg("SimG4Alg").outputs += [saveFluxTool]

# next, create the G4 algorithm, giving the list of names of tools ("XX/YY")
from Configurables import SimG4PrimariesFromEdmTool
particle_converter = SimG4PrimariesFromEdmTool("EdmConverter")
particle_converter.GenParticles.Path = "GenParticles"

from Configurables import SimG4Alg
geantsim = SimG4Alg("SimG4Alg")
geantsim.eventProvider = particle_converter
ApplicationMgr().TopAlg += [geantsim]

################ Output
from Configurables import PodioOutput
out = PodioOutput("out")
out.outputCommands = ["keep *"]
import uuid
out.filename = "output_fullCalo_SimAndDigi_" + uuid.uuid4().hex[:12] + ".root"
# out.filename = "output_fullCalo_SimAndDigi.root"
ApplicationMgr().TopAlg += [out]

#CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
# gen.AuditExecute = True
# hepmc_converter.AuditExecute = True
geantsim.AuditExecute = True
out.AuditExecute = True

from Configurables import EventCounter
event_counter = EventCounter('event_counter')
event_counter.Frequency = 1000
ApplicationMgr().TopAlg += [event_counter]
