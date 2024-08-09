import os
from Gaudi.Configuration import *
from GaudiKernel import SystemOfUnits as units
from GaudiKernel import PhysicalConstants as constants
from GaudiKernel.SystemOfUnits import MeV, GeV, tesla

from Configurables import ApplicationMgr
ApplicationMgr().EvtSel = 'NONE'
ApplicationMgr().EvtMax = 2
ApplicationMgr().OutputLevel = INFO
ApplicationMgr().StopOnSignal = True
ApplicationMgr().ExtSvc += ['RndmGenSvc']

from Configurables import k4DataSvc
## Data service
podioevent = k4DataSvc("EventDataSvc")
ApplicationMgr().ExtSvc += [podioevent]

from Configurables import MomentumRangeParticleGun
guntool = MomentumRangeParticleGun()
guntool.ThetaMin = 45 * constants.pi / 180.
guntool.ThetaMax = 135 * constants.pi / 180.
guntool.PhiMin = 0.
guntool.PhiMax = 2. * constants.pi
guntool.MomentumMin = 10. *units.GeV
guntool.MomentumMax = 10. *units.GeV
guntool.PdgCodes = [11] # 11 electron, 13 muon, 22 photon, 111 pi0, 211 pi+

from Configurables import GenAlg
gen = GenAlg()
gen.SignalProvider=guntool
gen.hepmc.Path = "hepmc"
ApplicationMgr().TopAlg += [gen]

from Configurables import HepMCToEDMConverter
## reads an HepMC::GenEvent from the data service and writes a collection of EDM Particles
hepmc_converter = HepMCToEDMConverter("Converter")
hepmc_converter.hepmc.Path="hepmc"
hepmc_converter.GenParticles.Path="GenParticles"
ApplicationMgr().TopAlg += [hepmc_converter]



################## Simulation setup
# Detector geometry
from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc")
# if FCC_DETECTORS is empty, this should use relative path to working directory
path_to_detectors = os.environ.get("K4GEO", "")
detectors = [
        'FCCee/ALLEGRO/compact/ALLEGRO_o1_v01/ALLEGRO_o1_v01.xml'
]
# prefix all xmls with path_to_detectors
for det in detectors:
    geoservice.detectors += [os.path.join(path_to_detectors, det)]
geoservice.OutputLevel = INFO
ApplicationMgr().ExtSvc += [geoservice]


# Geant4 service
# Configures the Geant simulation: geometry, physics list and user actions
from Configurables import SimG4Svc
# giving the names of tools will initialize the tools of that type
geantservice = SimG4Svc("SimG4Svc")
geantservice.detector =     "SimG4DD4hepDetector"
geantservice.physicslist =  "SimG4FtfpBert"
geantservice.actions =      "SimG4FullSimActions"
# Fixed seed to have reproducible results, change it for each job if you split one production into several jobs
# Mind that if you leave Gaudi handle random seed and some job start within the same second (very likely) you will have duplicates
geantservice.randomNumbersFromGaudi = False
geantservice.seedValue = 4242
# Range cut
geantservice.g4PreInitCommands += ["/run/setCut 0.1 mm"]
ApplicationMgr().ExtSvc += [geantservice]



# Magnetic field
from Configurables import SimG4ConstantMagneticFieldTool
field = SimG4ConstantMagneticFieldTool("SimG4ConstantMagneticFieldTool")
field.FieldComponentZ = -2 * units.tesla
field.FieldOn = True
field.IntegratorStepper="ClassicalRK4"

# Geant4 algorithm
# Translates EDM to G4Event, passes the event to G4, writes out outputs via tools
# and a tool that saves the calorimeter hits


from Configurables import CreateTruthJet
from Configurables import FilterTruthParticlesForGenJets

filterMCParticles = FilterTruthParticlesForGenJets(
    "filterMCParticles",
    InputMCParticleCollection=["GenParticles"],
    OutputMCParticleCollection=["FilteredGenParticles"]
)
ApplicationMgr().TopAlg += [filterMCParticles]

createJets = CreateTruthJet(
    "createTruthJets",
    InputMCParticleCollection=["FilteredGenParticles"],
    OutputJetCollection=["TruthJets"],
    OutputAssociationsCollection=["TruthJetsAssociations"],
    MinPt=5
)
ApplicationMgr().TopAlg += [createJets]


#CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
gen.AuditExecute = True
hepmc_converter.AuditExecute = True
#out.AuditExecute = True

