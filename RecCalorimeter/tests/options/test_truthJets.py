#
# IMPORTS
#
from Gaudi.Configuration import INFO

#
# SETTINGS
#
inputfile = "truthJets_sim.root"   # input file produced with ddsim
outputfile = "truthJets_rec.root"  # output file to be produced

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


# Input/Output handling
from k4FWCore import IOSvc
from Configurables import EventDataSvc
io_svc = IOSvc("IOSvc")
io_svc.Input = inputfile
io_svc.Output = outputfile
ExtSvc += [EventDataSvc("EventDataSvc")]


# create truth jets
from Configurables import CreateTruthJet
from Configurables import FilterTruthParticlesForGenJets

filterMCParticles = FilterTruthParticlesForGenJets(
    "filterMCParticles",
    InputMCParticleCollection=["MCParticles"],
    OutputMCParticleCollection=["FilteredMCParticles"]
)
TopAlg += [filterMCParticles]

createJets = CreateTruthJet(
    "createTruthJets",
    InputMCParticleCollection=["FilteredMCParticles"],
    OutputJetCollection=["TruthJets"],
    OutputAssociationsCollection=["TruthJetsAssociations"],
    MinPt=5
)
TopAlg += [createJets]



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

# for algo in applicationMgr.TopAlg:
#     algo.AuditExecute = True
