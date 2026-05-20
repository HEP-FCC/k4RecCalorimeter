from Gaudi.Configuration import *
from Configurables import EventDataSvc
from k4FWCore import ApplicationMgr, IOSvc

iosvc = IOSvc()
iosvc.Input = "IDEA_o2_v01.root"
iosvc.CollectionNames = [
    "DRBTScin",  # collections from IDEA_o2 hadronic calo
    "DRBTScinContributions",
    "DRBTCher",
    "DRBTCherContributions",
    "DRETScinLeft",
    "DRETScinLeftContributions",
    "DRETScinRight",
    "DRETScinRightContributions",
    "DRETCherLeft",
    "DRETCherLeftContributions",
    "DRETCherRight",
    "DRETCherRightContributions",
]

dataservice = EventDataSvc("EventDataSvc")

# SimulateSiPMwithContrib applied sipm digitization
# over the output format of the IDEA_o2 hadronic calorimeter.
# In this format we do not store photons but photo-electrons (fired cells)
# therefore the wavelength and sipmEfficiency is dummy:
# all photo-electrons are converted.
from Configurables import SimulateSiPMwithContrib

sipmContribBTScin = SimulateSiPMwithContrib(
    "sipmContribBTScin",
    OutputLevel=DEBUG,
    inputHitCollection="DRBTScin",
    outputHitCollection="DRBTScin_digi",
    outputLinkCollection="DRBTScin_links",
    # readoutName = "DRTcaloSiPMreadout",
    isCherenkov=False,
    # wavelength in nm (decreasing order)
    wavelength=[900.0, 300.0],
    sipmEfficiency=[1.0, 1.0],
)

sipmContribETScinLeft = SimulateSiPMwithContrib(
    "sipmContribETScinLeft",
    OutputLevel=DEBUG,
    inputHitCollection="DRETScinLeft",
    outputHitCollection="DRETScinLeft_digi",
    outputLinkCollection="DRETScinLeft_links",
    # readoutName = "DRTcaloSiPMreadout",
    isCherenkov=False,
    # wavelength in nm (decreasing order)
    wavelength=[900.0, 300.0],
    sipmEfficiency=[1.0, 1.0],
)

sipmContribETScinRight = SimulateSiPMwithContrib(
    "sipmContribETScinRight",
    OutputLevel=DEBUG,
    inputHitCollection="DRETScinRight",
    outputHitCollection="DRETScinRight_digi",
    outputLinkCollection="DRETScinRight_links",
    # readoutName = "DRTcaloSiPMreadout",
    isCherenkov=False,
    # wavelength in nm (decreasing order)
    wavelength=[900.0, 300.0],
    sipmEfficiency=[1.0, 1.0],
)

sipmContribBTCher = SimulateSiPMwithContrib(
    "sipmContribBTCher",
    OutputLevel=DEBUG,
    inputHitCollection="DRBTCher",
    outputHitCollection="DRBTCher_digi",
    outputLinkCollection="DRBTCher_links",
    # readoutName = "DRTcaloSiPMreadout",
    isCherenkov=True,
    # wavelength in nm (decreasing order)
    wavelength=[900.0, 300.0],
    sipmEfficiency=[1.0, 1.0],
)

sipmContribETCherLeft = SimulateSiPMwithContrib(
    "sipmContribETCherLeft",
    OutputLevel=DEBUG,
    inputHitCollection="DRETCherLeft",
    outputHitCollection="DRETCherLeft_digi",
    outputLinkCollection="DRETCherLeft_links",
    # readoutName = "DRTcaloSiPMreadout",
    isCherenkov=True,
    # wavelength in nm (decreasing order)
    wavelength=[900.0, 300.0],
    sipmEfficiency=[1.0, 1.0],
)

sipmContribETCherRight = SimulateSiPMwithContrib(
    "sipmContribETCherRight",
    OutputLevel=DEBUG,
    inputHitCollection="DRETCherRight",
    outputHitCollection="DRETCherRight_digi",
    outputLinkCollection="DRETCherRight_links",
    # readoutName = "DRTcaloSiPMreadout",
    isCherenkov=True,
    # wavelength in nm (decreasing order)
    wavelength=[900.0, 300.0],
    sipmEfficiency=[1.0, 1.0],
)

iosvc.Output = "IDEA_o2_v01_digi.root"
iosvc.outputCommands = ["keep *"]

from Configurables import HepRndm__Engine_CLHEP__RanluxEngine_ as RndmEngine

rndmEngine = RndmEngine(
    "RndmGenSvc.Engine",
    SetSingleton=True,
    Seeds=[2345678],  # default seed is 1234567
)

from Configurables import RndmGenSvc

rndmGenSvc = RndmGenSvc("RndmGenSvc", Engine=rndmEngine.name())

ApplicationMgr(
    TopAlg=[
        sipmContribBTScin,
        sipmContribETScinLeft,
        sipmContribETScinRight,
        sipmContribBTCher,
        sipmContribETCherLeft,
        sipmContribETCherRight,
    ],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[rndmEngine, rndmGenSvc, dataservice],
)
