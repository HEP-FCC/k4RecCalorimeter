from Gaudi.Configuration import *
from Configurables import ApplicationMgr

from Configurables import k4DataSvc
dataservice = k4DataSvc("EventDataSvc", input="IDEA_o2_v01.root")

from Configurables import PodioInput
podioinput = PodioInput("PodioInput",
    collections = [
         "DRBTScin",
         "DRBTScinContributions"
    ],
    OutputLevel = DEBUG
)

from Configurables import SimulateSiPMwithContrib
sipmContribBTScin = SimulateSiPMwithContrib("SimulateSiPMwithContrib",
    OutputLevel=DEBUG,
    inputHitCollection = "DRBTScin",
    outputHitCollection = "DRBTScin_digi",
    readoutName = "DRTcaloSiPMreadout",
    # wavelength in nm (decreasing order)
    wavelength = [900.,300.],
    sipmEfficiency = [1.0,1.0]
)

from Configurables import PodioOutput
podiooutput = PodioOutput("PodioOutput", filename = "IDEA_o2_v01_digi.root", OutputLevel = DEBUG)
podiooutput.outputCommands = ["keep *"]

from Configurables import HepRndm__Engine_CLHEP__RanluxEngine_ as RndmEngine
rndmEngine = RndmEngine('RndmGenSvc.Engine',
  SetSingleton = True,
  Seeds = [ 2345678 ] # default seed is 1234567
)

from Configurables import RndmGenSvc
rndmGenSvc = RndmGenSvc("RndmGenSvc",
  Engine = rndmEngine.name()
)

ApplicationMgr(
    TopAlg = [
        podioinput,
        sipmContribBTScin,
        podiooutput
    ],
    EvtSel = 'NONE',
    EvtMax = -1,
    ExtSvc = [rndmEngine,rndmGenSvc,dataservice]
)
