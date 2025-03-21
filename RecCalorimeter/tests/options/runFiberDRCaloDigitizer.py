from Gaudi.Configuration import *
from Configurables import ApplicationMgr

from Configurables import k4DataSvc
podioevent = k4DataSvc("EventDataSvc", input="/afs/cern.ch/user/s/sungwon/private/KEY4HEP/2025/for_digitization/k4geo/example/testIDEA_o1_v03_100Evt.root")

from Configurables import PodioInput
podioinput = PodioInput("PodioInput", collections = ["DRcaloSiPMreadoutTimeStruct", "DRcaloSiPMreadoutRawHit"], OutputLevel = DEBUG)

from Configurables import FiberDRCaloDigitizer
digi = FiberDRCaloDigitizer("FiberDRCaloDigitizer", OutputLevel=DEBUG)

from Configurables import PodioOutput
podiooutput = PodioOutput("PodioOutput", filename = "testDigiIDEA_o1_v03_100Evt.root", OutputLevel = DEBUG)
podiooutput.outputCommands = ["keep *"]

ApplicationMgr(
    TopAlg = [
        podioinput,
        digi,
        podiooutput
    ],
    EvtSel = 'NONE',
    EvtMax = 100,
    ExtSvc = [podioevent]
)
