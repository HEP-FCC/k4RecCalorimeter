from Gaudi.Configuration import *
from Configurables import ApplicationMgr

from Configurables import k4DataSvc
podioevent = k4DataSvc("EventDataSvc", input="/afs/cern.ch/user/s/sungwon/private/KEY4HEP/2025/for_digitization/k4geo/example/testIDEA_o1_v03_100Evt.root")

from Configurables import PodioInput
podioinput = PodioInput("PodioInput", collections = ["DRcaloSiPMreadoutTimeStruct", "DRcaloSiPMreadoutRawHit", "DRcaloSiPMreadoutWaveLen"], OutputLevel = DEBUG)

from Configurables import SimulateSiPMwithOpticalPhoton
digi = SimulateSiPMwithOpticalPhoton("SimulateSiPMwithOpticalPhoton",
    OutputLevel=DEBUG,
    wavelength = [
        900., 850., 800., 750., 725.,
        700., 675., 650., 625., 600.,
        590., 580., 570., 560., 550.,
        540., 530., 520., 510., 500.,
        490., 480., 470., 460., 450.,
        440., 430., 420., 400., 350.,
        300.
    ],
    sipmEfficiency = [
        0.02, 0.025, 0.045, 0.06, 0.0675,
        0.075, 0.0925, 0.11, 0.125, 0.14,
        0.146, 0.152, 0.158, 0.164, 0.17,
        0.173, 0.176, 0.178, 0.179, 0.18,
        0.181, 0.182, 0.183, 0.184, 0.18,
        0.173, 0.166, 0.158, 0.15, 0.12,
        0.05
    ]
)

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
