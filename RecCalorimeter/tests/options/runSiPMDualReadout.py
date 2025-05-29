from Gaudi.Configuration import *
from Configurables import ApplicationMgr

from Configurables import k4DataSvc
dataservice = k4DataSvc("EventDataSvc", input="testIDEA_o1_v03.root")

from Configurables import PodioInput
podioinput = PodioInput("PodioInput",
    collections = [
        "DRcaloSiPMreadout_scint",
        "DRcaloSiPMreadout_scintContributions",
        "DRcaloSiPMreadoutSimHit",
        "DRcaloSiPMreadoutTimeStruct",
        "DRcaloSiPMreadoutWaveLen"
    ],
    OutputLevel = DEBUG
)

from Configurables import SimulateSiPMwithEdep
sipmEdep = SimulateSiPMwithEdep("SimulateSiPMwithEdep",
    OutputLevel=DEBUG,
    inputHitCollection = "DRcaloSiPMreadout_scint",
    outputHitCollection = "DRcaloSiPMreadoutDigiHit_scint",
    outputTimeStructCollection = "DRcaloSiPMreadoutDigiWaveform_scint",
    # wavelength in nm (decreasing order)
    wavelength = [
        900., 850., 800., 750., 725.,
        700., 675., 650., 625., 600.,
        590., 580., 570., 560., 550.,
        540., 530., 520., 510., 500.,
        490., 480., 470., 460., 450.,
        440., 430., 420., 400., 350.,
        300.
    ],
    # Hamamatsu S14160-1310PS
    sipmEfficiency = [
        0.02, 0.025, 0.045, 0.06, 0.0675,
        0.075, 0.0925, 0.11, 0.125, 0.14,
        0.146, 0.152, 0.158, 0.164, 0.17,
        0.173, 0.176, 0.178, 0.179, 0.18,
        0.181, 0.182, 0.183, 0.184, 0.18,
        0.173, 0.166, 0.158, 0.15, 0.12,
        0.05
    ],
    # Kuraray SCSF-78
    scintSpectrum = [
        0., 0., 0., 0., 0.,
        0., 0., 0.0003, 0.0008, 0.0032,
        0.0057, 0.0084, 0.0153, 0.0234, 0.0343,
        0.0604, 0.0927, 0.1398, 0.2105, 0.2903,
        0.4122, 0.5518, 0.7086, 0.8678, 1.,
        0.8676, 0.2311, 0.0033, 0.0012, 0.,
        0.
    ],
    # Kuraray SCSF-78
    absorptionLength = [
        2.714, 3.619, 5.791, 4.343, 7.896,
        5.429, 36.19, 17.37, 36.19, 5.429,
        13., 14.5, 16., 18., 16.5,
        17., 14., 16., 15., 14.5,
        13., 12., 10., 8., 7.238,
        4., 1.2, 0.5, 0.2, 0.2,
        0.1
    ],
    # Kodak Wratten 9
    filterEfficiency = [
        0.903, 0.903, 0.903, 0.903, 0.903,
        0.903, 0.902, 0.901, 0.898, 0.895,
        0.893, 0.891, 0.888, 0.883, 0.87,
        0.838, 0.76, 0.62, 0.488, 0.345,
        0.207, 0.083, 0.018, 0., 0.,
        0., 0., 0., 0., 0.,
        0.
    ],
    # empirical value to keep (scint npe / ceren npe) =~ 5
    scintYield = 2.5
)

from Configurables import SimulateSiPMwithOpticalPhoton
sipmOptical = SimulateSiPMwithOpticalPhoton("SimulateSiPMwithOpticalPhoton",
    OutputLevel=DEBUG,
    inputHitCollection = "DRcaloSiPMreadoutSimHit",
    inputTimeStructCollection = "DRcaloSiPMreadoutTimeStruct",
    inputWavlenCollection = "DRcaloSiPMreadoutWaveLen",
    outputHitCollection = "DRcaloSiPMreadoutDigiHit",
    outputTimeStructCollection = "DRcaloSiPMreadoutDigiWaveform",
    # wavelength in nm (decreasing order)
    wavelength = [
        900., 850., 800., 750., 725.,
        700., 675., 650., 625., 600.,
        590., 580., 570., 560., 550.,
        540., 530., 520., 510., 500.,
        490., 480., 470., 460., 450.,
        440., 430., 420., 400., 350.,
        300.
    ],
    # Hamamatsu S14160-1310PS
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
podiooutput = PodioOutput("PodioOutput", filename = "testIDEA_o1_v03_digi.root", OutputLevel = DEBUG)
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
        sipmEdep,
        sipmOptical,
        podiooutput
    ],
    EvtSel = 'NONE',
    EvtMax = 10,
    ExtSvc = [rndmEngine,rndmGenSvc,dataservice]
)
