from Gaudi.Configuration import INFO
from Configurables import CaloDigitizerFunc, CaloFilterFunc
from k4FWCore import ApplicationMgr, IOSvc

Nevts = 10                              # -1 means all events
SamplingInterval = 25.0                 # ns
ecalBarrelInputName = "ECalBarrelModuleThetaMerged"  # name of the ECal barrel readout in input file
ecalBarrelSignalShapePath = "SignalPulseShapes.root"  # path to the root file with the signal pulse shapes


io_svc = IOSvc()

io_svc.Input = "OutDDSim/sample_10GeV_theta90_theta90.root"

io_svc.Output = "output_k4test_exampledata_transformer.root"
# The collections that we don't drop will also be present in the output file
io_svc.outputCommands = ["drop Lumi*", 
                         "drop Vertex*", 
                         "drop DriftChamber_simHits*", 
                         "drop MuonTagger*", 
                         "drop *SiWr*", 
                         "drop ECalEndcap*", 
                         "drop HCal*"]

CaloDigitizer = CaloDigitizerFunc("CaloDigitizerFunc",
                                signalFileName=ecalBarrelSignalShapePath,
                                treename="Signal_shape",
                                InputCollection=[ecalBarrelInputName],
                                samplingInterval = SamplingInterval,
                                OutputCollection=["ECalBarrelDigitized"],)

CaloFilter = CaloFilterFunc("CaloFilterFunc",
                            InputCollection = ["ECalBarrelDigitized"],
                            OutputCollection = ["ECalBarrelFiltered"],
                            filterName = "MatchedDirac")

ApplicationMgr(TopAlg=[CaloDigitizer, 
                       CaloFilter,
                       ],
               EvtSel="NONE",
               EvtMax=Nevts,
               ExtSvc=[],
               OutputLevel=INFO,
               )