from Gaudi.Configuration import INFO
from Configurables import CaloDigitizerFunc, CaloFilterFunc
from k4FWCore import ApplicationMgr, IOSvc

Nevts = 10                              # -1 means all events
DigitInitTime = 0.0
DigitEndTime = 775.0                  # time range of the digitization
PulseSampleLen = 31                   # number of samples in the signal pulse shape
ecalBarrelInputName = "ECalBarrelModuleThetaMerged"  # name of the ECal barrel readout in input file
ecalBarrelSignalShapePath = "SignalPulseShapes.root"  # path to the root file with the signal pulse shapes
PulseShapeName = "Gaussian"             # name of the signal pulse shape
GaussianMean = 100.0
GaussianSigma = 20.0
FilterSize = 7


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
                                pulseInitTime = DigitInitTime,
                                pulseEndTime = DigitEndTime,
                                pulseSamplingLength = PulseSampleLen,
                                pulseType = PulseShapeName,
                                mu = GaussianMean,
                                sigma = GaussianSigma,
                                OutputCollection=["ECalBarrelDigitized"],)

CaloFilter = CaloFilterFunc("CaloFilterFunc",
                            InputCollection = ["ECalBarrelDigitized"],
                            OutputCollectionFilteredPulse = ["ECalBarrelMatchedFilterPulse"], 
                            OutputCollectionMatchedSampleIdx = ["ECalBarrelMatchedFilterSampleIdx"],
                            OutputCollectionMatchedSampleEnergy = ["ECalBarrelMatchedFilterSampleEnergy"],
                            
                            filterName = "Matched_%s" % PulseShapeName,
                            pulseInitTime = DigitInitTime,
                            pulseEndTime = DigitEndTime,
                            pulseSamplingLength = PulseSampleLen,
                            filterTemplateSize = FilterSize,
                            mu = GaussianMean,
                            sigma = GaussianSigma,
                            )

ApplicationMgr(TopAlg=[CaloDigitizer, 
                       CaloFilter,
                       ],
               EvtSel="NONE",
               EvtMax=Nevts,
               ExtSvc=[],
               OutputLevel=INFO,
               )