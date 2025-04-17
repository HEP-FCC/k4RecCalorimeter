from Gaudi.Configuration import INFO
from Configurables import CaloDigitizerFunc, CaloFilterFunc, CaloAddNoise2Digits
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
FilterSize = 5

NoiseEnergy = 0.0 # 1 MeV
NoiseWidth = 1.0


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
                                signalFileName=ecalBarrelSignalShapePath, # Path to a file that stores all cell IDs
                                treename="Signal_shape", # Treename to access the cell IDs in the above file
                                InputCollection=[ecalBarrelInputName], # Name of input collection
                                pulseInitTime = 0.0, # Time of pulse start [ns]
                                pulseEndTime = 775.0, # Time of pulse ending [ns]
                                pulseSamplingLength = 31, # Number of samples in the signal pulse shape
                                pulseType = "Gaussian", # Name of the signal pulse shape
                                mu = 100.0, # Mean of the Gaussian pulse shape
                                sigma = 20.0, # Sigma of the Gaussian pulse shape
                                OutputCollection=["ECalBarrelDigitized"], # Name of output collection
                                )

CaloAddNoise = CaloAddNoise2Digits("CaloAddNoise",
                                   InputCollection=["ECalBarrelDigitized"],
                                   OutputCollection=["ECalBarrelDigitizedWithNoise"],
                                   noiseEnergy=NoiseEnergy,
                                   noiseWidth=NoiseWidth,
                                   )

CaloFilter = CaloFilterFunc("CaloFilterFunc",
                            InputCollection = ["ECalBarrelDigitizedWithNoise"], # Name of input collection
                            OutputCollectionFilteredPulse = ["ECalBarrelMatchedFilterPulse"], # Name of output collection
                            OutputCollectionMatchedSampleIdx = ["ECalBarrelMatchedFilterSampleIdx"], # Name of output collection
                            OutputCollectionMatchedSampleEnergy = ["ECalBarrelMatchedFilterSampleEnergy"], # Name of output collection
                            
                            filterName = "Matched_Gaussian", # Name of the filter template
                            pulseInitTime = 0.0, # Time of pulse start [ns]
                            pulseEndTime = 775.0, # Time of pulse ending [ns]
                            pulseSamplingLength = 31, # Number of samples in the signal pulse shape
                            filterTemplateSize = 5, # Number of samples in the filter template
                            mu = 100.0, # Mean of the Gaussian pulse shape
                            sigma = 20.0, # Sigma of the Gaussian pulse shape
                            )

ApplicationMgr(TopAlg=[CaloDigitizer, 
                       CaloAddNoise,
                       CaloFilter,
                       ],
               EvtSel="NONE",
               EvtMax=Nevts,
               ExtSvc=[],
               OutputLevel=INFO,
               )