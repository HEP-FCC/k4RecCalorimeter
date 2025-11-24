from Gaudi.Configuration import INFO
from Configurables import CaloDigitizerFunc, CaloFilterFunc, CaloAddNoise2Digits, CaloWhitening
from k4FWCore import ApplicationMgr, IOSvc


#Additional information on this algorithm can be found in https://indico.cern.ch/event/1580025/contributions/6686602/attachments/3133485/5559196/FCCDigitization4BNLWorkshopEndOfWeekUpdate.pdf

Nevts = 10                              # -1 means all events
DigitInitTime = 0.0                     # Defining the initial time for digitization
DigitEndTime = 775.0                  # time range of the digitization
PulseSampleLen = 31                   # number of samples in the signal pulse shape
ecalBarrelInputName = "ECalBarrelModuleThetaMerged"  # name of the ECal barrel readout in input file
PulseShapeName = "Gaussian"             # name of the signal pulse shape
GaussianMean = 100.0                    # Mean of the Gaussian pulse shape used to represent a "signal"
GaussianSigma = 20.0                    # Standard deviation of the Gaussian pulse shape used to represent a "signal"
FilterSize = 5                          # Size of the matched filter to consider

NumberOfNoiseSamplesToSimulate = 2000   # Number of noise samples to simulate
NoiseEnergy = 0.001                     # Noise mean value to consider (simulation is Gaussian and units are in GeV)
NoiseWidth = 0.0001                     # Width of noise to consider (units in GeV)
NoiseSampleSimulationFName="NoiseInfoTest_New_Modded4InvCorr.root" # File name to save the correlation matrix after simulating noise

WhiteningFilterName2Apply = "ZCA"       # Algorithm to calculate whitening filter

io_svc = IOSvc()

io_svc.Input = "OutDDSim/sample_10GeV_theta90_theta90.root" # Input filename

io_svc.Output = "output_k4test_exampledata_transformer.root" # Output filename

# The collections that we don't drop will also be present in the output file
io_svc.outputCommands = ["drop Lumi*", 
                         "drop Vertex*", 
                         "drop DriftChamber_simHits*", 
                         "drop MuonTagger*", 
                         "drop *SiWr*", 
                         "drop ECalEndcap*", 
                         "drop HCal*"]

CaloDigitizer = CaloDigitizerFunc("CaloDigitizerFunc",
                                InputCollection=[ecalBarrelInputName], # Name of input collection
                                pulseInitTime = DigitInitTime, # Time of pulse start [ns]
                                pulseEndTime = DigitEndTime, # Time of pulse ending [ns]
                                pulseSamplingLength = PulseSampleLen, # Number of samples in the signal pulse shape
                                pulseType = PulseShapeName, # Name of the signal pulse shape
                                mu = GaussianMean, # Mean of the Gaussian pulse shape
                                sigma = GaussianSigma, # Sigma of the Gaussian pulse shape
                                OutputCollection=["ECalBarrelDigitized"], # Name of output collection
                                )

CaloAddNoise = CaloAddNoise2Digits("CaloAddNoise",
                                   InputCollection=["ECalBarrelDigitized"],
                                   OutputCollection=["ECalBarrelDigitizedWithNoise"],
                                   noiseEnergy=NoiseEnergy,
                                   noiseWidth=NoiseWidth,
                                   noiseSimSamples=NumberOfNoiseSamplesToSimulate,
                                   noiseInfoFileName=NoiseSampleSimulationFName,
                                   pulseSamplingLength = PulseSampleLen,
                                   )

CaloWhiteningFilter = CaloWhitening("CaloWhitening",
                                   InputCollection=["ECalBarrelDigitizedWithNoise"],
                                   OutputCollection=["ECalBarrelWhitenedDigits"],
                                   noiseInfoFileName = NoiseSampleSimulationFName,
                                   invCorrMatName = "InvertedCorrelationMatrix",
                                   muVecName = "NoiseSampleMat_MeansVector",
                                   filterName = WhiteningFilterName2Apply,
                                   )

CaloFilter = CaloFilterFunc("CaloFilterFunc",
                            InputCollection = ["ECalBarrelWhitenedDigits"], # Name of input collection
                            OutputCollectionFilteredPulse = ["ECalBarrelMatchedFilterPulse"], # Name of output collection
                            OutputCollectionMatchedSampleIdx = ["ECalBarrelMatchedFilterSampleIdx"], # Name of output collection
                            OutputCollectionMatchedSampleEnergy = ["ECalBarrelMatchedFilterSampleEnergy"], # Name of output collection
                            
                            filterName = "Matched_Gaussian", # Name of the filter template
                            pulseInitTime = DigitInitTime, # Time of pulse start [ns]
                            pulseEndTime = DigitEndTime, # Time of pulse ending [ns]
                            pulseSamplingLength = PulseSampleLen, # Number of samples in the signal pulse shape
                            filterTemplateSize = FilterSize, # Number of samples in the filter template
                            mu = GaussianMean, # Mean of the Gaussian pulse shape
                            sigma = GaussianSigma, # Sigma of the Gaussian pulse shape
                            noiseInfoFileName = NoiseSampleSimulationFName,
                            invCorrMatName = "InvertedCorrelationMatrix",
                            )

ApplicationMgr(TopAlg=[
                       CaloDigitizer, 
                       CaloAddNoise,
                       CaloWhiteningFilter,
                       CaloFilter,
                       ],
               EvtSel="NONE",
               EvtMax=Nevts,
               ExtSvc=[],
               OutputLevel=INFO,
               )
