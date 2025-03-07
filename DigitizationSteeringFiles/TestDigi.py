from Configurables import GeoSvc
from Configurables import ApplicationMgr
from Configurables import CreateFCCeeCaloXTalkNeighbours
import os
from Gaudi.Configuration import INFO, DEBUG

inputfile = "ALLEGRO_sim_ee_z_qq.root"  # input file produced with ddsim
Nevts = -1                              # -1 means all events
SamplingInterval = 25.0                 # ns
ecalBarrelInputName = "ECalBarrelModuleThetaMerged"  # name of the ECal barrel readout in input file
ecalBarrelSignalShapePath = "SignalPulseShapes.root"  # path to the root file with the signal pulse shapes



# ECAL barrel parameters for digitisation
ecalBarrelSamplingFraction = [0.3800493723322256] * 1 + [0.13494147915064658] * 1 + [0.142866851721152] * 1 + [0.14839315921940666] * 1 + [0.15298362570665006] * 1 + [0.15709704561942747] * 1 + [0.16063717490147533] * 1 + [0.1641723795419055] * 1 + [0.16845490287689746] * 1 + [0.17111520115997653] * 1 + [0.1730605163148862] * 1
ecalBarrelUpstreamParameters = [[0.028158491043365624, -1.564259408365951, -76.52312805346982, 0.7442903558010191, -34.894692961350195, -74.19340877431723]]
ecalBarrelDownstreamParameters = [[0.00010587711361028165, 0.0052371999097777355, 0.69906696456064, -0.9348243433360095, -0.0364714212117143, 8.360401126995626]]

ecalBarrelThetaWeights = [-1, 3.0, 3.0, 3.0, 4.25, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0]  # to be recalculated for V03, separately for topo and calo clusters...

from Configurables import k4DataSvc, PodioInput
podioevent = k4DataSvc('EventDataSvc')
podioevent.input = inputfile
input_reader = PodioInput('InputReader')


# Detector geometry
geoservice = GeoSvc("GeoSvc")
# if K4GEO is empty, this should use relative path to working directory
path_to_detector = os.environ.get("K4GEO", "")
print(path_to_detector)
detectors_to_use = [
    'FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml'
]

print("Using the following detectors:", [os.path.join(path_to_detector, _det) for _det in detectors_to_use])

# prefix all xmls with path_to_detector
geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]
geoservice.OutputLevel = INFO

# - Calibrate cells to EM scale
from Configurables import CalibrateInLayersTool
# This algorithm loops over Hits, which looks to be a vector of <cell ID, energy deposit> pairs. It loops over these pairs, figures out the layer by the ID and then applies the sampling fraction (by dividing by the sampling fraction for that layer). This algorithm overwrites the energy deposit in the cell with the calibrated energy deposit!!

calibEcalBarrel = CalibrateInLayersTool("CalibrateECalBarrel",
                        samplingFraction=ecalBarrelSamplingFraction, # hard coded above
                        readoutName=ecalBarrelInputName, # Can be found from the geant xml files under <readouts> inECalBarrel_thetamodulemerged.xml
                        layerFieldName="layer") # I assume this is linked to the GEANT files as well somehow....

from Configurables import CreatePositionedCaloCells
"""
This algorithm creates and stores digits (Energy of hit * pulse shape) as a timeseries object for each cell ID. The digitizer tool and sampling interval are given as inputs to the algorithm.
"""
from Configurables import CaloDigitizer
ECALBarrelDigits = CaloDigitizer("CaloDigitizer",
                        hits=ecalBarrelInputName,
                        digitizerTool = DigitzerReaderEcalBarrel,
                        SamplingInterval = SamplingInterval,
                        OutputLevel=INFO)


# Output
from Configurables import PodioOutput
out = PodioOutput("out",
                  OutputLevel=INFO)
out.filename = "ALLEGRO_sim_digi_reco_test.root"

out.outputCommands = ["keep *"
                      ]

# drop lumi, vertex, DCH, Muons (unless want to keep for event display)
out.outputCommands.append("drop Lumi*")
out.outputCommands.append("drop Vertex*")
out.outputCommands.append("drop DriftChamber_simHits*")
out.outputCommands.append("drop MuonTagger*")
out.outputCommands.append("drop *SiWr*")

#Get rid of HCAL and ECAL (only barrel is kept)
out.outputCommands.append("drop *ECalEndcap*")
out.outputCommands.append("drop *HCal*")

# CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
out.AuditExecute = True

# Configure list of external services
ExtSvc = [geoservice, podioevent, audsvc]
ECALBarrelDigits.AuditExecute = True

# Setup alg sequence
TopAlg = [
    input_reader,
    ECALBarrelDigits,
    out,
]

# ApplicationMgr
ApplicationMgr(TopAlg=TopAlg,
               EvtSel='NONE',
               EvtMax=Nevts,
               # order is important, as GeoSvc is needed by G4SimSvc
               ExtSvc=ExtSvc,
               StopOnSignal=True,
               OutputLevel=INFO
               )
