from Configurables import ApplicationMgr
from Configurables import AuditorSvc, ChronoAuditor
from Configurables import PodioOutput
from Configurables import CorrectCaloClusters
from Configurables import CreateCaloClustersSlidingWindow
from Configurables import CaloTowerTool
from Configurables import CreateEmptyCaloCellsCollection
from Configurables import CreateCaloCellPositionsFCCee
from Configurables import CellPositionsECalBarrelModuleThetaSegTool
from Configurables import RedoSegmentation
from Configurables import CreateCaloCells
from Configurables import CalibrateCaloHitsTool
from Configurables import CalibrateInLayersTool
from Configurables import SimG4Alg
from Configurables import SimG4PrimariesFromEdmTool
from Configurables import SimG4SaveCalHits
from Configurables import SimG4ConstantMagneticFieldTool
from Configurables import SimG4Svc
from Configurables import SimG4FullSimActions
from Configurables import SimG4SaveParticleHistory
from Configurables import GeoSvc
from Configurables import HepMCToEDMConverter
from Configurables import GenAlg
from Configurables import FCCDataSvc
from Configurables import CaloTopoClusterInputTool
from Configurables import TopoCaloNeighbours
from Configurables import TopoCaloNoisyCells
from Configurables import CaloTopoClusterFCCee

from Gaudi.Configuration import *
import os

from GaudiKernel.SystemOfUnits import GeV, tesla

use_pythia = False
addNoise = False
dumpGDML = False
runHCal = False

# Input for simulations (momentum is expected in GeV!)
# Parameters for the particle gun simulations, dummy if use_pythia is set
# to True
# theta from 80 to 100 degrees corresponds to -0.17 < eta < 0.17
# reminder: cell granularity in theta = 0.5625 degrees
# (in strips: 0.5625/4=0.14)

# Nevts = 20000
Nevts = 100
# Nevts = 1
# Nevts=1000
# momentum = 100 # in GeV
# momentum = 50  # in GeV
momentum = 10 # in GeV
_pi = 3.14159
thetaMin = 40 # degrees
thetaMax = 140 # degrees
# thetaMin = 89
# thetaMax = 91
# thetaMin = 90 # degrees
# thetaMax = 90 # degrees
# phiMin = _pi/2.
# phiMax = _pi/2.
phiMin = 0
phiMax = 2 * _pi

# particle type: 11 electron, 13 muon, 22 photon, 111 pi0, 211 pi+
pdgCode = 11
# pdgCode = 22
# pdgCode = 111
# pdgCode = 211

# Set to true if history from Geant4 decays is needed (e.g. to get the
# photons from pi0)
saveG4Hist = False
if (pdgCode == 111):
    saveG4Hist = True

magneticField = False


podioevent = FCCDataSvc("EventDataSvc")

# Particle gun setup

genAlg = GenAlg()
if use_pythia:
    from Configurables import PythiaInterface
    pythia8gentool = PythiaInterface()
    pythia8gentool.pythiacard = os.path.join(
        os.environ.get('PWD', ''),
        "MCGeneration/ee_Zgamma_inclusive.cmd"
    )
    # pythia8gentool.pythiacard = "MCGeneration/ee_Z_ee.cmd"
    pythia8gentool.printPythiaStatistics = False
    pythia8gentool.pythiaExtraSettings = [""]
    genAlg.SignalProvider = pythia8gentool
    # to smear the primary vertex position:
    # from Configurables import GaussSmearVertex
    # smeartool = GaussSmearVertex()
    # smeartool.xVertexSigma =   0.5*units.mm
    # smeartool.yVertexSigma =   0.5*units.mm
    # smeartool.zVertexSigma =  40.0*units.mm
    # smeartool.tVertexSigma = 180.0*units.picosecond
    # genAlg.VertexSmearingTool = smeartool
else:
    from Configurables import MomentumRangeParticleGun
    pgun = MomentumRangeParticleGun("ParticleGun")
    pgun.PdgCodes = [pdgCode]
    pgun.MomentumMin = momentum * GeV
    pgun.MomentumMax = momentum * GeV
    pgun.PhiMin = phiMin
    pgun.PhiMax = phiMax
    pgun.ThetaMin = thetaMin * _pi / 180.
    pgun.ThetaMax = thetaMax * _pi / 180.
    genAlg.SignalProvider = pgun

genAlg.hepmc.Path = "hepmc"

hepmc_converter = HepMCToEDMConverter()
hepmc_converter.hepmc.Path = "hepmc"
genParticlesOutputName = "genParticles"
hepmc_converter.GenParticles.Path = genParticlesOutputName
hepmc_converter.hepmcStatusList = []
hepmc_converter.OutputLevel = INFO

# Simulation setup
# Detector geometry
geoservice = GeoSvc("GeoSvc")
# if K4GEO is empty, this should use relative path to working directory
path_to_detector = os.environ.get("K4GEO", "")
print(path_to_detector)
detectors_to_use = [
    #    'Detector/DetFCCeeIDEA-LAr/compact/FCCee_DectMaster_thetamodulemerged.xml',
    'FCCee/ALLEGRO/compact/ALLEGRO_o1_v02/ALLEGRO_o1_v02.xml'
]
# prefix all xmls with path_to_detector
geoservice.detectors = [
    os.path.join(path_to_detector, _det) for _det in detectors_to_use
]
geoservice.OutputLevel = INFO

# Geant4 service
# Configures the Geant simulation: geometry, physics list and user actions
actions = SimG4FullSimActions()

if saveG4Hist:
    actions.enableHistory = True
    actions.energyCut = 1.0 * GeV
    saveHistTool = SimG4SaveParticleHistory("saveHistory")

geantservice = SimG4Svc(
    "SimG4Svc",
    detector='SimG4DD4hepDetector',
    physicslist="SimG4FtfpBert",
    actions=actions
)

# from Configurables import GeoToGdmlDumpSvc
if dumpGDML:
    from Configurables import GeoToGdmlDumpSvc
    gdmldumpservice = GeoToGdmlDumpSvc("GeoToGdmlDumpSvc")

# Fixed seed to have reproducible results, change it for each job if you
# split one production into several jobs
# Mind that if you leave Gaudi handle random seed and some job start within
# the same second (very likely) you will have duplicates
geantservice.randomNumbersFromGaudi = False
geantservice.seedValue = 4242

# Range cut
geantservice.g4PreInitCommands += ["/run/setCut 0.1 mm"]

# Magnetic field
if magneticField == 1:
    field = SimG4ConstantMagneticFieldTool(
        "SimG4ConstantMagneticFieldTool",
        FieldComponentZ=-2*tesla,
        FieldOn=True,
        IntegratorStepper="ClassicalRK4"
    )
else:
    field = SimG4ConstantMagneticFieldTool(
        "SimG4ConstantMagneticFieldTool",
        FieldOn=False
    )

# Geant4 algorithm
# Translates EDM to G4Event, passes the event to G4, writes out outputs
# via tools and a tool that saves the calorimeter hits

# Detector readouts
# ECAL
ecalBarrelReadoutName = "ECalBarrelModuleThetaMerged"
ecalBarrelReadoutName2 = "ECalBarrelModuleThetaMerged2"
ecalEndcapReadoutName = "ECalEndcapPhiEta"
# HCAL
if runHCal:
    hcalBarrelReadoutName = "HCalBarrelReadout"
    hcalEndcapReadoutName = "HCalEndcapReadout"
else:
    hcalBarrelReadoutName = ""
    hcalEndcapReadoutName = ""

# Configure saving of calorimeter hits
ecalBarrelHitsName = "ECalBarrelPositionedHits"
saveECalBarrelTool = SimG4SaveCalHits(
    "saveECalBarrelHits",
    readoutName=ecalBarrelReadoutName,
    OutputLevel=INFO
)
saveECalBarrelTool.CaloHits.Path = ecalBarrelHitsName

saveECalEndcapTool = SimG4SaveCalHits(
    "saveECalEndcapHits",
    readoutName=ecalEndcapReadoutName
)
saveECalEndcapTool.CaloHits.Path = "ECalEndcapHits"

if runHCal:
    hcalBarrelHitsName = "HCalBarrelPositionedHits"
    saveHCalTool = SimG4SaveCalHits(
        "saveHCalBarrelHits",
        readoutName=hcalBarrelReadoutName
    )
    saveHCalTool.CaloHits.Path = hcalBarrelHitsName

    # saveHCalEndcapTool = SimG4SaveCalHits(
    #    "saveHCalEndcapHits",
    #    readoutName = hcalEndcapReadoutName
    # )
    # saveHCalEndcapTool.CaloHits.Path = "HCalEndcapHits"

# next, create the G4 algorithm, giving the list of names of tools ("XX/YY")
particle_converter = SimG4PrimariesFromEdmTool("EdmConverter")
particle_converter.GenParticles.Path = genParticlesOutputName

outputTools = [
    saveECalBarrelTool,
    saveECalEndcapTool
]
if runHCal:
    outputTools += [
        saveHCalTool,
        # saveHCalEndcapTool
    ]

if saveG4Hist:
    outputTools += [saveHistTool]

geantsim = SimG4Alg("SimG4Alg",
                    outputs=outputTools,
                    eventProvider=particle_converter,
                    OutputLevel=INFO)

# Digitization (Merging hits into cells, EM scale calibration)
# EM scale calibration (sampling fraction)
calibEcalBarrel = CalibrateInLayersTool("CalibrateECalBarrel",
                                        samplingFraction=[0.36599110182660616] * 1 + [0.1366222373338866] * 1 + [0.1452035173747207] * 1 + [0.1504319190969367] * 1 + [0.15512713637727382] * 1 + [0.1592916726494782] * 1 + [
                                            0.16363478857307595] * 1 + [0.1674697333180323] * 1 + [0.16998205747422343] * 1 + [0.1739146363733975] * 1 + [0.17624609543603845] * 1 + [0.1768613530850488] * 1,
                                        readoutName=ecalBarrelReadoutName,
                                        layerFieldName="layer")

calibEcalEndcap = CalibrateCaloHitsTool(
    "CalibrateECalEndcap", invSamplingFraction="4.27")
if runHCal:
    calibHcells = CalibrateCaloHitsTool(
        "CalibrateHCal", invSamplingFraction="31.4")
    calibHcalEndcap = CalibrateCaloHitsTool(
        "CalibrateHCalEndcap", invSamplingFraction="31.7")

# Create cells in ECal barrel
# 1. step - merge hits into cells with theta and module segmentation
# (module is a 'physical' cell i.e. lead + LAr + PCB + LAr +lead)
# 2. step - rewrite the cellId using the merged theta-module segmentation
# (merging several modules and severla theta readout cells).
# Add noise at this step if you derived the noise already assuming merged cells

# Step 1: merge hits into cells according to initial segmentation
ecalBarrelCellsName = "ECalBarrelCells"
createEcalBarrelCells = CreateCaloCells("CreateECalBarrelCells",
                                        doCellCalibration=True,
                                        calibTool=calibEcalBarrel,
                                        addCellNoise=False,
                                        filterCellNoise=False,
                                        addPosition=True,
                                        OutputLevel=INFO,
                                        hits=ecalBarrelHitsName,
                                        cells=ecalBarrelCellsName)

# Step 2a: compute new cellID of cells based on new readout
# (merged module-theta segmentation with variable merging vs layer)
resegmentEcalBarrel = RedoSegmentation("ReSegmentationEcal",
                                       # old bitfield (readout)
                                       oldReadoutName=ecalBarrelReadoutName,
                                       # specify which fields are going to be altered (deleted/rewritten)
                                       oldSegmentationIds=["module", "theta"],
                                       # new bitfield (readout), with new segmentation (merged modules and theta cells)
                                       newReadoutName=ecalBarrelReadoutName2,
                                       OutputLevel=INFO,
                                       debugPrint=200,
                                       inhits=ecalBarrelCellsName,
                                       outhits="ECalBarrelCellsMerged")

# Step 2b: merge new cells with same cellID together
# do not apply cell calibration again since cells were already
# calibrated in Step 1
ecalBarrelCellsName2 = "ECalBarrelCells2"
createEcalBarrelCells2 = CreateCaloCells("CreateECalBarrelCells2",
                                         doCellCalibration=False,
                                         addCellNoise=False,
                                         filterCellNoise=False,
                                         OutputLevel=INFO,
                                         hits="ECalBarrelCellsMerged",
                                         cells=ecalBarrelCellsName2)

# Add to Ecal barrel cells the position information
# (good for physics, all coordinates set properly)

cellPositionEcalBarrelTool = CellPositionsECalBarrelModuleThetaSegTool(
    "CellPositionsECalBarrel",
    readoutName=ecalBarrelReadoutName,
    OutputLevel=INFO
)
ecalBarrelPositionedCellsName = "ECalBarrelPositionedCells"
createEcalBarrelPositionedCells = CreateCaloCellPositionsFCCee(
    "CreateECalBarrelPositionedCells",
    OutputLevel=INFO
)
createEcalBarrelPositionedCells.positionsTool = cellPositionEcalBarrelTool
createEcalBarrelPositionedCells.hits.Path = ecalBarrelCellsName
createEcalBarrelPositionedCells.positionedHits.Path = ecalBarrelPositionedCellsName

cellPositionEcalBarrelTool2 = CellPositionsECalBarrelModuleThetaSegTool(
    "CellPositionsECalBarrel2",
    readoutName=ecalBarrelReadoutName2,
    OutputLevel=INFO
)
createEcalBarrelPositionedCells2 = CreateCaloCellPositionsFCCee(
    "CreateECalBarrelPositionedCells2",
    OutputLevel=INFO
)
createEcalBarrelPositionedCells2.positionsTool = cellPositionEcalBarrelTool2
createEcalBarrelPositionedCells2.hits.Path = ecalBarrelCellsName2
createEcalBarrelPositionedCells2.positionedHits.Path = "ECalBarrelPositionedCells2"


# Create cells in ECal endcap
createEcalEndcapCells = CreateCaloCells("CreateEcalEndcapCaloCells",
                                        doCellCalibration=True,
                                        calibTool=calibEcalEndcap,
                                        addCellNoise=False,
                                        filterCellNoise=False,
                                        OutputLevel=INFO)
createEcalEndcapCells.hits.Path = "ECalEndcapHits"
createEcalEndcapCells.cells.Path = "ECalEndcapCells"

if runHCal:
    # Create cells in HCal
    # 1. step - merge hits into cells with the default readout
    hcalBarrelCellsName = "HCalBarrelCells"
    createHcalBarrelCells = CreateCaloCells("CreateHCalBarrelCells",
                                            doCellCalibration=True,
                                            calibTool=calibHcells,
                                            addCellNoise=False,
                                            filterCellNoise=False,
                                            addPosition=True,
                                            hits=hcalBarrelHitsName,
                                            cells=hcalBarrelCellsName,
                                            OutputLevel=INFO)

    # createHcalEndcapCells = CreateCaloCells("CreateHcalEndcapCaloCells",
    #                                    doCellCalibration=True,
    #                                    calibTool=calibHcalEndcap,
    #                                    addCellNoise=False,
    #                                    filterCellNoise=False,
    #                                    OutputLevel=INFO)
    # createHcalEndcapCells.hits.Path="HCalEndcapHits"
    # createHcalEndcapCells.cells.Path="HCalEndcapCells"

    from Configurables import CellPositionsHCalBarrelPhiThetaSegTool
    cellPositionHcalBarrelTool = CellPositionsHCalBarrelPhiThetaSegTool(
        "CellPositionsHCalBarrel",
        readoutName=hcalBarrelReadoutName,
        OutputLevel=INFO
    )
    hcalBarrelPositionedCellsName = "HCalBarrelPositionedCells"    
    createHcalBarrelPositionedCells = CreateCaloCellPositionsFCCee(
        "CreateHcalBarrelPositionedCells",
        OutputLevel=INFO
    )
    createHcalBarrelPositionedCells.positionsTool = cellPositionHcalBarrelTool
    createHcalBarrelPositionedCells.hits.Path = hcalBarrelCellsName
    createHcalBarrelPositionedCells.positionedHits.Path = hcalBarrelPositionedCellsName
else:
    hcalBarrelCellsName = "emptyCaloCells"
    hcalBarrelPositionedCellsName = "emptyCaloCells"
    cellPositionHcalBarrelTool = None
    
# Empty cells for parts of calorimeter not implemented yet
createemptycells = CreateEmptyCaloCellsCollection("CreateEmptyCaloCells")
createemptycells.cells.Path = "emptyCaloCells"

# Produce sliding window clusters
towers = CaloTowerTool("towers",
                       deltaEtaTower=0.01, deltaPhiTower=2*_pi/768,
                       radiusForPosition=2160 + 40 / 2.0,
                       # ecalBarrelReadoutName = ecalBarrelReadoutNamePhiEta,
                       ecalBarrelReadoutName=ecalBarrelReadoutName,
                       ecalEndcapReadoutName=ecalEndcapReadoutName,
                       ecalFwdReadoutName="",
                       hcalBarrelReadoutName=hcalBarrelReadoutName,
                       hcalExtBarrelReadoutName="",
                       hcalEndcapReadoutName="",
                       hcalFwdReadoutName="",
                       OutputLevel=INFO)
towers.ecalBarrelCells.Path = ecalBarrelCellsName
towers.ecalEndcapCells.Path = "ECalEndcapCells"
towers.ecalFwdCells.Path = "emptyCaloCells"

towers.hcalBarrelCells.Path = hcalBarrelCellsName
towers.hcalExtBarrelCells.Path = "emptyCaloCells"
towers.hcalEndcapCells.Path = "emptyCaloCells"
towers.hcalFwdCells.Path = "emptyCaloCells"

# Cluster variables
windE = 9
windP = 17
posE = 5
posP = 11
dupE = 7
dupP = 13
finE = 9
finP = 17
# Minimal energy to create a cluster in GeV (FCC-ee detectors have to reconstruct low energy particles)
threshold = 0.040

createClusters = CreateCaloClustersSlidingWindow("CreateClusters",
                                                 towerTool=towers,
                                                 nEtaWindow=windE, nPhiWindow=windP,
                                                 nEtaPosition=posE, nPhiPosition=posP,
                                                 nEtaDuplicates=dupE, nPhiDuplicates=dupP,
                                                 nEtaFinal=finE, nPhiFinal=finP,
                                                 energyThreshold=threshold,
                                                 energySharingCorrection=False,
                                                 attachCells=True,
                                                 OutputLevel=INFO
                                                 )
createClusters.clusters.Path = "CaloClusters"
createClusters.clusterCells.Path = "CaloClusterCells"

createEcalBarrelPositionedCaloClusterCells = CreateCaloCellPositionsFCCee(
    "ECalBarrelPositionedCaloClusterCells",
    OutputLevel=INFO
)
createEcalBarrelPositionedCaloClusterCells.positionsTool = cellPositionEcalBarrelTool
createEcalBarrelPositionedCaloClusterCells.hits.Path = "CaloClusterCells"
createEcalBarrelPositionedCaloClusterCells.positionedHits.Path = "PositionedCaloClusterCells"

correctCaloClusters = CorrectCaloClusters("correctCaloClusters",
                                          inClusters=createClusters.clusters.Path,
                                          outClusters="Corrected"+createClusters.clusters.Path,
                                          numLayers=[0],
                                          firstLayerIDs=[0],
                                          lastLayerIDs=[-1],
                                          # readoutNames = [ecalBarrelReadoutNamePhiEta],
                                          readoutNames=[ecalBarrelReadoutName],
                                          # upstreamParameters = [[0.02729094887360858, -1.378665489864182, -68.40424543618059, 3.6930827214130053, -5528.714729126099, -1630.7911298009794]],
                                          upstreamParameters=[
                                              [0.02729094887360858, -1.378665489864182, -68.40424543618059, 3.6930827214130053, -5528.714729126099, -1630.7911298009794]],
                                          upstreamFormulas=[
                                              ['[0]+[1]/(x-[2])', '[0]+[1]/(x-[2])']],
                                          # downstreamParameters = [[-0.0032351643028483354, 0.006597484738888312, 0.8972024981692965, -1.0207168610322181, 0.017878133854084398, 9.108099243443101]],
                                          downstreamParameters=[
                                              [-0.0032351643028483354, 0.006597484738888312, 0.8972024981692965, -1.0207168610322181, 0.017878133854084398, 9.108099243443101]],
                                          downstreamFormulas=[
                                              ['[0]+[1]*x', '[0]+[1]/sqrt(x)', '[0]+[1]/x']],
                                          OutputLevel=INFO
                                          )

# TOPO CLUSTERS PRODUCTION
createTopoInput = CaloTopoClusterInputTool("CreateTopoInput",
                                           ecalBarrelReadoutName=ecalBarrelReadoutName,
                                           ecalEndcapReadoutName="",
                                           ecalFwdReadoutName="",
                                           hcalBarrelReadoutName=hcalBarrelReadoutName,
                                           hcalExtBarrelReadoutName="",
                                           hcalEndcapReadoutName="",
                                           hcalFwdReadoutName="",
                                           OutputLevel=INFO)

createTopoInput.ecalBarrelCells.Path = ecalBarrelPositionedCellsName
# createTopoInput.ecalBarrelCells.Path = "ECalBarrelPositionedCells2"
createTopoInput.ecalEndcapCells.Path = "emptyCaloCells"
createTopoInput.ecalFwdCells.Path = "emptyCaloCells"
createTopoInput.hcalBarrelCells.Path = hcalBarrelPositionedCellsName
createTopoInput.hcalExtBarrelCells.Path = "emptyCaloCells"
createTopoInput.hcalEndcapCells.Path = "emptyCaloCells"
createTopoInput.hcalFwdCells.Path = "emptyCaloCells"
cellPositionHcalBarrelNoSegTool = None
cellPositionHcalExtBarrelTool = None

readNeighboursMap = TopoCaloNeighbours("ReadNeighboursMap",
                                       fileName=os.environ['FCCBASEDIR'] +
                                       "/LAr_scripts/data/neighbours_map_barrel_thetamodulemerged.root",
                                       #"/LAr_scripts/data/neighbours_map_HCalBarrel.root",
                                       OutputLevel=INFO)

# Noise levels per cell
readNoisyCellsMap = TopoCaloNoisyCells("ReadNoisyCellsMap",
                                       fileName=os.environ['FCCBASEDIR'] +
                                       "/LAr_scripts/data/cellNoise_map_electronicsNoiseLevel_thetamodulemerged.root",
                                       OutputLevel=INFO)

createTopoClusters = CaloTopoClusterFCCee("CreateTopoClusters",
                                          TopoClusterInput=createTopoInput,
                                          # expects neighbours map from cellid->vec < neighbourIds >
                                          neigboursTool=readNeighboursMap,
                                          # tool to get noise level per cellid
                                          noiseTool=readNoisyCellsMap,
                                          # cell positions tools for all sub - systems
                                          positionsECalBarrelTool=cellPositionEcalBarrelTool,
                                          positionsHCalBarrelTool=cellPositionHcalBarrelTool,
                                          positionsHCalBarrelNoSegTool=cellPositionHcalBarrelNoSegTool,
                                          positionsHCalExtBarrelTool=cellPositionHcalExtBarrelTool,
                                          # positionsHCalExtBarrelTool = HCalExtBcells,
                                          # positionsEMECTool = EMECcells,
                                          # positionsHECTool = HECcells,
                                          # positionsEMFwdTool = ECalFwdcells,
                                          # positionsHFwdTool = HCalFwdcells,
                                          noSegmentationHCal=False,
                                          seedSigma=4,
                                          neighbourSigma=2,
                                          lastNeighbourSigma=0,
                                          OutputLevel=INFO)
createTopoClusters.clusters.Path = "CaloTopoClusters"
createTopoClusters.clusterCells.Path = "CaloTopoClusterCells"

createEcalBarrelPositionedCaloTopoClusterCells = CreateCaloCellPositionsFCCee(
    "ECalBarrelPositionedCaloTopoClusterCells",
    OutputLevel=INFO
)
# createEcalBarrelPositionedCaloTopoClusterCells.positionsTool = cellPositionEcalBarrelTool2
createEcalBarrelPositionedCaloTopoClusterCells.positionsTool = cellPositionEcalBarrelTool
createEcalBarrelPositionedCaloTopoClusterCells.hits.Path = "CaloTopoClusterCells"
createEcalBarrelPositionedCaloTopoClusterCells.positionedHits.Path = "PositionedCaloTopoClusterCells"

correctCaloTopoClusters = CorrectCaloClusters(
    "correctCaloTopoClusters",
    inClusters=createTopoClusters.clusters.Path,
    outClusters="Corrected"+createTopoClusters.clusters.Path,
    numLayers=[0],
    firstLayerIDs=[0],
    lastLayerIDs=[-1],
    # readoutNames = [ecalBarrelReadoutNamePhiEta],
    readoutNames=[ecalBarrelReadoutName],
    # upstreamParameters = [[0.02729094887360858, -1.378665489864182, -68.40424543618059, 3.6930827214130053, -5528.714729126099, -1630.7911298009794]],
    upstreamParameters=[[0.02729094887360858, -1.378665489864182, -68.40424543618059,
                         3.6930827214130053, -5528.714729126099, -1630.7911298009794]],
    upstreamFormulas=[['[0]+[1]/(x-[2])', '[0]+[1]/(x-[2])']],
    # downstreamParameters = [[-0.0032351643028483354, 0.006597484738888312, 0.8972024981692965, -1.0207168610322181, 0.017878133854084398, 9.108099243443101]],
    downstreamParameters=[[-0.0032351643028483354, 0.006597484738888312,
                           0.8972024981692965, -1.0207168610322181, 0.017878133854084398, 9.108099243443101]],
    downstreamFormulas=[['[0]+[1]*x', '[0]+[1]/sqrt(x)', '[0]+[1]/x']],
    OutputLevel=INFO
)

# Output
out = PodioOutput("out",
                  OutputLevel=INFO)

# out.outputCommands = ["keep *"]
# out.outputCommands = ["keep *", "drop ECalBarrelHits", "drop HCal*", "drop ECalBarrelCellsStep*", "drop ECalBarrelPositionedHits", "drop emptyCaloCells", "drop CaloClusterCells"]
# out.outputCommands = ["keep *", "drop ECalBarrelHits", "drop HCal*", "drop ECalBarrelCellsStep*", "drop ECalBarrelPositionedHits", "drop emptyCaloCells", "drop CaloClusterCells", "drop %s"%ecalBarrelCellsName, "drop %s"%createEcalBarrelPositionedCells.positionedHits.Path]
# out.outputCommands = ["keep *", "drop ECalBarrelHits", "drop HCal*", "drop *ells*", "drop ECalBarrelPositionedHits", "drop emptyCaloCells"]
# out.outputCommands = ["keep *", "drop HCal*", "drop ECalBarrel*", "drop emptyCaloCells"]
if runHCal:
    out.outputCommands = ["keep *", "drop emptyCaloCells"]
else:
    out.outputCommands = ["keep *", "drop HCal*", "drop emptyCaloCells"]

# out.filename = "root/output_fullCalo_SimAndDigi_withTopoCluster_MagneticField_"+str(magneticField)+"_pMin_"+str(momentum*1000)+"_MeV"+"_ThetaMinMax_"+str(thetaMin)+"_"+str(thetaMax)+"_pdgId_"+str(pdgCode)+"_pythia"+str(use_pythia)+"_Noise"+str(addNoise)+".root"
out.filename = "./root_merge/output_evts_"+str(Nevts)+"_pdg_"+str(pdgCode)+"_"+str(momentum)+"_GeV"+"_ThetaMinMax_"+str(thetaMin)+"_"+str(
    thetaMax)+"_PhiMinMax_"+str(phiMin)+"_"+str(phiMax)+"_MagneticField_"+str(magneticField)+"_Noise"+str(addNoise)+".root"

# CPU information
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
genAlg.AuditExecute = True
hepmc_converter.AuditExecute = True
geantsim.AuditExecute = True
createEcalBarrelCells.AuditExecute = True
createEcalBarrelPositionedCells.AuditExecute = True
if runHCal:
    createHcalBarrelCells.AuditExecute = True
createTopoClusters.AuditExecute = True
out.AuditExecute = True

ExtSvc = [geoservice, podioevent, geantservice, audsvc]
if dumpGDML:
    ExtSvc += [gdmldumpservice]

TopAlg = [
    genAlg,
    hepmc_converter,
    geantsim,
    createEcalBarrelCells,
    createEcalBarrelPositionedCells,
    resegmentEcalBarrel,
    createEcalBarrelCells2,
    createEcalBarrelPositionedCells2,
    createEcalEndcapCells
]
if runHCal:
    TopAlg += [
        createHcalBarrelCells,
        createHcalBarrelPositionedCells,
        # createHcalEndcapCells
    ]
TopAlg += [
    createemptycells,
    # createClusters,
    # createEcalBarrelPositionedCaloClusterCells,
    # correctCaloClusters,
    createTopoClusters,
    createEcalBarrelPositionedCaloTopoClusterCells,
    correctCaloTopoClusters,
    out
]

ApplicationMgr(
    TopAlg=TopAlg,
    EvtSel='NONE',
    EvtMax=Nevts,
    ExtSvc=ExtSvc,
    StopOnSignal=True,
)
