from Configurables import ApplicationMgr
from Configurables import AuditorSvc, ChronoAuditor
from Configurables import PodioOutput
from Configurables import CaloTowerToolFCCee
from Configurables import CreateCaloClustersSlidingWindowFCCee
from Configurables import CorrectCaloClusters
from Configurables import CalibrateCaloClusters
from Configurables import AugmentClustersFCCee
from Configurables import CreateEmptyCaloCellsCollection
from Configurables import CreateCaloCellPositionsFCCee
from Configurables import CellPositionsECalBarrelModuleThetaSegTool
from Configurables import CellPositionsECalEndcapTurbineSegTool
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
from Configurables import RewriteBitfield
from Configurables import ReadCaloCrosstalkMap
from Gaudi.Configuration import INFO

import os

from GaudiKernel.SystemOfUnits import GeV, tesla, mm
from GaudiKernel.PhysicalConstants import pi, halfpi, twopi
from math import cos, sin, tan

use_pythia = False
addNoise = False
dumpGDML = False
runHCal = False
# for big productions, save significant space removing hits and cells
# however, hits and cluster cells might be wanted for small productions for detailed event displays
# also, cluster cells are needed for the MVA training
saveHits = False
saveCells = False
saveClusterCells = False

# clustering
doClustering = True
# NOTE: since topoclustering requires root files with noise and neighbours that are
# not in the release, topoclusters are disabled
# cluster energy corrections
# simple parametrisations of up/downstream losses
applyUpDownstreamCorrections = True
# calculate cluster energy and barycenter per layer and save it as extra parameters
addShapeParameters = True

# Input for simulations (momentum is expected in GeV!)
# Parameters for the particle gun simulations, dummy if use_pythia is set
# to True
# theta from 80 to 100 degrees corresponds to -0.17 < eta < 0.17
# reminder: cell granularity in theta = 0.5625 degrees
# (in strips: 0.5625/4=0.14)

# Nevts = 20000
Nevts = 2
# Nevts = 1
# Nevts=1000

# particle momentum and direction
# momentum = 100  # in GeV
momentum = 10  # in GeV
# momentum = 10  # in GeV
thetaMin = 45  # degrees
thetaMax = 135  # degrees
# thetaMin = 89
# thetaMax = 91
# thetaMin = 90  # degrees
# thetaMax = 90  # degrees
# phiMin = halfpi
# phiMax = halfpi
phiMin = 0
phiMax = twopi

# particle origin
# origR = 1000.0*mm
origR = 0.0 * mm
origTheta = halfpi
origPhi = 0.0

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
else:
    from Configurables import MomentumRangeParticleGun
    pgun = MomentumRangeParticleGun("ParticleGun")
    pgun.PdgCodes = [pdgCode]
    pgun.MomentumMin = momentum * GeV
    pgun.MomentumMax = momentum * GeV
    pgun.PhiMin = phiMin
    pgun.PhiMax = phiMax
    pgun.ThetaMin = thetaMin * pi / 180.
    pgun.ThetaMax = thetaMax * pi / 180.
    genAlg.SignalProvider = pgun

genAlg.hepmc.Path = "hepmc"

# smear/shift vertex
if origR > 0.0:
    origX = origR * cos(origPhi)
    origY = origR * sin(origPhi)
    origZ = origR / tan(origTheta)
    print("Particle gun will be moved to %f %f %f" % (origX, origY, origZ))
    from Configurables import GaussSmearVertex
    vertexSmearAndShiftTool = GaussSmearVertex()
    vertexSmearAndShiftTool.xVertexSigma = 0.
    vertexSmearAndShiftTool.yVertexSigma = 0.
    vertexSmearAndShiftTool.tVertexSigma = 0.
    vertexSmearAndShiftTool.xVertexMean = origX
    vertexSmearAndShiftTool.yVertexMean = origY
    vertexSmearAndShiftTool.zVertexMean = origZ
    vertexSmearAndShiftTool.tVertexMean = 0.
    genAlg.VertexSmearingTool = vertexSmearAndShiftTool

# hepMC -> EDM converter
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
    'FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml'
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
        FieldComponentZ=-2 * tesla,
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
ecalEndcapReadoutName = "ECalEndcapTurbine"
# HCAL
if runHCal:
    hcalBarrelReadoutName = "HCalBarrelReadout"
    hcalBarrelReadoutName2 = "BarHCal_Readout_phitheta"
    hcalEndcapReadoutName = "HCalEndcapReadout"
else:
    hcalBarrelReadoutName = ""
    hcalBarrelReadoutName2 = ""
    hcalEndcapReadoutName = ""

# Configure saving of positioned calorimeter hits
ecalBarrelHitsName = "ECalBarrelPositionedHits"
saveECalBarrelTool = SimG4SaveCalHits(
    "saveECalBarrelHits",
    readoutName=ecalBarrelReadoutName,
    OutputLevel=INFO
)
saveECalBarrelTool.CaloHits.Path = ecalBarrelHitsName

ecalEndcapHitsName = "ECalEndcapPositionedHits"
saveECalEndcapTool = SimG4SaveCalHits(
    "saveECalEndcapHits",
    readoutName=ecalEndcapReadoutName
)
saveECalEndcapTool.CaloHits.Path = ecalEndcapHitsName

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
                                        samplingFraction=[0.3775596654349802] * 1 + [0.13400227700041234] * 1 + [0.14390509963164044] * 1 + [0.14998482026270935] * 1 + [0.15457673722531148] * 1 + [0.15928098152159675] * 1 + [0.1635367867767212] * 1 + [0.16801070646031507] * 1 + [0.1713409944779989] * 1 + [0.17580195406064622] * 1 + [0.17966699467772812] * 1,
                                        readoutName=ecalBarrelReadoutName,
                                        layerFieldName="layer")

#calibEcalEndcap = CalibrateCaloHitsTool(
#    "CalibrateECalEndcap", invSamplingFraction="4.27")
calibEcalEndcap = CalibrateInLayersTool("CalibrateECalEndcap",
                                        samplingFraction = [0.16419] * 1 + [0.192898] * 1 + [0.18783] * 1 + [0.193203] * 1 + [0.193928] * 1 + [0.192286] * 1 + [0.199959] * 1 + [0.200153] * 1 + [0.212635] * 1 + [0.180345] * 1 + [0.18488] * 1 + [0.194762] * 1 + [0.197775] * 1 + [0.200504] * 1 + [0.205555] * 1 + [0.203601] * 1 + [0.210877] * 1 + [0.208376] * 1 + [0.216345] * 1 + [0.201452] * 1 + [0.202134] * 1 + [0.207566] * 1 + [0.208152] * 1 + [0.209889] * 1 + [0.211743] * 1 + [0.213188] * 1 + [0.215864] * 1 + [0.22972] * 1 + [0.192515] * 1 + [0.0103233] * 1,
                                        readoutName=ecalEndcapReadoutName,
                                        layerFieldName="layer")

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

# read the crosstalk map
readCrosstalkMap = ReadCaloCrosstalkMap("ReadCrosstalkMap",
                                       fileName="https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/xtalk_neighbours_map_ecalB_thetamodulemerged.root",
                                       OutputLevel=INFO)

# Step 1: merge hits into cells according to initial segmentation
ecalBarrelCellsName = "ECalBarrelCells"
createEcalBarrelCells = CreateCaloCells("CreateECalBarrelCells",
                                        doCellCalibration=True,
                                        calibTool=calibEcalBarrel,
                                        crosstalksTool=readCrosstalkMap,
                                        addCrosstalk=False,
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
ecalEndcapCellsName = "ECalEndcapCells"
createEcalEndcapCells = CreateCaloCells("CreateEcalEndcapCaloCells",
                                        doCellCalibration=True,
                                        calibTool=calibEcalEndcap,
                                        addCellNoise=False,
                                        filterCellNoise=False,
                                        OutputLevel=INFO,
                                        hits=ecalEndcapHitsName,
                                        cells=ecalEndcapCellsName)

# Add to Ecal endcap cells the position information
# (good for physics, all coordinates set properly)
cellPositionEcalEndcapTool = CellPositionsECalEndcapTurbineSegTool(
    "CellPositionsECalEndcap",
    readoutName=ecalEndcapReadoutName,
    OutputLevel=INFO
)
ecalEndcapPositionedCellsName = "ECalEndcapPositionedCells"
createEcalEndcapPositionedCells = CreateCaloCellPositionsFCCee(
    "CreateECalEndcapPositionedCells",
    OutputLevel=INFO
)
createEcalEndcapPositionedCells.positionsTool = cellPositionEcalEndcapTool
createEcalEndcapPositionedCells.hits.Path = ecalEndcapCellsName
createEcalEndcapPositionedCells.positionedHits.Path = ecalEndcapPositionedCellsName

if runHCal:
    # Create cells in HCal
    # 1 - merge hits into cells with the default readout
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

    # 2 - attach positions to the cells
    from Configurables import CellPositionsHCalPhiThetaSegTool
    cellPositionHcalBarrelTool = CellPositionsHCalPhiThetaSegTool(
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

    # 3 - compute new cellID of cells based on new readout - removing row information
    hcalBarrelCellsName2 = "HCalBarrelCells2"
    rewriteHCalBarrel = RewriteBitfield("RewriteHCalBarrel",
                                        # old bitfield (readout)
                                        oldReadoutName=hcalBarrelReadoutName,
                                        # specify which fields are going to be deleted
                                        removeIds=["row"],
                                        # new bitfield (readout), with new segmentation
                                        newReadoutName=hcalBarrelReadoutName2,
                                        debugPrint=10,
                                        OutputLevel=INFO)
    # clusters are needed, with deposit position and cellID in bits
    rewriteHCalBarrel.inhits.Path = hcalBarrelCellsName
    rewriteHCalBarrel.outhits.Path = hcalBarrelCellsName2

    # 4 - attach positions to the new cells
    from Configurables import CellPositionsHCalPhiThetaSegTool
    hcalBarrelPositionedCellsName2 = "HCalBarrelPositionedCells2"
    cellPositionHcalBarrelTool2 = CellPositionsHCalPhiThetaSegTool(
        "CellPositionsHCalBarrel2",
        readoutName=hcalBarrelReadoutName2,
        OutputLevel=INFO
    )
    createHcalBarrelPositionedCells2 = CreateCaloCellPositionsFCCee(
        "CreateHCalBarrelPositionedCells2",
        OutputLevel=INFO
    )
    createHcalBarrelPositionedCells2.positionsTool = cellPositionHcalBarrelTool2
    createHcalBarrelPositionedCells2.hits.Path = hcalBarrelCellsName2
    createHcalBarrelPositionedCells2.positionedHits.Path = hcalBarrelPositionedCellsName2

    # createHcalEndcapCells = CreateCaloCells("CreateHcalEndcapCaloCells",
    #                                    doCellCalibration=True,
    #                                    calibTool=calibHcalEndcap,
    #                                    addCellNoise=False,
    #                                    filterCellNoise=False,
    #                                    OutputLevel=INFO)
    # createHcalEndcapCells.hits.Path="HCalEndcapHits"
    # createHcalEndcapCells.cells.Path="HCalEndcapCells"

else:
    hcalBarrelCellsName = "emptyCaloCells"
    hcalBarrelPositionedCellsName = "emptyCaloCells"
    hcalBarrelCellsName2 = "emptyCaloCells"
    hcalBarrelPositionedCellsName2 = "emptyCaloCells"
    cellPositionHcalBarrelTool = None
    cellPositionHcalBarrelTool2 = None

# Empty cells for parts of calorimeter not implemented yet
createemptycells = CreateEmptyCaloCellsCollection("CreateEmptyCaloCells")
createemptycells.cells.Path = "emptyCaloCells"

# Produce sliding window clusters (ECAL only)
towers = CaloTowerToolFCCee("towers",
                            deltaThetaTower=4 * 0.009817477/4, deltaPhiTower=2 * 2 * pi / 1536.,
                            ecalBarrelReadoutName=ecalBarrelReadoutName,
                            ecalEndcapReadoutName=ecalEndcapReadoutName,
                            ecalFwdReadoutName="",
                            hcalBarrelReadoutName="",
                            hcalExtBarrelReadoutName="",
                            hcalEndcapReadoutName="",
                            hcalFwdReadoutName="",
                            OutputLevel=INFO)
towers.ecalBarrelCells.Path = ecalBarrelPositionedCellsName
towers.ecalEndcapCells.Path = ecalEndcapPositionedCellsName
towers.ecalFwdCells.Path = "emptyCaloCells"

towers.hcalBarrelCells.Path = "emptyCaloCells"
towers.hcalExtBarrelCells.Path = "emptyCaloCells"
towers.hcalEndcapCells.Path = "emptyCaloCells"
towers.hcalFwdCells.Path = "emptyCaloCells"

# Cluster variables
windT = 9
windP = 17
posT = 5
posP = 11
dupT = 7
dupP = 13
finT = 9
finP = 17
# Minimal energy to create a cluster in GeV (FCC-ee detectors have to reconstruct low energy particles)
threshold = 0.040

createClusters = CreateCaloClustersSlidingWindowFCCee("CreateClusters",
                                                      towerTool=towers,
                                                      nThetaWindow=windT, nPhiWindow=windP,
                                                      nThetaPosition=posT, nPhiPosition=posP,
                                                      nThetaDuplicates=dupT, nPhiDuplicates=dupP,
                                                      nThetaFinal=finT, nPhiFinal=finP,
                                                      energyThreshold=threshold,
                                                      energySharingCorrection=False,
                                                      attachCells=True,
                                                      OutputLevel=INFO
                                                      )
createClusters.clusters.Path = "CaloClusters"
createClusters.clusterCells.Path = "CaloClusterCells"

correctCaloClusters = CorrectCaloClusters("correctCaloClusters",
                                          inClusters=createClusters.clusters.Path,
                                          outClusters="Corrected" + createClusters.clusters.Path,
                                          systemIDs=[4],
                                          numLayers=[11],
                                          firstLayerIDs=[0],
                                          lastLayerIDs=[10],
                                          readoutNames=[ecalBarrelReadoutName],
                                          # do not split the following line or it will break scripts that update the values of the corrections
                                          upstreamParameters = [[0.025582045561310333, -0.9524128168665387, -53.10089405478649, 1.283851527438571, -295.30650178662637, -284.8945817377308]],
                                          upstreamFormulas=[
                                              ['[0]+[1]/(x-[2])', '[0]+[1]/(x-[2])']],
                                          # do not split the following line or it will break scripts that update the values of the corrections
                                          downstreamParameters = [[0.0018280333929494054, 0.004932212590963076, 0.8409676097173655, -1.2676690014715288, 0.005347798049886769, 4.161741293789687]],
                                          downstreamFormulas=[
                                              ['[0]+[1]*x', '[0]+[1]/sqrt(x)', '[0]+[1]/x']],
                                          OutputLevel=INFO
                                          )

augmentCaloClusters = AugmentClustersFCCee("augmentCaloClusters",
                                           inClusters=createClusters.clusters.Path,
                                           outClusters="Augmented" + createClusters.clusters.Path,
                                           systemIDs=[4],
                                           systemNames=["EMB"],
                                           numLayers=[11],
                                           readoutNames=[ecalBarrelReadoutName],
                                           layerFieldNames=["layer"],
                                           thetaRecalcWeights=[
                                               # [-1, 3.0, 3.0, 3.0, 4.25, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0]  # -1 : use linear weights
                                               [-1, 3.0, 3.0, 3.0, 4.25, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0]  # -1 : use linear weights
                                           ],
                                           OutputLevel=INFO
                                           )


# Output
out = PodioOutput("out",
                  OutputLevel=INFO)

if runHCal:
    out.outputCommands = ["keep *", "drop emptyCaloCells"]
else:
    out.outputCommands = ["keep *", "drop HCal*", "drop emptyCaloCells"]

if not saveCells:
    out.outputCommands.append("drop ECal*Cells*")
if not saveClusterCells:
    out.outputCommands.append("drop *ClusterCells*")
if not saveHits:
    out.outputCommands.append("drop ECal*Hits*")

out.filename = "./output_evts_" + str(Nevts) + "_pdg_" + str(pdgCode) + "_" + str(momentum) + "_GeV" + "_ThetaMinMax_" + str(thetaMin) + "_" + str(
    thetaMax) + "_PhiMinMax_" + str(phiMin) + "_" + str(phiMax) + "_MagneticField_" + str(magneticField) + "_Noise" + str(addNoise) + ".root"

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
    createEcalEndcapCells,
    createEcalEndcapPositionedCells,
]

if runHCal:
    TopAlg += [
        createHcalBarrelCells,
        createHcalBarrelPositionedCells,
        rewriteHCalBarrel,
        createHcalBarrelPositionedCells2,
        # createHcalEndcapCells
    ]

if doClustering:
    TopAlg += [
        createemptycells,
        createClusters,
    ]

    if applyUpDownstreamCorrections:
        TopAlg += [
            correctCaloClusters,
        ]

    if addShapeParameters:
        TopAlg += [
            augmentCaloClusters,
        ]

TopAlg += [
    out
]

ApplicationMgr(
    TopAlg=TopAlg,
    EvtSel='NONE',
    EvtMax=Nevts,
    ExtSvc=ExtSvc,
    StopOnSignal=True,
)
