#
# IMPORTS
#
from Configurables import ApplicationMgr
#from Configurables import EventCounter
from Configurables import AuditorSvc, ChronoAuditor
# Input/output
from Configurables import k4DataSvc, PodioInput
from Configurables import PodioOutput
# Geometry
from Configurables import GeoSvc
# Create cells
from Configurables import CreateCaloCells
from Configurables import CreateEmptyCaloCellsCollection
# Cell positioning tools
from Configurables import CreateCaloCellPositionsFCCee
from Configurables import CellPositionsECalBarrelModuleThetaSegTool
from Configurables import CellPositionsECalEndcapTurbineSegTool
# Redo segmentation for ECAL
from Configurables import RedoSegmentation
# Change HCAL segmentation
from Configurables import RewriteBitfield
# Apply sampling fraction corrections
from Configurables import CalibrateCaloHitsTool
from Configurables import CalibrateInLayersTool
# Up/down stream correction
from Configurables import CorrectCaloClusters
# SW clustering
from Configurables import CaloTowerToolFCCee
from Configurables import CreateCaloClustersSlidingWindowFCCee
# Topo clustering
from Configurables import CaloTopoClusterInputTool
from Configurables import TopoCaloNeighbours
from Configurables import TopoCaloNoisyCells
from Configurables import CaloTopoClusterFCCee
# Decorate clusters with shower shape parameters
from Configurables import AugmentClustersFCCee
# Read crosstalk map
from Configurables import ReadCaloCrosstalkMap
# Logger
from Gaudi.Configuration import INFO, VERBOSE, DEBUG
# units and physical constants
from GaudiKernel.SystemOfUnits import GeV, tesla, mm
from GaudiKernel.PhysicalConstants import pi, halfpi, twopi
# python libraries
import os
from math import cos, sin, tan

#
# SETTINGS
#

# - general settings
#
inputfile = "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/ALLEGRO_sim.root"
# note - this file probably contains the old ecal endcap segmentation so we disable the endcap digitisation later
Nevts = 50  # -1 means all events
dumpGDML = False

# - what to save in output file
#
# for big productions, save significant space removing hits and cells
# however, hits and cluster cells might be wanted for small productions for detailed event displays
# also, cluster cells are needed for the MVA training
# saveHits = False
# saveCells = False
saveHits = True
saveCells = True
saveClusterCells = True
doCrosstalk = True # switch on/off the crosstalk

# ECAL barrel parameters for digitisation
samplingFraction=[0.37586625991994105] * 1 + [0.13459486704309379] * 1 + [0.142660085165352] * 1 + [0.14768106642302886] * 1 + [0.15205230356024715] * 1 + [0.15593671843591686] * 1 + [0.15969313426201745] * 1 + [0.16334257010426537] * 1 + [0.16546584993953908] * 1 + [0.16930439771304764] * 1 + [0.1725913708958098] * 1
upstreamParameters = [[0.025582045561310333, -0.9524128168665387, -53.10089405478649, 1.283851527438571, -295.30650178662637, -284.8945817377308]]
downstreamParameters = [[0.0018280333929494054, 0.004932212590963076, 0.8409676097173655, -1.2676690014715288, 0.005347798049886769, 4.161741293789687]]
    
ecalBarrelLayers = len(samplingFraction)
resegmentECalBarrel = False

# - parameters for clustering
#
doSWClustering = True
doTopoClustering = True

# calculate cluster energy and barycenter per layer and save it as extra parameters
addShapeParameters = True
ecalBarrelThetaWeights = [-1, 3.0, 3.0, 3.0, 4.25, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0]  # to be recalculated for V03, separately for topo and calo clusters...

#
# ALGORITHMS AND SERVICES SETUP
#

# Input: load the output of the SIM step
evtsvc = k4DataSvc('EventDataSvc')
evtsvc.input = inputfile
input_reader = PodioInput('InputReader')
podioevent = k4DataSvc("EventDataSvc")

# Detector geometry
# prefix all xmls with path_to_detector
# if K4GEO is empty, this should use relative path to working directory
geoservice = GeoSvc("GeoSvc")
path_to_detector = os.environ.get("K4GEO", "")
detectors_to_use = [
    'FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml'
]
geoservice.detectors = [
    os.path.join(path_to_detector, _det) for _det in detectors_to_use
]
geoservice.OutputLevel = INFO

# GDML dump of detector model
if dumpGDML:
    from Configurables import GeoToGdmlDumpSvc
    gdmldumpservice = GeoToGdmlDumpSvc("GeoToGdmlDumpSvc")

# Digitisation (merging hits into cells, EM scale calibration via sampling fractions)

# - ECAL readouts
ecalBarrelReadoutName = "ECalBarrelModuleThetaMerged"
ecalBarrelReadoutName2 = "ECalBarrelModuleThetaMerged2"
ecalEndcapReadoutName = "ECalEndcapTurbine"

hcalBarrelReadoutName = ""
hcalBarrelReadoutName2 = ""
hcalEndcapReadoutName = ""

# - EM scale calibration (sampling fraction)
#   * ECAL barrel
calibEcalBarrel = CalibrateInLayersTool("CalibrateECalBarrel",
                                        samplingFraction=samplingFraction,
                                        readoutName=ecalBarrelReadoutName,
                                        layerFieldName="layer")
#   * ECAL endcap
calibEcalEndcap = CalibrateCaloHitsTool(
        "CalibrateECalEndcap", invSamplingFraction="4.27")  # FIXME: to be updated for ddsim

# Create cells in ECal barrel (needed if one wants to apply cell calibration,
# which is not performed by ddsim)

# read the crosstalk map
readCrosstalkMap = ReadCaloCrosstalkMap("ReadCrosstalkMap",
                                       fileName="https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/xtalk_neighbours_map_ecalB_thetamodulemerged.root",
                                       OutputLevel=INFO)

# - merge hits into cells according to initial segmentation
ecalBarrelCellsName = "ECalBarrelCells"
createEcalBarrelCells = CreateCaloCells("CreateECalBarrelCells",
                                        doCellCalibration=True,
                                        calibTool=calibEcalBarrel,
                                        crosstalksTool=readCrosstalkMap,
                                        addCrosstalk=doCrosstalk,
                                        addCellNoise=False,
                                        filterCellNoise=False,
                                        addPosition=True,
                                        OutputLevel=INFO,
                                        hits=ecalBarrelReadoutName,
                                        cells=ecalBarrelCellsName)

# - add to Ecal barrel cells the position information
#   (good for physics, all coordinates set properly)
cellPositionEcalBarrelTool = CellPositionsECalBarrelModuleThetaSegTool(
    "CellPositionsECalBarrel",
    readoutName=ecalBarrelReadoutName,
    OutputLevel=INFO
)
ecalBarrelPositionedCellsName = ecalBarrelReadoutName + "Positioned"
createEcalBarrelPositionedCells = CreateCaloCellPositionsFCCee(
    "CreateECalBarrelPositionedCells",
    OutputLevel=INFO
)
createEcalBarrelPositionedCells.positionsTool = cellPositionEcalBarrelTool
createEcalBarrelPositionedCells.hits.Path = ecalBarrelCellsName
createEcalBarrelPositionedCells.positionedHits.Path = ecalBarrelPositionedCellsName

# Create cells in ECal endcap (needed if one wants to apply cell calibration,
# which is not performed by ddsim)
#ecalEndcapCellsName = "ECalEndcapCells"
#createEcalEndcapCells = CreateCaloCells("CreateEcalEndcapCaloCells",
#                                        doCellCalibration=True,
#                                        calibTool=calibEcalEndcap,
#                                        addCellNoise=False,
#                                        filterCellNoise=False,
#                                        OutputLevel=INFO,
#                                        hits=ecalEndcapReadoutName,
#                                        cells=ecalEndcapCellsName)

# Add to Ecal endcap cells the position information
# (good for physics, all coordinates set properly)
#cellPositionEcalEndcapTool = CellPositionsECalEndcapTurbineSegTool(
#    "CellPositionsECalEndcap",
#    readoutName=ecalEndcapReadoutName,
#    OutputLevel=INFO
#)
#ecalEndcapPositionedCellsName = "ECalEndcapPositionedCells"
#createEcalEndcapPositionedCells = CreateCaloCellPositionsFCCee(
#    "CreateECalEndcapPositionedCells",
#    OutputLevel=INFO
#)
#createEcalEndcapPositionedCells.positionsTool = cellPositionEcalEndcapTool
#createEcalEndcapPositionedCells.hits.Path = ecalEndcapCellsName
#createEcalEndcapPositionedCells.positionedHits.Path = ecalEndcapPositionedCellsName

hcalBarrelCellsName = "emptyCaloCells"
hcalBarrelPositionedCellsName = "emptyCaloCells"
hcalBarrelCellsName2 = "emptyCaloCells"
hcalBarrelPositionedCellsName2 = "emptyCaloCells"
cellPositionHcalBarrelTool = None
cellPositionHcalBarrelTool2 = None

# Empty cells for parts of calorimeter not implemented yet
createemptycells = CreateEmptyCaloCellsCollection("CreateEmptyCaloCells")
createemptycells.cells.Path = "emptyCaloCells"

# Produce sliding window clusters
if doSWClustering:
    towers = CaloTowerToolFCCee("towers",
                                deltaThetaTower=4 * 0.009817477 / 4, deltaPhiTower=2 * 2 * pi / 1536.,
                                ecalBarrelReadoutName=ecalBarrelReadoutName,
                                #ecalEndcapReadoutName=ecalEndcapReadoutName,
                                ecalEndcapReadoutName="",
                                ecalFwdReadoutName="",
                                hcalBarrelReadoutName=hcalBarrelReadoutName2,
                                hcalExtBarrelReadoutName="",
                                hcalEndcapReadoutName="",
                                hcalFwdReadoutName="",
                                OutputLevel=INFO)
    towers.ecalBarrelCells.Path = ecalBarrelPositionedCellsName
    #towers.ecalEndcapCells.Path = ecalEndcapPositionedCellsName
    towers.ecalEndcapCells.Path = "emptyCaloCells"
    towers.ecalFwdCells.Path = "emptyCaloCells"
    towers.hcalBarrelCells.Path = hcalBarrelPositionedCellsName2
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

    if addShapeParameters:
        augmentCaloClusters = AugmentClustersFCCee("augmentCaloClusters",
                                                   inClusters=createClusters.clusters.Path,
                                                   outClusters="Augmented" + createClusters.clusters.Path,
                                                   systemIDs=[4],
                                                   systemNames=["EMB"],
                                                   numLayers=[ecalBarrelLayers],
                                                   readoutNames=[ecalBarrelReadoutName],
                                                   layerFieldNames=["layer"],
                                                   thetaRecalcWeights=[ecalBarrelThetaWeights],
                                                   OutputLevel=INFO
                                                   )
 
if doTopoClustering:
    # Produce topoclusters (ECAL only or ECAL+HCAL)
    createTopoInput = CaloTopoClusterInputTool("CreateTopoInput",
                                               ecalBarrelReadoutName=ecalBarrelReadoutName,
                                               ecalEndcapReadoutName="",
                                               ecalFwdReadoutName="",
                                               hcalBarrelReadoutName=hcalBarrelReadoutName2,
                                               hcalExtBarrelReadoutName="",
                                               hcalEndcapReadoutName="",
                                               hcalFwdReadoutName="",
                                               OutputLevel=INFO)

    createTopoInput.ecalBarrelCells.Path = ecalBarrelPositionedCellsName
    createTopoInput.ecalEndcapCells.Path = "emptyCaloCells"
    createTopoInput.ecalFwdCells.Path = "emptyCaloCells"
    createTopoInput.hcalBarrelCells.Path = hcalBarrelPositionedCellsName2
    createTopoInput.hcalExtBarrelCells.Path = "emptyCaloCells"
    createTopoInput.hcalEndcapCells.Path = "emptyCaloCells"
    createTopoInput.hcalFwdCells.Path = "emptyCaloCells"
    cellPositionHcalBarrelNoSegTool = None
    cellPositionHcalExtBarrelTool = None

    neighboursMap = "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/neighbours_map_ecalB_thetamodulemerged.root"
    noiseMap = "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/cellNoise_map_electronicsNoiseLevel_ecalB_thetamodulemerged.root"

    readNeighboursMap = TopoCaloNeighbours("ReadNeighboursMap",
                                           fileName=neighboursMap,
                                           OutputLevel=INFO)

    # Noise levels per cell
    readNoisyCellsMap = TopoCaloNoisyCells("ReadNoisyCellsMap",
                                           fileName=noiseMap,
                                           OutputLevel=INFO)

    createTopoClusters = CaloTopoClusterFCCee("CreateTopoClusters",
                                              TopoClusterInput=createTopoInput,
                                              # expects neighbours map from cellid->vec < neighbourIds >
                                              neigboursTool=readNeighboursMap,
                                              # tool to get noise level per cellid
                                              noiseTool=readNoisyCellsMap,
                                              # cell positions tools for all sub - systems
                                              positionsECalBarrelTool=cellPositionEcalBarrelTool,
                                              positionsHCalBarrelTool=cellPositionHcalBarrelTool2,
                                              # positionsHCalBarrelNoSegTool=cellPositionHcalBarrelNoSegTool,
                                              # positionsHCalExtBarrelTool=cellPositionHcalExtBarrelTool,
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


    # Correction below is for EM-only clusters
    # Need something different for EM+HCAL
    if addShapeParameters:
        augmentCaloTopoClusters = AugmentClustersFCCee("augmentCaloTopoClusters",
                                                       inClusters=createTopoClusters.clusters.Path,
                                                       outClusters="Augmented" + createTopoClusters.clusters.Path,
                                                       systemIDs=[4],
                                                       systemNames=["EMB"],
                                                       numLayers=[ecalBarrelLayers],
                                                       readoutNames=[ecalBarrelReadoutName],
                                                       layerFieldNames=["layer"],
                                                       thetaRecalcWeights=[ecalBarrelThetaWeights],
                                                       OutputLevel=INFO)
# Output
out = PodioOutput("out",
                  OutputLevel=INFO)
out.filename = "ALLEGRO_sim_digi_reco.root"

# drop the unpositioned ECal barrel cells
out.outputCommands = ["keep *", "drop HCal*", "drop emptyCaloCells", "drop ECalBarrelCells*"]
out.outputCommands.append("drop %s" % ecalBarrelReadoutName)
out.outputCommands.append("drop %s" % ecalBarrelReadoutName2)
out.outputCommands.append("drop ECalBarrelCellsMerged")

# drop lumi, vertex, DCH, Muons (unless want to keep for event display)
out.outputCommands.append("drop Lumi*")
# out.outputCommands.append("drop Vertex*")
# out.outputCommands.append("drop DriftChamber_simHits*")
out.outputCommands.append("drop MuonTagger*")

# drop hits/positioned cells/cluster cells if desired
if not saveHits:
    out.outputCommands.append("drop %s_contributions" % ecalBarrelReadoutName)
    out.outputCommands.append("drop %s_contributions" % ecalBarrelReadoutName2)
if not saveCells:
    out.outputCommands.append("drop %s" % ecalBarrelPositionedCellsName)
    out.outputCommands.append("drop %s" % ecalBarrelPositionedCellsName2)

if not saveClusterCells:
    out.outputCommands.append("drop Calo*ClusterCells*")

# CPU information
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
out.AuditExecute = True

# Event counter
#event_counter = EventCounter('event_counter')
#event_counter.Frequency = 10

# Configure list of external services
ExtSvc = [geoservice, podioevent, audsvc]
if dumpGDML:
    ExtSvc += [gdmldumpservice]

# Setup alg sequence
TopAlg = [
#    event_counter,
    input_reader,
    createEcalBarrelCells,
    createEcalBarrelPositionedCells,
#    createEcalEndcapCells,
#    createEcalEndcapPositionedCells
]
createEcalBarrelCells.AuditExecute = True
createEcalBarrelPositionedCells.AuditExecute = True
#createEcalEndcapCells.AuditExecute = True

if resegmentECalBarrel:
    TopAlg += [
        resegmentEcalBarrelTool,
        createEcalBarrelCells2,
        createEcalBarrelPositionedCells2,
    ]
    resegmentEcalBarrelTool.AuditExecute = True
    createEcalBarrelCells2.AuditExecute = True
    createEcalBarrelPositionedCells2.AuditExecute = True

if doSWClustering or doTopoClustering:
    TopAlg += [createemptycells]
    createemptycells.AuditExecute = True
    
    if doSWClustering:
        TopAlg += [createClusters]
        createClusters.AuditExecute = True

        if addShapeParameters:
            TopAlg += [augmentCaloClusters]
            augmentCaloClusters.AuditExecute = True

    if doTopoClustering:
        TopAlg += [createTopoClusters]
        createTopoClusters.AuditExecute = True

        if addShapeParameters:
            TopAlg += [augmentCaloTopoClusters]
            augmentCaloTopoClusters.AuditExecute = True

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

