#
# IMPORTS
#
from Configurables import ApplicationMgr
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
# Redo segmentation for ECAL and HCAL
from Configurables import RedoSegmentation
# Read noise values from file and generate noise in cells
from Configurables import NoiseCaloCellsVsThetaFromFileTool
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
inputfile = "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/ALLEGRO_sim.root"  # input file produced with ddsim
Nevts = 2                      # -1 means all events
addNoise = True                # add noise or not to the cell energy
dumpGDML = False                # create GDML file of detector model

# - what to save in output file
#
# for big productions, save significant space removing hits and cells
# however, hits and cluster cells might be wanted for small productions for detailed event displays
# cluster cells are not needed for the training of the MVA energy regression nor the photon ID since needed quantities are stored in cluster shapeParameters
# saveHits = False
# saveCells = False
# saveClusterCells = False
saveHits = True
saveCells = True
saveClusterCells = True

# ECAL barrel parameters for digitisation
samplingFraction=[0.3800493723322256] * 1 + [0.13494147915064658] * 1 + [0.142866851721152] * 1 + [0.14839315921940666] * 1 + [0.15298362570665006] * 1 + [0.15709704561942747] * 1 + [0.16063717490147533] * 1 + [0.1641723795419055] * 1 + [0.16845490287689746] * 1 + [0.17111520115997653] * 1 + [0.1730605163148862] * 1
upstreamParameters = [[0.028158491043365624, -1.564259408365951, -76.52312805346982, 0.7442903558010191, -34.894692961350195, -74.19340877431723]]
downstreamParameters = [[0.00010587711361028165, 0.0052371999097777355, 0.69906696456064, -0.9348243433360095, -0.0364714212117143, 8.360401126995626]]
    
ecalBarrelLayers = len(samplingFraction)
resegmentECalBarrel = False

# - parameters for clustering
#
doSWClustering = True
doTopoClustering = True

# cluster energy corrections
# simple parametrisations of up/downstream losses
# not to be applied for ECAL+HCAL clustering
applyUpDownstreamCorrections = False

# calculate cluster energy and barycenter per layer and save it as extra parameters
addShapeParameters = True
ecalBarrelThetaWeights = [-1, 3.0, 3.0, 3.0, 4.25, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0]  # to be recalculated for V03, separately for topo and calo clusters...

#
# ALGORITHMS AND SERVICES SETUP
#

# Input: load the output of the SIM step
podioevent = k4DataSvc('EventDataSvc')
podioevent.input = inputfile
input_reader = PodioInput('InputReader')


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
ecalBarrelReadoutName = "ECalBarrelModuleThetaMerged"     # barrel, original segmentation (baseline)
ecalBarrelReadoutName2 = "ECalBarrelModuleThetaMerged2"   # barrel, after re-segmentation (for optimisation studies)
ecalEndcapReadoutName = "ECalEndcapTurbine"               # endcap, turbine-like (baseline)
# - HCAL readouts
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
calibEcalEndcap = CalibrateInLayersTool("CalibrateECalEndcap",
                                        samplingFraction = [0.16419] * 1 + [0.192898] * 1 + [0.18783] * 1 + [0.193203] * 1 + [0.193928] * 1 + [0.192286] * 1 + [0.199959] * 1 + [0.200153] * 1 + [0.212635] * 1 + [0.180345] * 1 + [0.18488] * 1 + [0.194762] * 1 + [0.197775] * 1 + [0.200504] * 1 + [0.205555] * 1 + [0.203601] * 1 + [0.210877] * 1 + [0.208376] * 1 + [0.216345] * 1 + [0.201452] * 1 + [0.202134] * 1 + [0.207566] * 1 + [0.208152] * 1 + [0.209889] * 1 + [0.211743] * 1 + [0.213188] * 1 + [0.215864] * 1 + [0.22972] * 1 + [0.192515] * 1 + [0.0103233] * 1,
                                        readoutName=ecalEndcapReadoutName,
                                        layerFieldName="layer")

# Create cells in ECal barrel (needed if one wants to apply cell calibration,
# which is not performed by ddsim)
# - merge hits into cells according to initial segmentation
ecalBarrelCellsName = "ECalBarrelCells"
if not addNoise: # Create ECAL barrel cells without adding noise
    createEcalBarrelCells = CreateCaloCells("CreateECalBarrelCells",
                                            doCellCalibration=True,
                                            calibTool=calibEcalBarrel,
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

# -  now, if we want to also save cells with coarser granularity:
if resegmentECalBarrel:
    # 2. step - rewrite the cellId using the merged theta-module segmentation
    # (merging several modules and severla theta readout cells).
    # Add noise at this step if you derived the noise already assuming merged cells
    # Step 2a: compute new cellID of cells based on new readout
    # (merged module-theta segmentation with variable merging vs layer)
    resegmentEcalBarrelTool = RedoSegmentation("ReSegmentationEcal",
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
    createEcalBarrelPositionedCells2.positionedHits.Path = ecalBarrelReadoutName2 + "Positioned"


# Create cells in ECal endcap (needed if one wants to apply cell calibration,
# which is not performed by ddsim)
ecalEndcapCellsName = "ECalEndcapCells"
createEcalEndcapCells = CreateCaloCells("CreateEcalEndcapCaloCells",
                                        doCellCalibration=True,
                                        calibTool=calibEcalEndcap,
                                        addCellNoise=False,
                                        filterCellNoise=False,
                                        OutputLevel=INFO,
                                        hits=ecalEndcapReadoutName,
                                        cells=ecalEndcapCellsName)

# Add to Ecal endcap cells the position information
# (good for physics, all coordinates set properly)
cellPositionEcalEndcapTool = CellPositionsECalEndcapTurbineSegTool(
    "CellPositionsECalEndcap",
    readoutName=ecalEndcapReadoutName,
    OutputLevel=INFO
)
ecalEndcapPositionedCellsName = ecalEndcapReadoutName + "Positioned"
createEcalEndcapPositionedCells = CreateCaloCellPositionsFCCee(
    "CreateECalEndcapPositionedCells",
    OutputLevel=INFO
)
createEcalEndcapPositionedCells.positionsTool = cellPositionEcalEndcapTool
createEcalEndcapPositionedCells.hits.Path = ecalEndcapCellsName
createEcalEndcapPositionedCells.positionedHits.Path = ecalEndcapPositionedCellsName


if addNoise:
    # FIXME Input histograms of ECAL noise needs to be updated to match the latest ALLEGRO geometry
    ecalBarrelNoisePath = "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/elecNoise_ecalBarrelFCCee_theta.root"
    ecalBarrelNoiseRMSHistName = "h_elecNoise_fcc_"
    from Configurables import NoiseCaloCellsVsThetaFromFileTool
    noiseBarrel = NoiseCaloCellsVsThetaFromFileTool("NoiseBarrel",
                                                    cellPositionsTool=cellPositionEcalBarrelTool,
                                                    readoutName=ecalBarrelReadoutName,
                                                    noiseFileName=ecalBarrelNoisePath,
                                                    elecNoiseRMSHistoName=ecalBarrelNoiseRMSHistName,
                                                    setNoiseOffset=False,
                                                    activeFieldName="layer",
                                                    addPileup=False,
                                                    filterNoiseThreshold=0,
                                                    numRadialLayers=11,
                                                    scaleFactor=1 / 1000.,  # MeV to GeV
                                                    OutputLevel=INFO)

    # barrel geometry tool is migrated to match ALLEGRO v3 segmentation
    from Configurables import TubeLayerModuleThetaMergedCaloTool
    barrelGeometry = TubeLayerModuleThetaMergedCaloTool("EcalBarrelGeo",
                                             readoutName=ecalBarrelReadoutName,
                                             activeVolumeName="LAr_sensitive",
                                             activeFieldName="layer",
                                             activeVolumesNumber=11,
                                             fieldNames=["system"],
                                             fieldValues=[4])
    
    # cells with noise not filtered
    createEcalBarrelCells = CreateCaloCells("CreateECalBarrelCellsNoise",
                                                 doCellCalibration=True,
                                                 calibTool=calibEcalBarrel,
                                                 addCellNoise=True,
                                                 filterCellNoise=False,
                                                 noiseTool=noiseBarrel,
                                                 addPosition=True,
                                                 geometryTool=barrelGeometry,
                                                 OutputLevel=INFO,
                                                 hits=ecalBarrelReadoutName,
                                                 cells=ecalBarrelCellsName)

    # cells with noise filtered
    # createEcalBarrelCellsNoise = CreateCaloCells("CreateECalBarrelCellsNoise_filtered",
    #                                              doCellCalibration=False,
    #                                              addCellNoise=True,
    #                                              filterCellNoise=True,
    #                                              OutputLevel=INFO,
    #                                              hits="ECalBarrelCellsStep2",
    #                                              noiseTool=noiseBarrel,
    #                                              geometryTool=barrelGeometry,
    #                                              cells=EcalBarrelCellsName)

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
                                ecalEndcapReadoutName=ecalEndcapReadoutName,
                                ecalFwdReadoutName="",
                                hcalBarrelReadoutName=hcalBarrelReadoutName2,
                                hcalExtBarrelReadoutName="",
                                hcalEndcapReadoutName="",
                                hcalFwdReadoutName="",
                                OutputLevel=INFO)
    towers.ecalBarrelCells.Path = ecalBarrelPositionedCellsName
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

    if applyUpDownstreamCorrections:
        correctCaloClusters = CorrectCaloClusters("CorrectCaloClusters",
                                                  inClusters=createClusters.clusters.Path,
                                                  outClusters="Corrected" + createClusters.clusters.Path,
                                                  systemIDs=[4],
                                                  numLayers=[ecalBarrelLayers],
                                                  firstLayerIDs=[0],
                                                  lastLayerIDs=[ecalBarrelLayers-1],
                                                  readoutNames=[ecalBarrelReadoutName],
                                                  upstreamParameters=upstreamParameters,
                                                  upstreamFormulas=[
                                                      ['[0]+[1]/(x-[2])', '[0]+[1]/(x-[2])']],
                                                  downstreamParameters=downstreamParameters,
                                                  downstreamFormulas=[
                                                      ['[0]+[1]*x', '[0]+[1]/sqrt(x)', '[0]+[1]/x']],
                                                  OutputLevel=INFO
                                                  )

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
                                                   do_photon_shapeVar=False,
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
    if applyUpDownstreamCorrections:
        correctCaloTopoClusters = CorrectCaloClusters(
            "CorrectCaloTopoClusters",
            inClusters=createTopoClusters.clusters.Path,
            outClusters="Corrected" + createTopoClusters.clusters.Path,
            systemIDs=[4],
            numLayers=[ecalBarrelLayers],
            firstLayerIDs=[0],
            lastLayerIDs=[ecalBarrelLayers-1],
            readoutNames=[ecalBarrelReadoutName],
            # do not split the following line or it will break scripts that update the values of the corrections
            upstreamParameters=upstreamParameters,
            upstreamFormulas=[['[0]+[1]/(x-[2])', '[0]+[1]/(x-[2])']],
            # do not split the following line or it will break scripts that update the values of the corrections
            downstreamParameters=downstreamParameters,
            downstreamFormulas=[['[0]+[1]*x', '[0]+[1]/sqrt(x)', '[0]+[1]/x']],
            OutputLevel=INFO
        )

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
                                                       do_photon_shapeVar=False,
                                                       OutputLevel=INFO)

# Output
out = PodioOutput("out",
                  OutputLevel=INFO)
out.filename = "ALLEGRO_sim_digi_reco.root"

# drop the unpositioned ECal and HCal barrel and endcap cells
out.outputCommands = ["keep *",
                    "drop emptyCaloCells",
                    "drop %s" % ecalBarrelCellsName,
                    "drop %s" % ecalEndcapCellsName]

# drop the uncalibrated cells
out.outputCommands.append("drop %s" % ecalBarrelReadoutName)
out.outputCommands.append("drop %s" % ecalBarrelReadoutName2)
out.outputCommands.append("drop %s" % ecalEndcapReadoutName)

# drop the intermediate ecal barrel cells in case of a resegmentation
if resegmentECalBarrel:
    out.outputCommands.append("drop ECalBarrelCellsMerged")
    out.outputCommands.append("drop %s" % ecalBarrelCellsName2)

# drop lumi, vertex, DCH, Muons (unless want to keep for event display)
out.outputCommands.append("drop Lumi*")
# out.outputCommands.append("drop Vertex*")
# out.outputCommands.append("drop DriftChamber_simHits*")
out.outputCommands.append("drop MuonTagger*")

# drop hits/positioned cells/cluster cells if desired
if not saveHits:
    out.outputCommands.append("drop *%sContributions" % ecalBarrelReadoutName)
    out.outputCommands.append("drop *%sContributions" % ecalBarrelReadoutName2)
    out.outputCommands.append("drop *%sContributions" % ecalEndcapReadoutName)
if not saveCells:
    out.outputCommands.append("drop %s" % ecalBarrelPositionedCellsName)
    out.outputCommands.append("drop %s" % ecalEndcapPositionedCellsName)
    if resegmentECalBarrel:
        out.outputCommands.append("drop %s" % ecalBarrelPositionedCellsName2)
if not saveClusterCells:
    out.outputCommands.append("drop Calo*ClusterCells*")

# if we decorate the clusters, we can drop the non-decorated ones
if addShapeParameters:
    out.outputCommands.append("drop %s" % augmentCaloClusters.inClusters)

# CPU information
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
out.AuditExecute = True

# Configure list of external services
ExtSvc = [geoservice, podioevent, audsvc]
if dumpGDML:
    ExtSvc += [gdmldumpservice]

# Setup alg sequence
TopAlg = [
    input_reader,
    createEcalBarrelCells,
    createEcalBarrelPositionedCells,
]
createEcalBarrelCells.AuditExecute = True
createEcalBarrelPositionedCells.AuditExecute = True

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

        if applyUpDownstreamCorrections:
            TopAlg += [correctCaloClusters]
            correctCaloClusters.AuditExecute = True

        if addShapeParameters:
            TopAlg += [augmentCaloClusters]
            augmentCaloClusters.AuditExecute = True

    if doTopoClustering:
        TopAlg += [createTopoClusters]
        createTopoClusters.AuditExecute = True
        
        if applyUpDownstreamCorrections:
            TopAlg += [correctCaloTopoClusters]
            correctCaloTopoClusters.AuditExecute = True

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
