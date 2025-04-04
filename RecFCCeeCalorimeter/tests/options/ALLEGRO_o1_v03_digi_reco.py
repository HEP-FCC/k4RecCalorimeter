#
# COMMON IMPORTS
#

# Logger
from Gaudi.Configuration import INFO  # DEBUG, VERBOSE
# units and physical constants
from GaudiKernel.PhysicalConstants import pi


#
# SETTINGS
#

# - general settings
#
inputfile = "ALLEGRO_sim_ee_z_qq.root"  # input file produced with ddsim
Nevts = -1                              # -1 means all events
addNoise = True                         # add noise or not to the cell energy
addCrosstalk = True                     # switch on/off the crosstalk
dumpGDML = False                        # create GDML file of detector model
runHCal = True                          # if false, it will produce only ECAL clusters. if true, it will also produce ECAL+HCAL clusters

# - what to save in output file
#
# always drop uncalibrated cells, except for tests and debugging
# dropUncalibratedCells = True
dropUncalibratedCells = False

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
ecalBarrelSamplingFraction = [0.3800493723322256] * 1 + [0.13494147915064658] * 1 + [0.142866851721152] * 1 + [0.14839315921940666] * 1 + [0.15298362570665006] * 1 + [0.15709704561942747] * 1 + [0.16063717490147533] * 1 + [0.1641723795419055] * 1 + [0.16845490287689746] * 1 + [0.17111520115997653] * 1 + [0.1730605163148862] * 1
ecalBarrelUpstreamParameters = [[0.028158491043365624, -1.564259408365951, -76.52312805346982, 0.7442903558010191, -34.894692961350195, -74.19340877431723]]
ecalBarrelDownstreamParameters = [[0.00010587711361028165, 0.0052371999097777355, 0.69906696456064, -0.9348243433360095, -0.0364714212117143, 8.360401126995626]]

ecalBarrelLayers = len(ecalBarrelSamplingFraction)
resegmentECalBarrel = False

# - parameters for clustering
#
doSWClustering = True
doTopoClustering = True

# cluster energy corrections
# simple parametrisations of up/downstream losses for ECAL-only clusters
# not to be applied for ECAL+HCAL clustering
# superseded by MVA calibration, but turned on here for the purpose of testing that the code is not broken - will end up in separate cluster collection
applyUpDownstreamCorrections = True

# BDT regression from total cluster energy and fraction of energy in each layer (after correction for sampling fraction)
# not to be applied (yet) for ECAL+HCAL clustering (MVA trained only on ECAL so far)
applyMVAClusterEnergyCalibration = True

# calculate cluster energy and barycenter per layer and save it as extra parameters
addShapeParameters = True
ecalBarrelThetaWeights = [-1, 3.0, 3.0, 3.0, 4.25, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0]  # to be recalculated for V03, separately for topo and calo clusters...

# run photon ID algorithm
# not run by default in production, but to be turned on here for the purpose of testing that the code is not broken
# currently off till we provide the onnx files
runPhotonIDTool = False
logEWeightInPhotonID = False

# resolved pi0 reconstruction by cluster pairing
addPi0RecoTool = True

#
# ALGORITHMS AND SERVICES SETUP
#

# Input: load the output of the SIM step
from Configurables import k4DataSvc, PodioInput
podioevent = k4DataSvc('EventDataSvc')
podioevent.input = inputfile
input_reader = PodioInput('InputReader')


# Detector geometry
# prefix all xmls with path_to_detector
# if K4GEO is empty, this should use relative path to working directory
from Configurables import GeoSvc
import os
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
ecalBarrelReadoutName = "ECalBarrelModuleThetaMerged"      # barrel, original segmentation (baseline)
ecalBarrelReadoutName2 = "ECalBarrelModuleThetaMerged2"    # barrel, after re-segmentation (for optimisation studies)
ecalEndcapReadoutName = "ECalEndcapTurbine"                # endcap, turbine-like (baseline)
# - HCAL readouts
if runHCal:
    hcalBarrelReadoutName = "HCalBarrelReadout"            # barrel, original segmentation (HCalPhiTheta)
    hcalEndcapReadoutName = "HCalEndcapReadout"            # endcap, original segmentation (HCalPhiTheta)
else:
    hcalBarrelReadoutName = ""
    hcalEndcapReadoutName = ""

# - EM scale calibration (sampling fraction)
from Configurables import CalibrateInLayersTool
#   * ECAL barrel
calibEcalBarrel = CalibrateInLayersTool("CalibrateECalBarrel",
                                        samplingFraction=ecalBarrelSamplingFraction,
                                        readoutName=ecalBarrelReadoutName,
                                        layerFieldName="layer")
#   * ECAL endcap
# the turbine endcap has calibration "layers" in the both the z and radial
# directions, for each of the three wheels.  So the total number of layers
# is given by:
#
#   ECalEndcapNumCalibZLayersWheel1*ECalEndcapNumCalibRhoLayersWheel1
#  +ECalEndcapNumCalibZLayersWheel2*ECalEndcapNumCalibRhoLayersWheel2
#  +ECalEndcapNumCalibZLayersWheel3*ECalEndcapNumCalibRhoLayersWheel3
#
# which in the current design is 5*10+1*14+1*34 = 98
# NB some cells near the inner and outer edges of the calorimeter are difficult
# to calibrate as they are not part of the core of well-contained showers.
# The calibrated values can be <0 or >1 for such cells, so these nonsenical
# numbers are replaced by 1
calibEcalEndcap = CalibrateInLayersTool("CalibrateECalEndcap",

                                        samplingFraction = [0.0390585] * 1+ [0.463715] * 1+ [0.0853304] * 1+ [1] * 1+ [0.0442381] * 1+ [0.0552701] * 1+ [0.245704] * 1+ [0.137059] * 1+ [0.238013] * 1+ [0.147364] * 1+ [0.0788765] * 1+ [0.219239] * 1+ [0.151537] * 1+ [0.26711] * 1+ [0.0853741] * 1+ [0.0995584] * 1+ [0.206757] * 1+ [0.16663] * 1+ [0.238659] * 1+ [0.0978915] * 1+ [0.121322] * 1+ [0.204324] * 1+ [0.175887] * 1+ [0.244135] * 1+ [0.0994313] * 1+ [0.141767] * 1+ [0.205447] * 1+ [0.186003] * 1+ [0.237895] * 1+ [0.108545] * 1+ [0.16853] * 1+ [0.205409] * 1+ [0.195066] * 1+ [0.234766] * 1+ [0.109689] * 1+ [0.152743] * 1+ [0.202067] * 1+ [0.199587] * 1+ [0.22892] * 1+ [0.117163] * 1+ [0.247339] * 1+ [0.219485] * 1+ [0.217142] * 1+ [0.248389] * 1+ [0.125487] * 1+ [0.10275] * 1+ [0.185187] * 1+ [0.18939] * 1+ [0.171225] * 1+ [0.109632] * 1+ [0.159225] * 1+ [0.172072] * 1+ [0.174928] * 1+ [0.180492] * 1+ [0.184174] * 1+ [0.18843] * 1+ [0.191688] * 1+ [0.194472] * 1+ [0.198021] * 1+ [0.200064] * 1+ [0.20313] * 1+ [0.203984] * 1+ [0.208658] * 1+ [0.198391] * 1+ [0.153334] * 1+ [0.175712] * 1+ [0.169342] * 1+ [0.176852] * 1+ [0.175731] * 1+ [0.181481] * 1+ [0.17978] * 1+ [0.18426] * 1+ [0.185028] * 1+ [0.186753] * 1+ [0.190321] * 1+ [0.189175] * 1+ [0.193018] * 1+ [0.193594] * 1+ [0.195598] * 1+ [0.1969] * 1+ [0.197713] * 1+ [0.201379] * 1+ [0.193965] * 1+ [0.223362] * 1+ [0.166243] * 1+ [0.373141] * 1+ [0.112825] * 1+ [0.469126] * 1+ [1] * 1+ [0.0954407] * 1+ [0.0706573] * 1+ [1] * 1+ [1] * 1+ [0.112354] * 1+ [1] * 1+ [1] * 1+ [1] * 1+ [1] * 1,
                                        readoutName=ecalEndcapReadoutName,
                                        layerFieldName="layer")

if runHCal:
    from Configurables import CalibrateCaloHitsTool
    # HCAL barrel
    calibHCalBarrel = CalibrateCaloHitsTool(
        "CalibrateHCalBarrel", invSamplingFraction="29.4202")
    # HCAL endcap
    calibHCalEndcap = CalibrateCaloHitsTool(
        "CalibrateHCalEndcap", invSamplingFraction="29.4202")  # FIXME: to be updated for ddsim

# - cell positioning tools
from Configurables import CellPositionsECalBarrelModuleThetaSegTool
cellPositionEcalBarrelTool = CellPositionsECalBarrelModuleThetaSegTool(
    "CellPositionsECalBarrel",
    readoutName=ecalBarrelReadoutName,
    OutputLevel=INFO
)
# the noise tool needs the positioning tool, but if I reuse the previous one the code crashes..
cellPositionEcalBarrelToolForNoise = CellPositionsECalBarrelModuleThetaSegTool(
    "CellPositionsECalBarrelForNoise",
    readoutName=ecalBarrelReadoutName,
    OutputLevel=INFO
)
if resegmentECalBarrel:
    cellPositionEcalBarrelTool2 = CellPositionsECalBarrelModuleThetaSegTool(
        "CellPositionsECalBarrel2",
        readoutName=ecalBarrelReadoutName2,
        OutputLevel=INFO
    )

from Configurables import CellPositionsECalEndcapTurbineSegTool
cellPositionEcalEndcapTool = CellPositionsECalEndcapTurbineSegTool(
    "CellPositionsECalEndcap",
    readoutName=ecalEndcapReadoutName,
    OutputLevel=INFO
)

if runHCal:
    from Configurables import CellPositionsHCalPhiThetaSegTool
    cellPositionHCalBarrelTool = CellPositionsHCalPhiThetaSegTool(
        "CellPositionsHCalBarrel",
        readoutName=hcalBarrelReadoutName,
        detectorName="HCalBarrel",
        OutputLevel=INFO
    )
    cellPositionHCalEndcapTool = CellPositionsHCalPhiThetaSegTool(
        "CellPositionsHCalEndcap",
        readoutName=hcalEndcapReadoutName,
        detectorName="HCalThreePartsEndcap",
        numLayersHCalThreeParts=[6, 9, 22],
        OutputLevel=INFO
    )

# - crosstalk tool
if addCrosstalk:
    from Configurables import ReadCaloCrosstalkMap
    # read the crosstalk map
    readCrosstalkMap = ReadCaloCrosstalkMap("ReadCrosstalkMap",
                                            fileName="https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/xtalk_neighbours_map_ecalB_thetamodulemerged.root",
                                            OutputLevel=INFO)
else:
    readCrosstalkMap = None

# Create cells in ECal barrel (calibrated and positioned - optionally with xtalk and noise added)
# from uncalibrated cells (+cellID info) from ddsim
ecalBarrelPositionedCellsName = ecalBarrelReadoutName + "Positioned"
ecalBarrelLinks = ecalBarrelPositionedCellsName + "SimCaloHitLinks"
from Configurables import CreatePositionedCaloCells
createEcalBarrelCells = CreatePositionedCaloCells("CreatePositionedECalBarrelCells",
                                                  doCellCalibration=True,
                                                  positionsTool=cellPositionEcalBarrelTool,
                                                  calibTool=calibEcalBarrel,
                                                  crosstalkTool=readCrosstalkMap,
                                                  addCrosstalk=addCrosstalk,
                                                  addCellNoise=False,
                                                  filterCellNoise=False,
                                                  OutputLevel=INFO,
                                                  hits=ecalBarrelReadoutName,
                                                  cells=ecalBarrelPositionedCellsName,
                                                  links=ecalBarrelLinks
                                                  )

# -  now, if we want to also save cells with coarser granularity:
if resegmentECalBarrel:
    # rewrite the cellId using the merged theta-module segmentation
    # (merging several modules and severla theta readout cells).
    # Add noise at this step if you derived the noise already assuming merged cells
    # Step a: compute new cellID of cells based on new readout
    # (merged module-theta segmentation with variable merging vs layer)
    from Configurables import RedoSegmentation
    resegmentEcalBarrelTool = RedoSegmentation("ReSegmentationEcal",
                                               # old bitfield (readout)
                                               oldReadoutName=ecalBarrelReadoutName,
                                               # specify which fields are going to be altered (deleted/rewritten)
                                               oldSegmentationIds=["module", "theta"],
                                               # new bitfield (readout), with new segmentation (merged modules and theta cells)
                                               newReadoutName=ecalBarrelReadoutName2,
                                               OutputLevel=INFO,
                                               debugPrint=200,
                                               inhits=ecalBarrelPositionedCellsName,
                                               outhits="ECalBarrelCellsMerged")

    # Step b: merge new cells with same cellID together
    # do not apply cell calibration again since cells were already
    # calibrated in Step 1
    # noise and xtalk off assuming they were applied earlier
    ecalBarrelPositionedCellsName2 = ecalBarrelReadoutName2 + "Positioned"
    ecalBarrelLinks2 = ecalBarrelPositionedCellsName2 + "SimCaloHitLinks"
    createEcalBarrelCells2 = CreatePositionedCaloCells("CreatePositionedECalBarrelCells2",
                                                       doCellCalibration=False,
                                                       positionsTool=cellPositionEcalBarrelTool2,
                                                       calibTool=None,
                                                       crosstalkTool=None,
                                                       addCrosstalk=False,
                                                       addCellNoise=False,
                                                       filterCellNoise=False,
                                                       OutputLevel=INFO,
                                                       hits="ECalBarrelCellsMerged",
                                                       cells=ecalBarrelPositionedCellsName2,
                                                       links=ecalBarrelLinks2)

# Create cells in ECal endcap (needed if one wants to apply cell calibration,
# which is not performed by ddsim)
ecalEndcapPositionedCellsName = ecalEndcapReadoutName + "Positioned"
ecalEndcapLinks = ecalEndcapPositionedCellsName + "SimCaloHitLinks"
createEcalEndcapCells = CreatePositionedCaloCells("CreatePositionedECalEndcapCells",
                                                  doCellCalibration=True,
                                                  positionsTool=cellPositionEcalEndcapTool,
                                                  calibTool=calibEcalEndcap,
                                                  crosstalkTool=None,
                                                  addCrosstalk=False,
                                                  addCellNoise=False,
                                                  filterCellNoise=False,
                                                  OutputLevel=INFO,
                                                  hits=ecalEndcapReadoutName,
                                                  cells=ecalEndcapPositionedCellsName,
                                                  links=ecalEndcapLinks)

if addNoise:
    ecalBarrelNoisePath = "elecNoise_ecalBarrelFCCee_theta.root"
    ecalBarrelNoiseRMSHistName = "h_elecNoise_fcc_"
    from Configurables import NoiseCaloCellsVsThetaFromFileTool
    noiseBarrel = NoiseCaloCellsVsThetaFromFileTool("NoiseBarrel",
                                                    cellPositionsTool=cellPositionEcalBarrelToolForNoise,
                                                    readoutName=ecalBarrelReadoutName,
                                                    noiseFileName=ecalBarrelNoisePath,
                                                    elecNoiseRMSHistoName=ecalBarrelNoiseRMSHistName,
                                                    setNoiseOffset=False,
                                                    activeFieldName="layer",
                                                    addPileup=False,
                                                    filterNoiseThreshold=1,
                                                    numRadialLayers=11,
                                                    scaleFactor=1 / 1000.,  # MeV to GeV
                                                    OutputLevel=INFO)

    from Configurables import TubeLayerModuleThetaCaloTool
    barrelGeometry = TubeLayerModuleThetaCaloTool("EcalBarrelGeo",
                                                  readoutName=ecalBarrelReadoutName,
                                                  activeVolumeName="LAr_sensitive",
                                                  activeFieldName="layer",
                                                  activeVolumesNumber=11,
                                                  fieldNames=["system"],
                                                  fieldValues=[4],
                                                  OutputLevel=INFO)

    # cells with noise not filtered
    ecalBarrelCellsNoiseLinks = ecalBarrelPositionedCellsName + "WithNoise" + "SimCaloHitLinks"
    createEcalBarrelCellsNoise = CreatePositionedCaloCells("CreatePositionedECalBarrelCellsWithNoise",
                                                           doCellCalibration=True,
                                                           positionsTool=cellPositionEcalBarrelTool,
                                                           calibTool=calibEcalBarrel,
                                                           addCellNoise=True,
                                                           filterCellNoise=False,
                                                           noiseTool=noiseBarrel,
                                                           geometryTool=barrelGeometry,
                                                           OutputLevel=INFO,
                                                           hits=ecalBarrelReadoutName,  # uncalibrated & unpositioned cells without noise
                                                           cells=ecalBarrelPositionedCellsName + "WithNoise",
                                                           links=ecalBarrelCellsNoiseLinks)

    # cells with noise filtered
    ecalBarrelCellsNoiseFilteredLinks = ecalBarrelPositionedCellsName + "WithNoiseFiltered" + "SimCaloHitLinks"
    createEcalBarrelCellsNoiseFiltered = CreatePositionedCaloCells("CreateECalBarrelCellsWithNoiseFiltered",
                                                                   doCellCalibration=True,
                                                                   calibTool=calibEcalBarrel,
                                                                   positionsTool=cellPositionEcalBarrelTool,
                                                                   addCellNoise=True,
                                                                   filterCellNoise=True,
                                                                   noiseTool=noiseBarrel,
                                                                   geometryTool=barrelGeometry,
                                                                   OutputLevel=INFO,
                                                                   hits=ecalBarrelReadoutName,  # uncalibrated & unpositioned cells without noise
                                                                   cells=ecalBarrelPositionedCellsName + "WithNoiseFiltered",
                                                                   links=ecalBarrelCellsNoiseFilteredLinks
                                                                   )

if runHCal:
    # Apply calibration and positioning to cells in HCal barrel
    hcalBarrelPositionedCellsName = hcalBarrelReadoutName + "Positioned"
    hcalBarrelLinks = hcalBarrelPositionedCellsName + "SimCaloHitLinks"
    createHCalBarrelCells = CreatePositionedCaloCells("CreateHCalBarrelCells",
                                                      doCellCalibration=True,
                                                      calibTool=calibHCalBarrel,
                                                      positionsTool=cellPositionHCalBarrelTool,
                                                      addCellNoise=False,
                                                      filterCellNoise=False,
                                                      hits=hcalBarrelReadoutName,
                                                      cells=hcalBarrelPositionedCellsName,
                                                      links=hcalBarrelLinks,
                                                      OutputLevel=INFO)

    # Create cells in HCal endcap
    hcalEndcapPositionedCellsName = hcalEndcapReadoutName + "Positioned"
    hcalEndcapLinks = hcalEndcapPositionedCellsName + "SimCaloHitLinks"
    createHCalEndcapCells = CreatePositionedCaloCells("CreateHCalEndcapCells",
                                                      doCellCalibration=True,
                                                      calibTool=calibHCalEndcap,
                                                      addCellNoise=False,
                                                      filterCellNoise=False,
                                                      positionsTool=cellPositionHCalEndcapTool,
                                                      OutputLevel=INFO,
                                                      hits=hcalEndcapReadoutName,
                                                      cells=hcalEndcapPositionedCellsName,
                                                      links=hcalEndcapLinks)

else:
    hcalBarrelPositionedCellsName = "emptyCaloCells"
    hcalEndcapPositionedCellsName = "emptyCaloCells"
    hcalBarrelLinks = ""
    hcalEndcapLinks = ""
    cellPositionHCalBarrelTool = None
    cellPositionHCalEndcapTool = None

# Empty cells for parts of calorimeter not implemented yet
from Configurables import CreateEmptyCaloCellsCollection
createemptycells = CreateEmptyCaloCellsCollection("CreateEmptyCaloCells")
createemptycells.cells.Path = "emptyCaloCells"

if doSWClustering:

    # Produce sliding window clusters
    from Configurables import CaloTowerToolFCCee
    from Configurables import CreateCaloClustersSlidingWindowFCCee

    # Clustering parameters
    # - phi-theta window sizes
    windT = 9
    windP = 17
    posT = 5
    posP = 11
    dupT = 7
    dupP = 13
    finT = 9
    finP = 17
    # - minimal energy to create a cluster in GeV (FCC-ee detectors have to reconstruct low energy particles)
    threshold = 0.040

    # ECAL-only clusters
    ecalBarrelTowers = CaloTowerToolFCCee("CreateECalBarrelTowers",
                                          deltaThetaTower=4 * 0.009817477 / 4, deltaPhiTower=2 * 2 * pi / 1536.,
                                          ecalBarrelReadoutName=ecalBarrelReadoutName,
                                          ecalEndcapReadoutName="",
                                          ecalFwdReadoutName="",
                                          hcalBarrelReadoutName="",
                                          hcalExtBarrelReadoutName="",
                                          hcalEndcapReadoutName="",
                                          hcalFwdReadoutName="",
                                          OutputLevel=INFO)
    ecalBarrelTowers.ecalBarrelCells.Path = ecalBarrelPositionedCellsName
    ecalBarrelTowers.ecalEndcapCells.Path = "emptyCaloCells"
    ecalBarrelTowers.ecalFwdCells.Path = "emptyCaloCells"
    ecalBarrelTowers.hcalBarrelCells.Path = "emptyCaloCells"
    ecalBarrelTowers.hcalExtBarrelCells.Path = "emptyCaloCells"
    ecalBarrelTowers.hcalEndcapCells.Path = "emptyCaloCells"
    ecalBarrelTowers.hcalFwdCells.Path = "emptyCaloCells"

    createECalBarrelClusters = CreateCaloClustersSlidingWindowFCCee("CreateECalBarrelClusters",
                                                                    towerTool=ecalBarrelTowers,
                                                                    nThetaWindow=windT, nPhiWindow=windP,
                                                                    nThetaPosition=posT, nPhiPosition=posP,
                                                                    nThetaDuplicates=dupT, nPhiDuplicates=dupP,
                                                                    nThetaFinal=finT, nPhiFinal=finP,
                                                                    energyThreshold=threshold,
                                                                    energySharingCorrection=False,
                                                                    attachCells=True,
                                                                    OutputLevel=INFO
                                                                    )
    createECalBarrelClusters.clusters.Path = "EMBCaloClusters"
    createECalBarrelClusters.clusterCells.Path = "EMBCaloClusterCells"

    ecalEndcapTowers = CaloTowerToolFCCee("CreateECalEndcapTowers",
                                          deltaThetaTower=4 * 0.009817477 / 4, deltaPhiTower=2 * 2 * pi / 1536.,
                                          ecalBarrelReadoutName="",
                                          ecalEndcapReadoutName=ecalEndcapReadoutName,
                                          ecalFwdReadoutName="",
                                          hcalBarrelReadoutName="",
                                          hcalExtBarrelReadoutName="",
                                          hcalEndcapReadoutName="",
                                          hcalFwdReadoutName="",
                                          OutputLevel=INFO)
    ecalEndcapTowers.ecalBarrelCells.Path = "emptyCaloCells"
    ecalEndcapTowers.ecalEndcapCells.Path = ecalEndcapPositionedCellsName
    ecalEndcapTowers.ecalFwdCells.Path = "emptyCaloCells"
    ecalEndcapTowers.hcalBarrelCells.Path = "emptyCaloCells"
    ecalEndcapTowers.hcalExtBarrelCells.Path = "emptyCaloCells"
    ecalEndcapTowers.hcalEndcapCells.Path = "emptyCaloCells"
    ecalEndcapTowers.hcalFwdCells.Path = "emptyCaloCells"

    createECalEndcapClusters = CreateCaloClustersSlidingWindowFCCee("CreateECalEndcapClusters",
                                                                    towerTool=ecalEndcapTowers,
                                                                    nThetaWindow=windT, nPhiWindow=windP,
                                                                    nThetaPosition=posT, nPhiPosition=posP,
                                                                    nThetaDuplicates=dupT, nPhiDuplicates=dupP,
                                                                    nThetaFinal=finT, nPhiFinal=finP,
                                                                    energyThreshold=threshold,
                                                                    energySharingCorrection=False,
                                                                    attachCells=True,
                                                                    OutputLevel=INFO
                                                                    )
    createECalEndcapClusters.clusters.Path = "EMECCaloClusters"
    createECalEndcapClusters.clusterCells.Path = "EMECCaloClusterCells"

    if applyUpDownstreamCorrections:
        from Configurables import CorrectCaloClusters
        correctECalBarrelClusters = CorrectCaloClusters("CorrectECalBarrelClusters",
                                                        inClusters=createECalBarrelClusters.clusters.Path,
                                                        outClusters="Corrected" + createECalBarrelClusters.clusters.Path,
                                                        systemIDs=[4],
                                                        numLayers=[ecalBarrelLayers],
                                                        firstLayerIDs=[0],
                                                        lastLayerIDs=[ecalBarrelLayers - 1],
                                                        readoutNames=[ecalBarrelReadoutName],
                                                        upstreamParameters=ecalBarrelUpstreamParameters,
                                                        upstreamFormulas=[
                                                            ['[0]+[1]/(x-[2])', '[0]+[1]/(x-[2])']],
                                                        downstreamParameters=ecalBarrelDownstreamParameters,
                                                        downstreamFormulas=[
                                                            ['[0]+[1]*x', '[0]+[1]/sqrt(x)', '[0]+[1]/x']],
                                                        OutputLevel=INFO
                                                        )

    if addShapeParameters:
        from Configurables import AugmentClustersFCCee
        augmentECalBarrelClusters = AugmentClustersFCCee("AugmentECalBarrelClusters",
                                                         inClusters=createECalBarrelClusters.clusters.Path,
                                                         outClusters="Augmented" + createECalBarrelClusters.clusters.Path,
                                                         systemIDs=[4],
                                                         systemNames=["EMB"],
                                                         numLayers=[ecalBarrelLayers],
                                                         readoutNames=[ecalBarrelReadoutName],
                                                         layerFieldNames=["layer"],
                                                         thetaRecalcWeights=[ecalBarrelThetaWeights],
                                                         do_photon_shapeVar=runPhotonIDTool,
                                                         do_widthTheta_logE_weights=logEWeightInPhotonID,
                                                         OutputLevel=INFO
                                                         )

    if applyMVAClusterEnergyCalibration:
        inClusters = ""
        if addShapeParameters:
            inClusters = augmentECalBarrelClusters.outClusters.Path
        else:
            inClusters = createECalBarrelClusters.clusters.Path

        from Configurables import CalibrateCaloClusters
        calibrateECalBarrelClusters = CalibrateCaloClusters("calibrateECalBarrelClusters",
                                                            inClusters=inClusters,
                                                            outClusters="Calibrated" + createECalBarrelClusters.clusters.Path,
                                                            systemIDs=[4],
                                                            systemNames=["EMB"],
                                                            numLayers=[ecalBarrelLayers],
                                                            firstLayerIDs=[0],
                                                            readoutNames=[
                                                                ecalBarrelReadoutName],
                                                            layerFieldNames=["layer"],
                                                            calibrationFile="lgbm_calibration-CaloClusters.onnx",
                                                            OutputLevel=INFO
                                                            )

    if runPhotonIDTool:
        if not addShapeParameters:
            print("Photon ID tool cannot be run if shower shape parameters are not calculated")
            runPhotonIDTool = False
        else:
            inClusters = ""
            if applyMVAClusterEnergyCalibration:
                inClusters = calibrateECalBarrelClusters.outClusters.Path
            else:
                inClusters = augmentECalBarrelClusters.outClusters.Path

            from Configurables import PhotonIDTool
            photonIDECalBarrelClusters = PhotonIDTool("photonIDECalBarrelClusters",
                                                      inClusters=inClusters,
                                                      outClusters="PhotonID" + inClusters,
                                                      mvaModelFile="bdt-photonid-weights-CaloClusters.onnx",
                                                      mvaInputsFile="bdt-photonid-inputs-CaloClusters.json",
                                                      OutputLevel=INFO
                                                      )

    # ECAL + HCAL clusters
    if runHCal:
        towers = CaloTowerToolFCCee("towers",
                                    deltaThetaTower=4 * 0.009817477 / 4, deltaPhiTower=2 * 2 * pi / 1536.,
                                    ecalBarrelReadoutName=ecalBarrelReadoutName,
                                    ecalEndcapReadoutName=ecalEndcapReadoutName,
                                    ecalFwdReadoutName="",
                                    hcalBarrelReadoutName=hcalBarrelReadoutName,
                                    hcalExtBarrelReadoutName="",
                                    hcalEndcapReadoutName=hcalEndcapReadoutName,
                                    hcalFwdReadoutName="",
                                    OutputLevel=INFO)
        towers.ecalBarrelCells.Path = ecalBarrelPositionedCellsName
        towers.ecalEndcapCells.Path = ecalEndcapPositionedCellsName
        towers.ecalFwdCells.Path = "emptyCaloCells"
        towers.hcalBarrelCells.Path = hcalBarrelPositionedCellsName
        towers.hcalExtBarrelCells.Path = "emptyCaloCells"
        towers.hcalEndcapCells.Path = hcalEndcapPositionedCellsName
        towers.hcalFwdCells.Path = "emptyCaloCells"

        createClusters = CreateCaloClustersSlidingWindowFCCee("CreateCaloClusters",
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

        # add here E+H cluster calibration tool or anything else for E+H clusters


if doTopoClustering:

    # Produce topoclusters

    #  ECAL only

    # Neighbours map
    ecalBarrelNeighboursMap = "neighbours_map_ecalB_thetamodulemerged.root"
    from Configurables import TopoCaloNeighbours
    readECalBarrelNeighboursMap = TopoCaloNeighbours("ReadECalBarrelNeighboursMap",
                                                     fileName=ecalBarrelNeighboursMap,
                                                     OutputLevel=INFO)

    # Noise levels per cell
    ecalBarrelNoiseMap = "cellNoise_map_electronicsNoiseLevel_ecalB_thetamodulemerged.root"
    from Configurables import TopoCaloNoisyCells
    readECalBarrelNoisyCellsMap = TopoCaloNoisyCells("ReadECalBarrelNoisyCellsMap",
                                                     fileName=ecalBarrelNoiseMap,
                                                     OutputLevel=INFO)
    
    from Configurables import CaloTopoClusterFCCee
    createECalBarrelTopoClusters = CaloTopoClusterFCCee("CreateECalBarrelTopoClusters",
                                                        cells=[ecalBarrelPositionedCellsName],
                                                        clusters="EMBCaloTopoClusters",
                                                        clusterCells="EMBCaloTopoClusterCells",
                                                        neigboursTool=readECalBarrelNeighboursMap,
                                                        noiseTool=readECalBarrelNoisyCellsMap,
                                                        seedSigma=4,
                                                        neighbourSigma=2,
                                                        lastNeighbourSigma=0,
                                                        minClusterEnergy=0.2,
                                                        OutputLevel=INFO)

    # Neighbours map
    ecalEndcapNeighboursMap = "neighbours_map_ecalE_turbine.root"
    readECalEndcapNeighboursMap = TopoCaloNeighbours("ReadECalEndcapNeighboursMap",
                                                     fileName=ecalEndcapNeighboursMap,
                                                     OutputLevel=INFO)

    # Noise levels per cell
    ecalEndcapNoiseMap = "cellNoise_map_endcapTurbine_electronicsNoiseLevel.root"
    readECalEndcapNoisyCellsMap = TopoCaloNoisyCells("ReadECalEndcapNoisyCellsMap",
                                                     fileName=ecalEndcapNoiseMap,
                                                     OutputLevel=INFO)

    createECalEndcapTopoClusters = CaloTopoClusterFCCee("CreateECalEndcapTopoClusters",
                                                        cells=[ecalEndcapPositionedCellsName],
                                                        clusters="EMECaloTopoClusters",
                                                        clusterCells="EMECaloTopoClusterCells",
                                                        neigboursTool=readECalEndcapNeighboursMap,
                                                        noiseTool=readECalEndcapNoisyCellsMap,
                                                        seedSigma=4,
                                                        neighbourSigma=2,
                                                        lastNeighbourSigma=0,
                                                        OutputLevel=INFO)

    if applyUpDownstreamCorrections:
        from Configurables import CorrectCaloClusters
        correctECalBarrelTopoClusters = CorrectCaloClusters(
            "CorrectECalBarrelTopoClusters",
            inClusters=createECalBarrelTopoClusters.clusters.Path,
            outClusters="Corrected" + createECalBarrelTopoClusters.clusters.Path,
            systemIDs=[4],
            numLayers=[ecalBarrelLayers],
            firstLayerIDs=[0],
            lastLayerIDs=[ecalBarrelLayers - 1],
            readoutNames=[ecalBarrelReadoutName],
            # do not split the following line or it will break scripts that update the values of the corrections
            upstreamParameters=ecalBarrelUpstreamParameters,
            upstreamFormulas=[['[0]+[1]/(x-[2])', '[0]+[1]/(x-[2])']],
            # do not split the following line or it will break scripts that update the values of the corrections
            downstreamParameters=ecalBarrelDownstreamParameters,
            downstreamFormulas=[['[0]+[1]*x', '[0]+[1]/sqrt(x)', '[0]+[1]/x']],
            OutputLevel=INFO
        )

    if addShapeParameters:
        from Configurables import AugmentClustersFCCee
        augmentECalBarrelTopoClusters = AugmentClustersFCCee("augmentECalBarrelTopoClusters",
                                                             inClusters=createECalBarrelTopoClusters.clusters.Path,
                                                             outClusters="Augmented" + createECalBarrelTopoClusters.clusters.Path,
                                                             systemIDs=[4],
                                                             systemNames=["EMB"],
                                                             numLayers=[ecalBarrelLayers],
                                                             readoutNames=[ecalBarrelReadoutName],
                                                             layerFieldNames=["layer"],
                                                             thetaRecalcWeights=[ecalBarrelThetaWeights],
                                                             do_photon_shapeVar=runPhotonIDTool,
                                                             do_widthTheta_logE_weights=logEWeightInPhotonID,
                                                             OutputLevel=INFO)
        if addPi0RecoTool:
            from Configurables import PairCaloClustersPi0
            Pi0RecoAlg = PairCaloClustersPi0(
                "resolvedPi0FromClusterPair",
                inClusters="Augmented" + createECalBarrelTopoClusters.clusters.Path,
                unpairedClusters="Unpaired" + "Augmented" + createECalBarrelTopoClusters.clusters.Path,
                reconstructedPi0="ResolvedPi0Particle",
                massPeak=0.122201, # values determined from a dedicated study
                massLow=0.0754493,
                massHigh=0.153543,
                OutputLevel=INFO
            )

    if applyMVAClusterEnergyCalibration:
        inClusters = ""
        if addShapeParameters:
            inClusters = "Augmented" + createECalBarrelTopoClusters.clusters.Path
        else:
            inClusters = createECalBarrelTopoClusters.clusters.Path

        from Configurables import CalibrateCaloClusters
        calibrateECalBarrelTopoClusters = CalibrateCaloClusters("calibrateECalBarrelTopoClusters",
                                                                inClusters=inClusters,
                                                                outClusters="Calibrated" + createECalBarrelTopoClusters.clusters.Path,
                                                                systemIDs=[4],
                                                                systemNames=["EMB"],
                                                                numLayers=[ecalBarrelLayers],
                                                                firstLayerIDs=[0],
                                                                readoutNames=[
                                                                    ecalBarrelReadoutName],
                                                                layerFieldNames=["layer"],
                                                                calibrationFile="lgbm_calibration-CaloTopoClusters.onnx",
                                                                OutputLevel=INFO
                                                                )

    if runPhotonIDTool:
        if not addShapeParameters:
            print("Photon ID tool cannot be run if shower shape parameters are not calculated")
            runPhotonIDTool = False
        else:
            inClusters = ""
            if applyMVAClusterEnergyCalibration:
                inClusters = calibrateECalBarrelTopoClusters.outClusters.Path
            else:
                inClusters = augmentECalBarrelTopoClusters.outClusters.Path

            from Configurables import PhotonIDTool
            photonIDECalBarrelTopoClusters = PhotonIDTool("photonIDECalBarrelTopoClusters",
                                                          inClusters=inClusters,
                                                          outClusters="PhotonID" + inClusters,
                                                          mvaModelFile="bdt-photonid-weights-CaloTopoClusters.onnx",
                                                          mvaInputsFile="bdt-photonid-inputs-CaloTopoClusters.json",
                                                          OutputLevel=INFO)

    # ECAL + HCAL
    if runHCal:
        # Neighbours map
        neighboursMap = "neighbours_map_ecalB_thetamodulemerged_ecalE_turbine_hcalB_hcalEndcap_phitheta.root"
        readNeighboursMap = TopoCaloNeighbours("ReadNeighboursMap",
                                               fileName=neighboursMap,
                                               OutputLevel=INFO)

        # Noise levels per cell
        noiseMap = "cellNoise_map_electronicsNoiseLevel_ecalB_ECalBarrelModuleThetaMerged_ecalE_ECalEndcapTurbine_hcalB_HCalBarrelReadout_hcalE_HCalEndcapReadout.root"
        readNoisyCellsMap = TopoCaloNoisyCells("ReadNoisyCellsMap",
                                               fileName=noiseMap,
                                               OutputLevel=INFO)

        createTopoClusters = CaloTopoClusterFCCee("CreateTopoClusters",
                                                  cells=[ecalBarrelPositionedCellsName, ecalEndcapPositionedCellsName, hcalBarrelPositionedCellsName, hcalEndcapPositionedCellsName],
                                                  clusters="CaloTopoClusters",
                                                  clusterCells="CaloTopoClusterCells",
                                                  neigboursTool=readNeighboursMap,
                                                  noiseTool=readNoisyCellsMap,
                                                  seedSigma=4,
                                                  neighbourSigma=2,
                                                  lastNeighbourSigma=0,
                                                  OutputLevel=INFO)

# Output
from Configurables import PodioOutput
out = PodioOutput("out",
                  OutputLevel=INFO)
out.filename = "ALLEGRO_sim_digi_reco.root"

out.outputCommands = ["keep *",
                      "drop emptyCaloCells"]

# drop the uncalibrated cells
if dropUncalibratedCells:
    out.outputCommands.append("drop %s" % ecalBarrelReadoutName)
    out.outputCommands.append("drop %s" % ecalBarrelReadoutName2)
    out.outputCommands.append("drop %s" % ecalEndcapReadoutName)
    if runHCal:
        out.outputCommands.append("drop %s" % hcalBarrelReadoutName)
        out.outputCommands.append("drop %s" % hcalEndcapReadoutName)
    else:
        out.outputCommands += ["drop HCal*"]

    # drop the intermediate ecal barrel cells in case of a resegmentation
    if resegmentECalBarrel:
        out.outputCommands.append("drop ECalBarrelCellsMerged")
    # drop the intermediate hcal barrel cells before resegmentation
    if runHCal:
        out.outputCommands.append("drop %s" % hcalBarrelPositionedCellsName)
        out.outputCommands.append("drop %s" % hcalEndcapPositionedCellsName)

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
    if runHCal:
        out.outputCommands.append("drop %s" % hcalBarrelPositionedCellsName)
        out.outputCommands.append("drop %s" % hcalEndcapPositionedCellsName)
if not saveClusterCells:
    out.outputCommands.append("drop Calo*ClusterCells*")

# if we decorate the clusters, we can drop the non-decorated ones
# commented in tests, for debugging
# if addShapeParameters:
#     out.outputCommands.append("drop %s" % augmentECalBarrelClusters.inClusters)

# CPU information
from Configurables import AuditorSvc, ChronoAuditor
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
    createEcalEndcapCells,
]
createEcalBarrelCells.AuditExecute = True
createEcalEndcapCells.AuditExecute = True
if addNoise:
    TopAlg += [
        createEcalBarrelCellsNoise,
        createEcalBarrelCellsNoiseFiltered
    ]
    createEcalBarrelCellsNoise.AuditExecute = True
    createEcalBarrelCellsNoiseFiltered.AuditExecute = True

if resegmentECalBarrel:
    TopAlg += [
        resegmentEcalBarrelTool,
        createEcalBarrelCells2,
    ]
    resegmentEcalBarrelTool.AuditExecute = True
    createEcalBarrelCells2.AuditExecute = True

if runHCal:
    TopAlg += [
        createHCalBarrelCells,
        createHCalEndcapCells,
    ]
    createHCalBarrelCells.AuditExecute = True
    createHCalEndcapCells.AuditExecute = True

if doSWClustering or doTopoClustering:
    TopAlg += [createemptycells]
    createemptycells.AuditExecute = True

    if doSWClustering:
        TopAlg += [createECalBarrelClusters, createECalEndcapClusters]
        createECalBarrelClusters.AuditExecute = True
        createECalEndcapClusters.AuditExecute = True

        if applyUpDownstreamCorrections:
            TopAlg += [correctECalBarrelClusters]
            correctECalBarrelClusters.AuditExecute = True

        if addShapeParameters:
            TopAlg += [augmentECalBarrelClusters]
            augmentECalBarrelClusters.AuditExecute = True

        if applyMVAClusterEnergyCalibration:
            TopAlg += [calibrateECalBarrelClusters]
            calibrateECalBarrelClusters.AuditExecute = True

        if runPhotonIDTool:
            TopAlg += [photonIDECalBarrelClusters]
            photonIDECalBarrelClusters.AuditExecute = True

        if runHCal:
            TopAlg += [createClusters]
            createClusters.AuditExecute = True

    if doTopoClustering:
        TopAlg += [createECalBarrelTopoClusters]
        createECalBarrelTopoClusters.AuditExecute = True

        TopAlg += [createECalEndcapTopoClusters]
        createECalEndcapTopoClusters.AuditExecute = True
        
        if applyUpDownstreamCorrections:
            TopAlg += [correctECalBarrelTopoClusters]
            correctECalBarrelTopoClusters.AuditExecute = True

        if addShapeParameters:
            TopAlg += [augmentECalBarrelTopoClusters]
            augmentECalBarrelTopoClusters.AuditExecute = True

            if addPi0RecoTool:
                TopAlg += [Pi0RecoAlg]
                Pi0RecoAlg.AuditExecute = True

        if applyMVAClusterEnergyCalibration:
            TopAlg += [calibrateECalBarrelTopoClusters]
            calibrateECalBarrelTopoClusters.AuditExecute = True

        if runPhotonIDTool:
            TopAlg += [photonIDECalBarrelTopoClusters]
            photonIDECalBarrelTopoClusters.AuditExecute = True

        if runHCal:
            TopAlg += [createTopoClusters]
            createTopoClusters.AuditExecute = True

TopAlg += [
    out
]

from Configurables import ApplicationMgr
ApplicationMgr(
    TopAlg=TopAlg,
    EvtSel='NONE',
    EvtMax=Nevts,
    ExtSvc=ExtSvc,
    StopOnSignal=True,
)
