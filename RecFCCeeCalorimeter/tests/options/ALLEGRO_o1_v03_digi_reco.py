#
# COMMON IMPORTS
#

# Logger
from Gaudi.Configuration import INFO, DEBUG, VERBOSE
# units and physical constants
from GaudiKernel.PhysicalConstants import pi


#
# SETTINGS
#

# - general settings
#
inputfile = "ALLEGRO_sim_ee_z_qq.root"     # input file produced with ddsim
outputfile = "ALLEGRO_sim_digi_reco.root"  # output file to be produced
Nevts = -1                                 # -1 means all events
addNoise = True                            # add noise or not to the cell energy
addCrosstalk = True                        # switch on/off the crosstalk
dumpGDML = False                           # create GDML file of detector model
runHCal = True                             # if false, it will produce only ECAL clusters. if true, it will also produce ECAL+HCAL clusters

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
# ECAL endcap parameters for digitisation
# the turbine endcap has calibration "layers" in the both the z and radial
# directions, for each of the three wheels.  So the total number of layers
# is given by:
#
#   ECalEndcapNumCalibZLayersWheel1*ECalEndcapNumCalibRhoLayersWheel1
#  +ECalEndcapNumCalibZLayersWheel2*ECalEndcapNumCalibRhoLayersWheel2
#  +ECalEndcapNumCalibZLayersWheel3*ECalEndcapNumCalibRhoLayersWheel3
#
# which in the current design is 5*10+1*14+1*34 = 98
# NB these are the sampling fractions calculated directly from the detector simulation with
# 40 GeV photons.  An additional calibration will be required to exactly reproduce the
# true cluster energy on average
ecalEndcapSamplingFraction = [0.17550000000000002] * 1  + [0.1695] * 1  + [0.1655] * 1  + [0.1665] * 1  + [0.1635] * 1  + [0.1765] * 1  + [0.1715] * 1  + [0.1655] * 1  + [0.1565] * 1  + [0.1515] * 1  + [0.1845] * 1  + [0.1775] * 1  + [0.1715] * 1  + [0.1615] * 1  + [0.1495] * 1  + [0.1905] * 1  + [0.1845] * 1  + [0.1785] * 1  + [0.1665] * 1  + [0.1615] * 1  + [0.1965] * 1  + [0.1905] * 1  + [0.1845] * 1  + [0.1765] * 1  + [0.1595] * 1  + [0.2005] * 1  + [0.1955] * 1  + [0.1895] * 1  + [0.1825] * 1  + [0.17350000000000002] * 1  + [0.20450000000000002] * 1  + [0.2005] * 1  + [0.1925] * 1  + [0.1875] * 1  + [0.17250000000000001] * 1  + [0.2115] * 1  + [0.20350000000000001] * 1  + [0.1985] * 1  + [0.1925] * 1  + [0.17450000000000002] * 1  + [0.2125] * 1  + [0.20650000000000002] * 1  + [0.2005] * 1  + [0.1945] * 1  + [0.1825] * 1  + [0.2165] * 1  + [0.2095] * 1  + [0.20350000000000001] * 1  + [0.1955] * 1  + [0.1895] * 1  + [0.1675] * 1  + [0.1695] * 1  + [0.17550000000000002] * 1  + [0.1805] * 1  + [0.1855] * 1  + [0.1895] * 1  + [0.1935] * 1  + [0.1965] * 1  + [0.1995] * 1  + [0.2025] * 1  + [0.20550000000000002] * 1  + [0.20750000000000002] * 1  + [0.2095] * 1  + [0.2105] * 1  + [0.1715] * 1  + [0.17350000000000002] * 1  + [0.17450000000000002] * 1  + [0.1785] * 1  + [0.1805] * 1  + [0.1825] * 1  + [0.1845] * 1  + [0.1875] * 1  + [0.1895] * 1  + [0.1905] * 1  + [0.1925] * 1  + [0.1945] * 1  + [0.1955] * 1  + [0.1975] * 1  + [0.1985] * 1  + [0.2005] * 1  + [0.2015] * 1  + [0.20350000000000001] * 1  + [0.20450000000000002] * 1  + [0.20550000000000002] * 1  + [0.20650000000000002] * 1  + [0.20750000000000002] * 1  + [0.20850000000000002] * 1  + [0.2095] * 1  + [0.2105] * 1  + [0.2115] * 1  + [0.2115] * 1  + [0.2145] * 1  + [0.2145] * 1  + [0.2155] * 1  + [0.2155] * 1  + [0.2165] * 1  + [0.2165] * 1  + [0.2155] * 1

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
TopAlg = []  # alg sequence
ExtSvc = []  # list of external services


# Event counter
from Configurables import EventCounter
eventCounter = EventCounter("EventCounter",
                            OutputLevel=INFO,
                            Frequency=1)
TopAlg += [eventCounter]
# add a message sink service if you want a summary table at the end (not needed..)
# ExtSvc += ["Gaudi::Monitoring::MessageSvcSink"]

# CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
ExtSvc += [audsvc]


# Detector geometry
# prefix all xmls with path_to_detector
# if K4GEO is empty, this should use relative path to working directory
from Configurables import GeoSvc
import os
geoservice = GeoSvc("GeoSvc",
                    OutputLevel=INFO
                    # OutputLevel=DEBUG  # set to DEBUG to print dd4hep::DEBUG messages in k4geo C++ drivers
                    )

path_to_detector = os.environ.get("K4GEO", "") + "/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/"
detectors_to_use = [
    'ALLEGRO_o1_v03.xml'
]
geoservice.detectors = [
    os.path.join(path_to_detector, _det) for _det in detectors_to_use
]
ExtSvc += [geoservice]

# retrieve subdetector IDs
import xml.etree.ElementTree as ET
tree = ET.parse(path_to_detector + 'DectDimensions.xml')
root = tree.getroot()
IDs = {}
for constant in root.find('define').findall('constant'):
    if (
        constant.get('name') == 'DetID_VXD_Barrel'
        or constant.get('name') == 'DetID_VXD_Disks'
        or constant.get('name') == 'DetID_DCH'
        or constant.get('name') == 'DetID_SiWr_Barrel'
        or constant.get('name') == 'DetID_SiWr_Disks'
        or constant.get('name') == 'DetID_ECAL_Barrel'
        or constant.get('name') == 'DetID_ECAL_Endcap'
        or constant.get('name') == 'DetID_HCAL_Barrel'
        or constant.get('name') == 'DetID_HCAL_Endcap'
        or constant.get('name') == 'DetID_Muon_Barrel'
    ):
        IDs[constant.get("name")[6:]] = int(constant.get('value'))
    if (constant.get('name') == 'DetID_Muon_Endcap_1'):
        IDs[constant.get("name")[6:-2]] = int(constant.get('value'))
# debug
print("Subdetector IDs:")
print(IDs)

# Input/Output handling
from k4FWCore import IOSvc
from Configurables import EventDataSvc
io_svc = IOSvc("IOSvc")
io_svc.Input = inputfile
io_svc.Output = outputfile
ExtSvc += [EventDataSvc("EventDataSvc")]

# GDML dump of detector model
if dumpGDML:
    from Configurables import GeoToGdmlDumpSvc
    gdmldumpservice = GeoToGdmlDumpSvc("GeoToGdmlDumpSvc")
    ExtSvc += [gdmldumpservice]

# Calorimeter digitisation (merging hits into cells, EM scale calibration via sampling fractions)

# - ECAL readouts
ecalBarrelReadoutName = "ECalBarrelModuleThetaMerged"      # barrel, original segmentation (baseline)
ecalBarrelReadoutName2 = "ECalBarrelModuleThetaMerged2"    # barrel, after re-segmentation (for optimisation studies)
ecalEndcapReadoutName = "ECalEndcapTurbine"                # endcap, turbine-like (baseline)
# - HCAL readouts
if runHCal:
    hcalBarrelReadoutName = "HCalBarrelReadout"            # barrel, original segmentation (phi-theta)
    # hcalBarrelReadoutName = "HCalBarrelReadoutPhiRow"    # barrel, alternative segmentation (phi-row)
    hcalEndcapReadoutName = "HCalEndcapReadout"            # endcap, original segmentation
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
calibEcalEndcap = CalibrateInLayersTool("CalibrateECalEndcap",
                                        samplingFraction=ecalEndcapSamplingFraction,
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

# - noise tool
if addNoise:
    ecalBarrelNoisePath = "elecNoise_ecalBarrelFCCee_theta.root"
    ecalBarrelNoiseRMSHistName = "h_elecNoise_fcc_"
    from Configurables import NoiseCaloCellsVsThetaFromFileTool
    ecalBarrelNoiseTool = NoiseCaloCellsVsThetaFromFileTool("ecalBarrelNoiseTool",
                                                            cellPositionsTool=cellPositionEcalBarrelToolForNoise,
                                                            readoutName=ecalBarrelReadoutName,
                                                            noiseFileName=ecalBarrelNoisePath,
                                                            elecNoiseRMSHistoName=ecalBarrelNoiseRMSHistName,
                                                            setNoiseOffset=False,
                                                            activeFieldName="layer",
                                                            addPileup=False,
                                                            filterNoiseThreshold=1,
                                                            useAbsInFilter=True,
                                                            numRadialLayers=ecalBarrelLayers,
                                                            scaleFactor=1 / 1000.,  # MeV to GeV
                                                            OutputLevel=INFO)

    from Configurables import TubeLayerModuleThetaCaloTool
    ecalBarrelGeometryTool = TubeLayerModuleThetaCaloTool("ecalBarrelGeometryTool",
                                                          readoutName=ecalBarrelReadoutName,
                                                          activeVolumeName="LAr_sensitive",
                                                          activeFieldName="layer",
                                                          activeVolumesNumber=ecalBarrelLayers,
                                                          fieldNames=["system"],
                                                          fieldValues=[IDs["ECAL_Barrel"]],
                                                          OutputLevel=INFO)
else:
    ecalBarrelNoiseTool = None
    ecalBarrelGeometryTool = None

# Create cells in ECal barrel (calibrated and positioned - optionally with xtalk and noise added)
# from uncalibrated cells (+cellID info) from ddsim
ecalBarrelPositionedCellsName = ecalBarrelReadoutName + "Positioned"
ecalBarrelLinks = ecalBarrelPositionedCellsName + "SimCaloHitLinks"
from Configurables import CreatePositionedCaloCells
createEcalBarrelCells = CreatePositionedCaloCells("CreatePositionedECalBarrelCells",
                                                  doCellCalibration=True,
                                                  calibTool=calibEcalBarrel,
                                                  positionsTool=cellPositionEcalBarrelTool,
                                                  addCrosstalk=addCrosstalk,
                                                  crosstalkTool=readCrosstalkMap,
                                                  addCellNoise=False,
                                                  filterCellNoise=False,
                                                  OutputLevel=INFO,
                                                  hits=ecalBarrelReadoutName,
                                                  cells=ecalBarrelPositionedCellsName,
                                                  links=ecalBarrelLinks
                                                  )
TopAlg += [createEcalBarrelCells]

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
    TopAlg += [
        resegmentEcalBarrelTool,
        createEcalBarrelCells2,
    ]

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
TopAlg += [createEcalEndcapCells]

if addNoise:
    # cells with noise not filtered
    ecalBarrelCellsNoiseLinks = ecalBarrelPositionedCellsName + "WithNoise" + "SimCaloHitLinks"
    createEcalBarrelCellsNoise = CreatePositionedCaloCells("CreatePositionedECalBarrelCellsWithNoise",
                                                           doCellCalibration=True,
                                                           calibTool=calibEcalBarrel,
                                                           positionsTool=cellPositionEcalBarrelTool,
                                                           addCellNoise=True,
                                                           filterCellNoise=False,
                                                           noiseTool=ecalBarrelNoiseTool,
                                                           geometryTool=ecalBarrelGeometryTool,
                                                           OutputLevel=INFO,
                                                           hits=ecalBarrelReadoutName,
                                                           cells=ecalBarrelPositionedCellsName + "WithNoise",
                                                           links=ecalBarrelCellsNoiseLinks)
    TopAlg += [createEcalBarrelCellsNoise]

    # cells with noise filtered
    ecalBarrelCellsNoiseFilteredLinks = ecalBarrelPositionedCellsName + "WithNoiseFiltered" + "SimCaloHitLinks"
    createEcalBarrelCellsNoiseFiltered = CreatePositionedCaloCells("CreatePositionedECalBarrelCellsWithNoiseFiltered",
                                                                   doCellCalibration=True,
                                                                   calibTool=calibEcalBarrel,
                                                                   positionsTool=cellPositionEcalBarrelTool,
                                                                   addCellNoise=True,
                                                                   filterCellNoise=True,
                                                                   noiseTool=ecalBarrelNoiseTool,
                                                                   geometryTool=ecalBarrelGeometryTool,
                                                                   OutputLevel=INFO,
                                                                   hits=ecalBarrelReadoutName,  # uncalibrated & unpositioned cells without noise
                                                                   cells=ecalBarrelPositionedCellsName + "WithNoiseFiltered",
                                                                   links=ecalBarrelCellsNoiseFilteredLinks
                                                                   )
    TopAlg += [createEcalBarrelCellsNoiseFiltered]

if runHCal:
    # Apply calibration and positioning to cells in HCal barrel
    hcalBarrelPositionedCellsName = hcalBarrelReadoutName + "Positioned"
    hcalBarrelLinks = hcalBarrelPositionedCellsName + "SimCaloHitLinks"
    createHCalBarrelCells = CreatePositionedCaloCells("CreatePositionedHCalBarrelCells",
                                                      doCellCalibration=True,
                                                      calibTool=calibHCalBarrel,
                                                      positionsTool=cellPositionHCalBarrelTool,
                                                      addCellNoise=False,
                                                      filterCellNoise=False,
                                                      hits=hcalBarrelReadoutName,
                                                      cells=hcalBarrelPositionedCellsName,
                                                      links=hcalBarrelLinks,
                                                      OutputLevel=INFO)
    TopAlg += [createHCalBarrelCells]

    # Apply calibration and positioning to cells in HCal endcap
    hcalEndcapPositionedCellsName = hcalEndcapReadoutName + "Positioned"
    hcalEndcapLinks = hcalEndcapPositionedCellsName + "SimCaloHitLinks"
    createHCalEndcapCells = CreatePositionedCaloCells("CreatePositionedHCalEndcapCells",
                                                      doCellCalibration=True,
                                                      calibTool=calibHCalEndcap,
                                                      addCellNoise=False,
                                                      filterCellNoise=False,
                                                      positionsTool=cellPositionHCalEndcapTool,
                                                      OutputLevel=INFO,
                                                      hits=hcalEndcapReadoutName,
                                                      cells=hcalEndcapPositionedCellsName,
                                                      links=hcalEndcapLinks)
    TopAlg += [createHCalEndcapCells]

finalClusters = []
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
    cells = [ecalBarrelPositionedCellsName]
    caloIDs = [IDs["ECAL_Barrel"]]
    ecalBarrelTowers = CaloTowerToolFCCee("CreateECalBarrelTowers",
                                          deltaThetaTower=4 * 0.009817477 / 4, deltaPhiTower=2 * 2 * pi / 1536.,
                                          thetaMin=pi-2.55254,
                                          thetaMax=2.55254,
                                          phiMin=-pi,
                                          phiMax=pi,
                                          cells=cells,
                                          calorimeterIDs=caloIDs,
                                          nSubDetectors=3,  # just for test here, since there is only the ECAL..
                                          OutputLevel=INFO)

    createECalBarrelClusters = CreateCaloClustersSlidingWindowFCCee("CreateECalBarrelClusters",
                                                                    towerTool=ecalBarrelTowers,
                                                                    nThetaWindow=windT, nPhiWindow=windP,
                                                                    nThetaPosition=posT, nPhiPosition=posP,
                                                                    nThetaDuplicates=dupT, nPhiDuplicates=dupP,
                                                                    nThetaFinal=finT, nPhiFinal=finP,
                                                                    energyThreshold=threshold,
                                                                    energySharingCorrection=False,
                                                                    createClusterCellCollection=True,
                                                                    OutputLevel=INFO
                                                                    )
    createECalBarrelClusters.clusters.Path = "EMBCaloClusters"
    createECalBarrelClusters.clusterCells.Path = "EMBCaloClusterCells"
    TopAlg += [createECalBarrelClusters]
    finalClusters += ["EMBCaloClusters"]

    cells = [ecalEndcapPositionedCellsName]
    caloIDs = [IDs["ECAL_Endcap"]]
    ecalEndcapTowers = CaloTowerToolFCCee("CreateECalEndcapTowers",
                                          deltaThetaTower=4 * 0.009817477 / 4, deltaPhiTower=2 * 2 * pi / 1536.,
                                          nSubDetectors=0,
                                          cells=cells,
                                          calorimeterIDs=caloIDs,
                                          OutputLevel=INFO)

    createECalEndcapClusters = CreateCaloClustersSlidingWindowFCCee("CreateECalEndcapClusters",
                                                                    towerTool=ecalEndcapTowers,
                                                                    nThetaWindow=windT, nPhiWindow=windP,
                                                                    nThetaPosition=posT, nPhiPosition=posP,
                                                                    nThetaDuplicates=dupT, nPhiDuplicates=dupP,
                                                                    nThetaFinal=finT, nPhiFinal=finP,
                                                                    energyThreshold=threshold,
                                                                    energySharingCorrection=False,
                                                                    createClusterCellCollection=True,
                                                                    OutputLevel=INFO
                                                                    )
    createECalEndcapClusters.clusters.Path = "EMECCaloClusters"
    createECalEndcapClusters.clusterCells.Path = "EMECCaloClusterCells"
    TopAlg += [createECalEndcapClusters]
    finalClusters += ["EMECCaloClusters"]

    if applyUpDownstreamCorrections:
        from Configurables import CorrectCaloClusters
        inClusters = createECalBarrelClusters.clusters.Path
        outClusters = "Corrected" + createECalBarrelClusters.clusters.Path
        correctECalBarrelClusters = CorrectCaloClusters("CorrectECalBarrelClusters",
                                                        inClusters=inClusters,
                                                        outClusters=outClusters,
                                                        systemIDs=[IDs["ECAL_Barrel"]],
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
        TopAlg += [correctECalBarrelClusters]

    if addShapeParameters:
        from Configurables import AugmentClustersFCCee
        inClusters=createECalBarrelClusters.clusters.Path
        outClusters="Augmented" + createECalBarrelClusters.clusters.Path
        augmentECalBarrelClusters = AugmentClustersFCCee("AugmentECalBarrelClusters",
                                                         inClusters=inClusters,
                                                         outClusters=outClusters,
                                                         systemIDs=[IDs["ECAL_Barrel"]],
                                                         systemNames=["EMB"],
                                                         numLayers=[ecalBarrelLayers],
                                                         readoutNames=[ecalBarrelReadoutName],
                                                         layerFieldNames=["layer"],
                                                         thetaRecalcWeights=[ecalBarrelThetaWeights],
                                                         do_photon_shapeVar=runPhotonIDTool,
                                                         do_widthTheta_logE_weights=logEWeightInPhotonID,
                                                         OutputLevel=INFO
                                                         )
        TopAlg += [augmentECalBarrelClusters]
        finalClusters.remove(inClusters)
        finalClusters += [outClusters]

    if applyMVAClusterEnergyCalibration:
        inClusters = ""
        if addShapeParameters:
            inClusters = augmentECalBarrelClusters.outClusters.Path
        else:
            inClusters = createECalBarrelClusters.clusters.Path
        outClusters="Calibrated" + createECalBarrelClusters.clusters.Path

        from Configurables import CalibrateCaloClusters
        calibrateECalBarrelClusters = CalibrateCaloClusters("calibrateECalBarrelClusters",
                                                            inClusters=inClusters,
                                                            outClusters=outClusters,
                                                            systemIDs=[IDs["ECAL_Barrel"]],
                                                            systemNames=["EMB"],
                                                            numLayers=[ecalBarrelLayers],
                                                            firstLayerIDs=[0],
                                                            readoutNames=[
                                                                ecalBarrelReadoutName],
                                                            layerFieldNames=["layer"],
                                                            calibrationFile="lgbm_calibration-CaloClusters.onnx",
                                                            OutputLevel=INFO
                                                            )
        TopAlg += [calibrateECalBarrelClusters]
        finalClusters.remove(inClusters)
        finalClusters += [outClusters]

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
            outClusters="PhotonID" + inClusters;

            from Configurables import PhotonIDTool
            photonIDECalBarrelClusters = PhotonIDTool("photonIDECalBarrelClusters",
                                                      inClusters=inClusters,
                                                      outClusters=outClusters,
                                                      mvaModelFile="bdt-photonid-weights-CaloClusters.onnx",
                                                      mvaInputsFile="bdt-photonid-inputs-CaloClusters.json",
                                                      OutputLevel=INFO
                                                      )
            TopAlg += [photonIDECalBarrelClusters]
            finalClusters.remove(inClusters)
            finalClusters += [outClusters]

    # ECAL + HCAL clusters
    if runHCal:
        cells = [ecalBarrelPositionedCellsName,
                 ecalEndcapPositionedCellsName,
                 hcalBarrelPositionedCellsName,
                 hcalEndcapPositionedCellsName]
        caloIDs = [IDs["ECAL_Barrel"],
                   IDs["ECAL_Endcap"],
                   IDs["HCAL_Barrel"],
                   IDs["HCAL_Endcap"]]

        towers = CaloTowerToolFCCee("towers",
                                    deltaThetaTower=4 * 0.009817477 / 4, deltaPhiTower=2 * 2 * pi / 1536.,
                                    cells=cells,
                                    calorimeterIDs=caloIDs,
                                    nSubDetectors = 3,
                                    OutputLevel=INFO)
        createClusters = CreateCaloClustersSlidingWindowFCCee("CreateCaloClusters",
                                                              towerTool=towers,
                                                              nThetaWindow=windT, nPhiWindow=windP,
                                                              nThetaPosition=posT, nPhiPosition=posP,
                                                              nThetaDuplicates=dupT, nPhiDuplicates=dupP,
                                                              nThetaFinal=finT, nPhiFinal=finP,
                                                              energyThreshold=threshold,
                                                              energySharingCorrection=False,
                                                              createClusterCellCollection=True,
                                                              OutputLevel=INFO
                                                              )
        createClusters.clusters.Path = "CaloClusters"
        createClusters.clusterCells.Path = "CaloClusterCells"
        TopAlg += [createClusters]
        finalClusters += ["CaloClusters"]

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
    caloIDs = [IDs["ECAL_Barrel"]]
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
                                                        calorimeterIDs=caloIDs,
                                                        createClusterCellCollection=True,
                                                        OutputLevel=INFO)
    TopAlg += [createECalBarrelTopoClusters]
    finalClusters += ["EMBCaloTopoClusters"]

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

    caloIDs = [IDs["ECAL_Endcap"]]
    createECalEndcapTopoClusters = CaloTopoClusterFCCee("CreateECalEndcapTopoClusters",
                                                        cells=[ecalEndcapPositionedCellsName],
                                                        clusters="EMECCaloTopoClusters",
                                                        clusterCells="EMECCaloTopoClusterCells",
                                                        neigboursTool=readECalEndcapNeighboursMap,
                                                        noiseTool=readECalEndcapNoisyCellsMap,
                                                        seedSigma=4,
                                                        neighbourSigma=2,
                                                        lastNeighbourSigma=0,
                                                        calorimeterIDs=caloIDs,
                                                        createClusterCellCollection=True,
                                                        OutputLevel=INFO)
    TopAlg += [createECalEndcapTopoClusters]
    finalClusters += ["EMECCaloTopoClusters"]

    if applyUpDownstreamCorrections:
        from Configurables import CorrectCaloClusters
        correctECalBarrelTopoClusters = CorrectCaloClusters(
            "CorrectECalBarrelTopoClusters",
            inClusters=createECalBarrelTopoClusters.clusters.Path,
            outClusters="Corrected" + createECalBarrelTopoClusters.clusters.Path,
            systemIDs=[IDs["ECAL_Barrel"]],
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
        TopAlg += [correctECalBarrelTopoClusters]

    if addShapeParameters:
        from Configurables import AugmentClustersFCCee
        inClusters=createECalBarrelTopoClusters.clusters.Path
        outClusters="Augmented" + createECalBarrelTopoClusters.clusters.Path
        augmentECalBarrelTopoClusters = AugmentClustersFCCee("augmentECalBarrelTopoClusters",
                                                             inClusters=inClusters,
                                                             outClusters=outClusters,
                                                             systemIDs=[IDs["ECAL_Barrel"]],
                                                             systemNames=["EMB"],
                                                             numLayers=[ecalBarrelLayers],
                                                             readoutNames=[ecalBarrelReadoutName],
                                                             layerFieldNames=["layer"],
                                                             thetaRecalcWeights=[ecalBarrelThetaWeights],
                                                             do_photon_shapeVar=runPhotonIDTool,
                                                             do_widthTheta_logE_weights=logEWeightInPhotonID,
                                                             OutputLevel=INFO)
        TopAlg += [augmentECalBarrelTopoClusters]
        finalClusters.remove(inClusters)
        finalClusters += [outClusters]

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
            TopAlg += [Pi0RecoAlg]

    if applyMVAClusterEnergyCalibration:
        inClusters = ""
        if addShapeParameters:
            inClusters = "Augmented" + createECalBarrelTopoClusters.clusters.Path
        else:
            inClusters = createECalBarrelTopoClusters.clusters.Path
        outClusters="Calibrated" + createECalBarrelTopoClusters.clusters.Path

        from Configurables import CalibrateCaloClusters
        calibrateECalBarrelTopoClusters = CalibrateCaloClusters("calibrateECalBarrelTopoClusters",
                                                                inClusters=inClusters,
                                                                outClusters=outClusters,
                                                                systemIDs=[IDs["ECAL_Barrel"]],
                                                                systemNames=["EMB"],
                                                                numLayers=[ecalBarrelLayers],
                                                                firstLayerIDs=[0],
                                                                readoutNames=[
                                                                    ecalBarrelReadoutName],
                                                                layerFieldNames=["layer"],
                                                                calibrationFile="lgbm_calibration-CaloTopoClusters.onnx",
                                                                OutputLevel=INFO
                                                                )
        TopAlg += [calibrateECalBarrelTopoClusters]
        finalClusters.remove(inClusters)
        finalClusters += [outClusters]

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
            outClusters="PhotonID" + inClusters

            from Configurables import PhotonIDTool
            photonIDECalBarrelTopoClusters = PhotonIDTool("photonIDECalBarrelTopoClusters",
                                                          inClusters=inClusters,
                                                          outClusters=outClusters,
                                                          mvaModelFile="bdt-photonid-weights-CaloTopoClusters.onnx",
                                                          mvaInputsFile="bdt-photonid-inputs-CaloTopoClusters.json",
                                                          OutputLevel=INFO)
            TopAlg += [photonIDECalBarrelTopoClusters]
            finalClusters.remove(inClusters)
            finalClusters += [outClusters]

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
        TopAlg += [createTopoClusters]
        finalClusters += ["CaloTopoClusters"]

# Create CaloHit<->MCParticle and Cluster<->MCParticle links
from Configurables import CreateTruthLinks
caloLinks = [ecalBarrelLinks, ecalEndcapLinks]
if runHCal:
    caloLinks += [hcalBarrelLinks, hcalEndcapLinks]
createTruthLinks = CreateTruthLinks("CreateTruthLinks",
                                    cell_hit_links=caloLinks,
                                    clusters=finalClusters,
                                    mcparticles="MCParticles",
                                    cell_mcparticle_links="CaloHitMCParticleLinks",
                                    cluster_mcparticle_links="ClusterMCParticleLinks",
                                    OutputLevel=INFO)
TopAlg += [createTruthLinks]


# Configure output
io_svc.outputCommands = ["keep *",
                         "drop emptyCaloCells"]

# drop the uncalibrated cells
if dropUncalibratedCells:
    io_svc.outputCommands.append("drop %s" % ecalBarrelReadoutName)
    io_svc.outputCommands.append("drop %s" % ecalBarrelReadoutName2)
    io_svc.outputCommands.append("drop %s" % ecalEndcapReadoutName)
    if runHCal:
        io_svc.outputCommands.append("drop %s" % hcalBarrelReadoutName)
        io_svc.outputCommands.append("drop %s" % hcalEndcapReadoutName)
    else:
        io_svc.outputCommands += ["drop HCal*"]

    # drop the intermediate ecal barrel cells in case of a resegmentation
    if resegmentECalBarrel:
        io_svc.outputCommands.append("drop ECalBarrelCellsMerged")

# drop lumi, vertex, DCH, Muons (unless want to keep for event display)
io_svc.outputCommands.append("drop Lumi*")
# io_svc.outputCommands.append("drop Vertex*")
# io_svc.outputCommands.append("drop DriftChamber_simHits*")
io_svc.outputCommands.append("drop MuonTagger*")

# drop hits/positioned cells/cluster cells if desired
if not saveHits:
    io_svc.outputCommands.append("drop *%sContributions" % ecalBarrelReadoutName)
    io_svc.outputCommands.append("drop *%sContributions" % ecalBarrelReadoutName2)
    io_svc.outputCommands.append("drop *%sContributions" % ecalEndcapReadoutName)
    if runHCal:
        io_svc.outputCommands.append("drop *%sContributions" % hcalBarrelReadoutName)
        io_svc.outputCommands.append("drop *%sContributions" % hcalEndcapReadoutName)
if not saveCells:
    io_svc.outputCommands.append("drop %s" % ecalBarrelPositionedCellsName)
    io_svc.outputCommands.append("drop %s" % ecalEndcapPositionedCellsName)
    if resegmentECalBarrel:
        io_svc.outputCommands.append("drop %s" % ecalBarrelPositionedCellsName2)
    if runHCal:
        io_svc.outputCommands.append("drop %s" % hcalBarrelPositionedCellsName)
        io_svc.outputCommands.append("drop %s" % hcalEndcapPositionedCellsName)
if not saveClusterCells:
    io_svc.outputCommands.append("drop *Calo*Cluster*Cells*")
# drop hits<->cells links if either of the two collections are not saved
if not saveHits or not saveCells:
    io_svc.outputCommands.append("drop *SimCaloHitLinks")

# if we decorate the clusters, we can drop the non-decorated ones
# commented in tests, for debugging
# if addShapeParameters:
#     io_svc.outputCommands.append("drop %s" % augmentECalBarrelClusters.inClusters)


# configure the application
print(TopAlg)
print(ExtSvc)
from k4FWCore import ApplicationMgr
applicationMgr = ApplicationMgr(
    TopAlg=TopAlg,
    EvtSel='NONE',
    EvtMax=Nevts,
    ExtSvc=ExtSvc,
    StopOnSignal=True,
)

for algo in applicationMgr.TopAlg:
    algo.AuditExecute = True
    # for debug
    # algo.OutputLevel = DEBUG
