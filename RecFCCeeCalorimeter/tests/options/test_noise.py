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
inputfile = "noise_sim.root"               # input file produced with ddsim
outputfile = "noise_digi.root"             # output file to be produced
Nevts = -1                                 # -1 means all events
addNoise = True                            # add noise or not to the cell energy
runHCal = False                            # off since no noise added to HCal so far..

# - what to save in output file
#
# always drop uncalibrated cells, except for tests and debugging
# dropUncalibratedCells = True
dropUncalibratedCells = False

# for big productions, save significant space removing hits and cells
# saveHits = False
# saveCells = False
saveHits = True
saveCells = True

# ECAL barrel parameters for digitisation
ecalBarrelSamplingFraction = [0.3800493723322256] * 1 + [0.13494147915064658] * 1 + [0.142866851721152] * 1 + [0.14839315921940666] * 1 + [0.15298362570665006] * 1 + [0.15709704561942747] * 1 + [0.16063717490147533] * 1 + [0.1641723795419055] * 1 + [0.16845490287689746] * 1 + [0.17111520115997653] * 1 + [0.1730605163148862] * 1
ecalBarrelUpstreamParameters = [[0.028158491043365624, -1.564259408365951, -76.52312805346982, 0.7442903558010191, -34.894692961350195, -74.19340877431723]]
ecalBarrelDownstreamParameters = [[0.00010587711361028165, 0.0052371999097777355, 0.69906696456064, -0.9348243433360095, -0.0364714212117143, 8.360401126995626]]
ecalBarrelLayers = len(ecalBarrelSamplingFraction)
# ECAL endcap parameters for digitisation
ecalEndcapLayers = 98
ecalEndcapSamplingFraction = [0.0897818] * 1+ [0.221318] * 1+ [0.0820002] * 1+ [0.994281] * 1+ [0.0414437] * 1+ [0.1148] * 1+ [0.178831] * 1+ [0.142449] * 1+ [0.181206] * 1+ [0.342843] * 1+ [0.137479] * 1+ [0.176479] * 1+ [0.153273] * 1+ [0.195836] * 1+ [0.0780405] * 1+ [0.150202] * 1+ [0.17846] * 1+ [0.164886] * 1+ [0.175758] * 1+ [0.10836] * 1+ [0.160243] * 1+ [0.183373] * 1+ [0.171818] * 1+ [0.194848] * 1+ [0.111899] * 1+ [0.170704] * 1+ [0.188455] * 1+ [0.178164] * 1+ [0.209113] * 1+ [0.105241] * 1+ [0.180637] * 1+ [0.192206] * 1+ [0.186096] * 1+ [0.211962] * 1+ [0.112019] * 1+ [0.180344] * 1+ [0.195684] * 1+ [0.190778] * 1+ [0.218259] * 1+ [0.118516] * 1+ [0.207786] * 1+ [0.204474] * 1+ [0.207048] * 1+ [0.225913] * 1+ [0.111325] * 1+ [0.147875] * 1+ [0.195625] * 1+ [0.173326] * 1+ [0.175449] * 1+ [0.104087] * 1+ [0.153645] * 1+ [0.161263] * 1+ [0.165499] * 1+ [0.171758] * 1+ [0.175789] * 1+ [0.180657] * 1+ [0.184563] * 1+ [0.187876] * 1+ [0.191762] * 1+ [0.19426] * 1+ [0.197959] * 1+ [0.199021] * 1+ [0.204428] * 1+ [0.195709] * 1+ [0.151751] * 1+ [0.171477] * 1+ [0.165509] * 1+ [0.172565] * 1+ [0.172961] * 1+ [0.175534] * 1+ [0.177989] * 1+ [0.18026] * 1+ [0.181898] * 1+ [0.183912] * 1+ [0.185654] * 1+ [0.187515] * 1+ [0.190408] * 1+ [0.188794] * 1+ [0.193699] * 1+ [0.192287] * 1+ [0.19755] * 1+ [0.190943] * 1+ [0.218553] * 1+ [0.161085] * 1+ [0.373086] * 1+ [0.122495] * 1+ [0.21103] * 1+ [1] * 1+ [0.138686] * 1+ [0.0545171] * 1+ [1] * 1+ [1] * 1+ [0.227945] * 1+ [0.0122872] * 1+ [0.00437334] * 1+ [0.00363533] * 1+ [1] * 1+ [1] * 1


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

# Calorimeter digitisation (merging hits into cells, EM scale calibration via sampling fractions)

# - ECAL readouts
ecalBarrelReadoutName = "ECalBarrelModuleThetaMerged"      # barrel, original segmentation (baseline)
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

from Configurables import CellPositionsECalEndcapTurbineSegTool
cellPositionEcalEndcapTool = CellPositionsECalEndcapTurbineSegTool(
    "CellPositionsECalEndcap",
    readoutName=ecalEndcapReadoutName,
    OutputLevel=INFO
)

# in principle not needed - have to make c++ code of noise tool not fail if position tool is None
# and fail only later when/if tool is used
cellPositionEcalEndcapToolForNoise = CellPositionsECalEndcapTurbineSegTool(
    "CellPositionsECalEndcapForNoise",
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

# - noise tool
if addNoise:
    ecalBarrelNoisePath = "elecNoise_ecalBarrelFCCee_theta.root"
    ecalBarrelNoiseRMSHistName = "h_elecNoise_fcc_"
    from Configurables import NoiseCaloCellsFromFileBarrelTool
    ecalBarrelNoiseTool = NoiseCaloCellsFromFileBarrelTool("ecalBarrelNoiseTool",
                                                           cellPositionsTool=cellPositionEcalBarrelToolForNoise,
                                                           readoutName=ecalBarrelReadoutName,
                                                           noiseFileName=ecalBarrelNoisePath,
                                                           elecNoiseRMSHistoName=ecalBarrelNoiseRMSHistName,
                                                           setNoiseOffset=False,
                                                           activeFieldName="layer",
                                                           addPileup=False,
                                                           filterNoiseThreshold=1,
                                                           useAbsInFilter=True,
                                                           numHistograms=ecalBarrelLayers,
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
    ecalEndcapNoisePath = "elecNoise_ecalendcap.root"
    ecalEndcapNoiseRMSHistName = "noise_endcap_wheel"

    # can enable only when geometry tool for ecal endcap is implemented
    # ecalEndcapGeometryTool = ..
    # from Configurables import NoiseCaloCellsFromFileTurbineEndcapTool
    # ecalEndcapNoiseTool = NoiseCaloCellsFromFileTurbineEndcapTool("ecalEndcapNoiseTool",
    #                                                               cellPositionsTool=cellPositionEcalEndcapToolForNoise,
    #                                                               readoutName=ecalEndcapReadoutName,
    #                                                               noiseFileName=ecalEndcapNoisePath,
    #                                                               elecNoiseRMSHistoName=ecalEndcapNoiseRMSHistName,
    #                                                               setNoiseOffset=False,
    #                                                               activeFieldName="wheel",
    #                                                               addPileup=False,
    #                                                               filterNoiseThreshold=1,
    #                                                               useAbsInFilter=True,
    #                                                               numHistograms=3,  # 3 wheels
    #                                                               scaleFactor=1 / 1000.,  # MeV to GeV
    #                                                               OutputLevel=INFO)
else:
    ecalBarrelNoiseTool = None
    ecalBarrelGeometryTool = None
    ecalEndcapNoiseTool = None

# Create cells in ECal barrel (calibrated and positioned - optionally with xtalk and noise added)
# from uncalibrated cells (+cellID info) from ddsim
ecalBarrelPositionedCellsName = ecalBarrelReadoutName + "Positioned"
ecalBarrelLinks = ecalBarrelPositionedCellsName + "SimCaloHitLinks"
from Configurables import CreatePositionedCaloCells
createEcalBarrelCells = CreatePositionedCaloCells("CreatePositionedECalBarrelCells",
                                                  doCellCalibration=True,
                                                  calibTool=calibEcalBarrel,
                                                  positionsTool=cellPositionEcalBarrelTool,
                                                  addCrosstalk=False,
                                                  crosstalkTool=None,
                                                  addCellNoise=False,
                                                  filterCellNoise=False,
                                                  OutputLevel=INFO,
                                                  hits=ecalBarrelReadoutName,
                                                  cells=ecalBarrelPositionedCellsName,
                                                  links=ecalBarrelLinks
                                                  )
TopAlg += [createEcalBarrelCells]

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

    # can enable only when geometry tool for ecal endcap is implemented
    # # cells with noise not filtered
    # ecalEndcapCellsNoiseLinks = ecalEndcapPositionedCellsName + "WithNoise" + "SimCaloHitLinks"
    # createEcalEndcapCellsNoise = CreatePositionedCaloCells("CreatePositionedECalEndcapCellsWithNoise",
    #                                                        doCellCalibration=True,
    #                                                        calibTool=calibEcalEndcap,
    #                                                        positionsTool=cellPositionEcalEndcapTool,
    #                                                        addCellNoise=True,
    #                                                        filterCellNoise=False,
    #                                                        noiseTool=ecalEndcapNoiseTool,
    #                                                        geometryTool=ecalEndcapGeometryTool,
    #                                                        OutputLevel=INFO,
    #                                                        hits=ecalEndcapReadoutName,
    #                                                        cells=ecalEndcapPositionedCellsName + "WithNoise",
    #                                                        links=ecalEndcapCellsNoiseLinks)
    # TopAlg += [createEcalEndcapCellsNoise]

    # # cells with noise filtered
    # ecalEndcapCellsNoiseFilteredLinks = ecalEndcapPositionedCellsName + "WithNoiseFiltered" + "SimCaloHitLinks"
    # createEcalEndcapCellsNoiseFiltered = CreatePositionedCaloCells("CreatePositionedECalEndcapCellsWithNoiseFiltered",
    #                                                                doCellCalibration=True,
    #                                                                calibTool=calibEcalEndcap,
    #                                                                positionsTool=cellPositionEcalEndcapTool,
    #                                                                addCellNoise=True,
    #                                                                filterCellNoise=True,
    #                                                                noiseTool=ecalEndcapNoiseTool,
    #                                                                geometryTool=ecalEndcapGeometryTool,
    #                                                                OutputLevel=INFO,
    #                                                                hits=ecalEndcapReadoutName,  # uncalibrated & unpositioned cells without noise
    #                                                                cells=ecalEndcapPositionedCellsName + "WithNoiseFiltered",
    #                                                                links=ecalEndcapCellsNoiseFilteredLinks
    #                                                                )
    # TopAlg += [createEcalEndcapCellsNoiseFiltered]

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


# Configure output
io_svc.outputCommands = ["keep *",
                         "drop emptyCaloCells"]

# drop the uncalibrated cells
if dropUncalibratedCells:
    io_svc.outputCommands.append("drop %s" % ecalBarrelReadoutName)
    io_svc.outputCommands.append("drop %s" % ecalEndcapReadoutName)
    if runHCal:
        io_svc.outputCommands.append("drop %s" % hcalBarrelReadoutName)
        io_svc.outputCommands.append("drop %s" % hcalEndcapReadoutName)
    else:
        io_svc.outputCommands += ["drop HCal*"]

# drop lumi, vertex, DCH, Muons (unless want to keep for event display)
io_svc.outputCommands.append("drop Lumi*")
io_svc.outputCommands.append("drop Vertex*")
io_svc.outputCommands.append("drop DriftChamber_simHits*")
io_svc.outputCommands.append("drop MuonTagger*")

# drop hits/positioned cells if desired
if not saveHits:
    io_svc.outputCommands.append("drop *%sContributions" % ecalBarrelReadoutName)
    io_svc.outputCommands.append("drop *%sContributions" % ecalEndcapReadoutName)
    if runHCal:
        io_svc.outputCommands.append("drop *%sContributions" % hcalBarrelReadoutName)
        io_svc.outputCommands.append("drop *%sContributions" % hcalEndcapReadoutName)
if not saveCells:
    io_svc.outputCommands.append("drop %s" % ecalBarrelPositionedCellsName)
    io_svc.outputCommands.append("drop %s" % ecalEndcapPositionedCellsName)
    if runHCal:
        io_svc.outputCommands.append("drop %s" % hcalBarrelPositionedCellsName)
        io_svc.outputCommands.append("drop %s" % hcalEndcapPositionedCellsName)
# drop hits<->cells links if either of the two collections are not saved
if not saveHits or not saveCells:
    io_svc.outputCommands.append("drop *SimCaloHitLinks")


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
