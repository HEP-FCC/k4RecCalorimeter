#
# IMPORTS
#
from Gaudi.Configuration import INFO
from GaudiKernel.PhysicalConstants import pi


#
# SETTINGS
#

# - general settings
#
inputfile = "ALLEGRO_sim.root"             # input file produced with ddsim
outputfile = "ALLEGRO_sim_digi_reco.root"  # output file to be produced
addNoise = False
dumpGDML = False
runHCal = True
# for big productions, save significant space removing hits and cells
# however, hits and cluster cells might be wanted for small productions for detailed event displays
# also, cluster cells are needed for the MVA training
saveHits = False
saveCells = False
saveClusterCells = False

# cluster energy corrections
# simple parametrisations of up/downstream losses
applyUpDownstreamBenchmarkCorrections = True


#
# ALGORITHMS AND SERVICES SETUP
#
TopAlg = []  # alg sequence
ExtSvc = []  # list of external services


# CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
ExtSvc += [audsvc]


# Detector geometry
# prefix all xmls with path_to_detector
# if K4GEO is empty, this should use relative path to working directory
import os
from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc",
                    OutputLevel=INFO
                    # OutputLevel=DEBUG  # set to DEBUG to print dd4hep::DEBUG messages in k4geo C++ drivers
                    )
path_to_detector = os.environ.get("K4GEO", "") + '/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/'
detectors_to_use = [
    'ALLEGRO_o1_v03.xml'
]
# prefix all xmls with path_to_detector
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
    if (constant.get('name') == 'DetID_VXD_Barrel' or
        constant.get('name') == 'DetID_VXD_Disks' or
        constant.get('name') == 'DetID_DCH' or
        constant.get('name') == 'DetID_SiWr_Barrel' or
        constant.get('name') == 'DetID_SiWr_Disks' or
        constant.get('name') == 'DetID_ECAL_Barrel' or
        constant.get('name') == 'DetID_ECAL_Endcap' or
        constant.get('name') == 'DetID_HCAL_Barrel' or
        constant.get('name') == 'DetID_HCAL_Endcap' or
        constant.get('name') == 'DetID_Muon_Barrel'):
        IDs[constant.get("name")[6:]] = int(constant.get('value'))
    if (constant.get('name') == 'DetID_Muon_Endcap_1'):
        IDs[constant.get("name")[6:-2]] = int(constant.get('value'))

# Input/Output handling
from k4FWCore import IOSvc
from Configurables import EventDataSvc
io_svc = IOSvc("IOSvc")
io_svc.Input = inputfile
io_svc.Output = outputfile
ExtSvc += [EventDataSvc("EventDataSvc")]

# from Configurables import GeoToGdmlDumpSvc
if dumpGDML:
    from Configurables import GeoToGdmlDumpSvc
    gdmldumpservice = GeoToGdmlDumpSvc("GeoToGdmlDumpSvc")
# 
# Calorimeter digitisation
# Detector readouts
# ECAL
ecalBarrelReadoutName = "ECalBarrelModuleThetaMerged"
ecalEndcapReadoutName = "ECalEndcapTurbine"
# HCAL
if runHCal:
    hcalBarrelReadoutName = "HCalBarrelReadout"
    hcalEndcapReadoutName = "HCalEndcapReadout"
else:
    hcalBarrelReadoutName = ""
    hcalEndcapReadoutName = ""

# EM scale calibration (sampling fraction)
from Configurables import CalibrateInLayersTool
calibEcalBarrel = CalibrateInLayersTool("CalibrateECalBarrel",
                                        samplingFraction=[0.3864252122990472] * 1 + [0.13597644835735828] * 1 + [0.14520427829645913] * 1 + [0.1510076084632846] * 1 + [0.1552347580991012] * 1 + [0.159694330729184] * 1 + [0.1632954482794191] * 1 + [0.16720711037339814] * 1 + [0.17047749048884808] * 1 + [0.17461698117974286] * 1 + [0.1798984163980135] * 1 + [0.17920355117405806] * 1,
                                        readoutName=ecalBarrelReadoutName,
                                        layerFieldName="layer")


calibEcalEndcap = CalibrateInLayersTool("CalibrateECalEndcap",
                                        samplingFraction = [0.16419] * 1 + [0.192898] * 1 + [0.18783] * 1 + [0.193203] * 1 + [0.193928] * 1 + [0.192286] * 1 + [0.199959] * 1 + [0.200153] * 1 + [0.212635] * 1 + [0.180345] * 1 + [0.18488] * 1 + [0.194762] * 1 + [0.197775] * 1 + [0.200504] * 1 + [0.205555] * 1 + [0.203601] * 1 + [0.210877] * 1 + [0.208376] * 1 + [0.216345] * 1 + [0.201452] * 1 + [0.202134] * 1 + [0.207566] * 1 + [0.208152] * 1 + [0.209889] * 1 + [0.211743] * 1 + [0.213188] * 1 + [0.215864] * 1 + [0.22972] * 1 + [0.192515] * 1 + [0.0103233] * 1,
                                        readoutName=ecalEndcapReadoutName,
                                        layerFieldName="layer")


if runHCal:
    from Configurables import CalibrateCaloHitsTool
    calibHcells = CalibrateCaloHitsTool(
        "CalibrateHCal", invSamplingFraction="30.4")
    calibHcalEndcap = CalibrateCaloHitsTool(
        "CalibrateHCalEndcap", invSamplingFraction="31.7")

# Create cells in ECal barrel
# 1. step - merge hits into cells with theta and module segmentation
# (module is a 'physical' cell i.e. lead + LAr + PCB + LAr +lead)
# 2. step - rewrite the cellId using the merged theta-module segmentation
# (merging several modules and severla theta readout cells).
# Add noise at this step if you derived the noise already assuming merged cells

# read the crosstalk map
from Configurables import ReadCaloCrosstalkMap
readCrosstalkMap = ReadCaloCrosstalkMap("ReadCrosstalkMap",
                                       fileName="https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/xtalk_neighbours_map_ecalB_thetamodulemerged.root",
                                       OutputLevel=INFO)

# Step 1: merge hits into cells according to initial segmentation
ecalBarrelCellsName = "ECalBarrelCells"
from Configurables import CreateCaloCells
createEcalBarrelCells = CreateCaloCells("CreateECalBarrelCells",
                                        doCellCalibration=True,
                                        calibTool=calibEcalBarrel,
                                        crosstalksTool=readCrosstalkMap,
                                        addCrosstalk=False,
                                        addCellNoise=False,
                                        filterCellNoise=False,
                                        addPosition=True,
                                        OutputLevel=INFO,
                                        hits=ecalBarrelReadoutName,
                                        cells=ecalBarrelCellsName)
TopAlg += [createEcalBarrelCells]


# Add to Ecal barrel cells the position information
# (good for physics, all coordinates set properly)
from Configurables import CellPositionsECalBarrelModuleThetaSegTool
cellPositionEcalBarrelTool = CellPositionsECalBarrelModuleThetaSegTool(
    "CellPositionsECalBarrel",
    readoutName=ecalBarrelReadoutName,
    OutputLevel=INFO
)
ecalBarrelPositionedCellsName = "ECalBarrelPositionedCells"
from Configurables import CreateCaloCellPositionsFCCee
createEcalBarrelPositionedCells = CreateCaloCellPositionsFCCee(
    "CreateECalBarrelPositionedCells",
    OutputLevel=INFO
)
TopAlg += [createEcalBarrelPositionedCells]
createEcalBarrelPositionedCells.positionsTool = cellPositionEcalBarrelTool
createEcalBarrelPositionedCells.hits.Path = ecalBarrelCellsName
createEcalBarrelPositionedCells.positionedHits.Path = ecalBarrelPositionedCellsName

# Create cells in ECal endcap
ecalEndcapCellsName = "ECalEndcapCells"
createEcalEndcapCells = CreateCaloCells("CreateEcalEndcapCaloCells",
                                        doCellCalibration=True,
                                        calibTool=calibEcalEndcap,
                                        addCellNoise=False,
                                        filterCellNoise=False,
                                        OutputLevel=INFO,
                                        hits=ecalEndcapReadoutName,
                                        cells=ecalEndcapCellsName)
TopAlg += [createEcalEndcapCells]

# Add to Ecal endcap cells the position information
# (good for physics, all coordinates set properly)
from Configurables import CellPositionsECalEndcapTurbineSegTool
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
TopAlg += [createEcalEndcapPositionedCells]
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
                                            hits=hcalBarrelReadoutName,
                                            cells=hcalBarrelCellsName,
                                            OutputLevel=INFO)
    TopAlg += [createHcalBarrelCells]

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
    TopAlg += [createHcalBarrelPositionedCells]
    createHcalBarrelPositionedCells.positionsTool = cellPositionHcalBarrelTool
    createHcalBarrelPositionedCells.hits.Path = hcalBarrelCellsName
    createHcalBarrelPositionedCells.positionedHits.Path = hcalBarrelPositionedCellsName

    # createHcalEndcapCells = CreateCaloCells("CreateHcalEndcapCaloCells",
    #                                    doCellCalibration=True,
    #                                    calibTool=calibHcalEndcap,
    #                                    addCellNoise=False,
    #                                    filterCellNoise=False,
    #                                    OutputLevel=INFO)
    # createHcalEndcapCells.hits.Path="HCalEndcapHits"
    # createHcalEndcapCells.cells.Path="HCalEndcapCells"
    # TopAlg+=[createHcalEndcapCells]
else:
    hcalBarrelCellsName = "emptyCaloCells"
    hcalBarrelPositionedCellsName = "emptyCaloCells"
    cellPositionHcalBarrelTool = None

# Empty cells for parts of calorimeter not implemented yet
from Configurables import CreateEmptyCaloCellsCollection
createemptycells = CreateEmptyCaloCellsCollection("CreateEmptyCaloCells")
createemptycells.cells.Path = "emptyCaloCells"
TopAlg += [createemptycells]


# Produce sliding window clusters (ECAL+HCAL)
from Configurables import CaloTowerToolFCCee
towers = CaloTowerToolFCCee("towers",
                            #deltaThetaTower=4 * 0.009817477 / 4, deltaPhiTower=2 * 2 * pi / 1536.,
                            deltaThetaTower=0.022180,
                            deltaPhiTower=2 * pi / 256.,
                            cells=[ecalBarrelPositionedCellsName,
                                   ecalEndcapPositionedCellsName,
                                   hcalBarrelPositionedCellsName],
                            calorimeterIDs=[
                                IDs["ECAL_Barrel"],
                                IDs["ECAL_Endcap"],
                                IDs["HCAL_Barrel"],
                                ],
                            OutputLevel=INFO)

# Cluster variables (not optimized)
windT = 18
windP = 34
posT = 10
posP = 22
dupT = 14
dupP = 26
finT = 18
finP = 34
# Minimal energy to create a cluster in GeV (FCC-ee detectors have to reconstruct low energy particles)
threshold = 0.5

from Configurables import CreateCaloClustersSlidingWindowFCCee
createClusters = CreateCaloClustersSlidingWindowFCCee("CreateClusters",
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

from Configurables import CorrectCaloClusters
correctCaloClusters = CorrectCaloClusters("correctCaloClusters",
                                          inClusters=createClusters.clusters.Path,
                                          outClusters="Corrected" + createClusters.clusters.Path,
                                          systemIDs=[4,8],
                                          numLayers=[12,13],
                                          firstLayerIDs=[0,0],
                                          lastLayerIDs=[11,12],
                                          readoutNames=[ecalBarrelReadoutName,hcalBarrelReadoutName],
                                          # do not split the following line or it will break scripts that update the values of the corrections
                                          upstreamParameters = [[0.03900891447361534, -4.322941016402328, -139.1811369546787, 0.498342628339746, -3.3545078429754813, -13.99996971344221],[]],
                                          upstreamFormulas=[['[0]+[1]/(x-[2])', '[0]+[1]/(x-[2])'],[]],
                                          # do not split the following line or it will break scripts that update the values of the corrections
                                          downstreamParameters = [[-0.0027661744480442195, 0.006059143775380306, 0.9788596364251927, -1.4951749409378743, -0.08491999337012696, 16.017621428757778],[]],
                                          downstreamFormulas=[['[0]+[1]*x', '[0]+[1]/sqrt(x)', '[0]+[1]/x'],[]],
                                          ## Approximate parameters for the first energy estimate obtained from the benchmark calib for 100 GeV pion
                                          benchmarkParamsApprox = [1.22, 1, 1.04, -0.0019, 25.69, 0.],
                                          ## below parameters and formulas for two functions with ene switch at 9~GeV, the constant term is set to zero 
                                          #benchmarkParametrization = [[1.0109, 2.4768, 1., 1.0837, -5.278, -1.9938, 0.0006, -0.2715, 45.29, -5674.31, 126.04, 0.],[2.0432, 0.2902, -0.0709, 0.004, 1., 1.2099, -432.2, 13.73, -0.0023, -0.2572, 1.99, -0.968, 0.259, -0.0152, 0.]], 
                                          #benchmarkFormulas = [['[0]+[1]/sqrt(x)', '[0]', '[0]+[1]/(x+[2])', '[0]+[1]/x', '[0]+[1]/(x+[2])', '[0]'],['[0]+[1]*x+[2]*x*x+[3]*x*x*x',  '[0]',  '[0]+[1]/(x+[2])**2', '[0]+[1]/x', '[0]+[1]*x+[2]*x*x+[3]*x*x*x', '[0]']],
                                          ## below parameters and formulas for one function (redireved in Nov 2023)
                                          benchmarkParametrization = [[2.04, 1.0626, 1., -3.5, 1.0182, -0.26, 0.0004, 45.9, -5906.11, -129.27, 0.],[]],
                                          benchmarkFormulas = [['[0]/sqrt(x)+[1]', '[0]', '[0]/x+[1]', '[0]/x+[1]', '[0]+[1]/(x-[2])', '[0]'],[]],
                                          benchmarkEneSwitch = -1.,
                                          upstreamCorr = False,
                                          downstreamCorr = False, 
                                          benchmarkCorr = True, 
                                          OutputLevel=INFO
                                          )
TopAlg += [correctCaloClusters]


# Configure output
io_svc.outputCommands = ["keep *",
                         "drop emptyCaloCells"]

if not runHCal:
    io_svc.outputCommands += ["drop HCal*"]

if not saveCells:
    io_svc.outputCommands.append("drop ECal*Cells*")
if not saveClusterCells:
    io_svc.outputCommands.append("drop *ClusterCells*")
if not saveHits:
    io_svc.outputCommands.append("drop ECal*Hits*")
    io_svc.outputCommands.append("drop HCal*Hits*")


# configure the application
print(TopAlg)
print(ExtSvc)
from k4FWCore import ApplicationMgr
applicationMgr = ApplicationMgr(
    TopAlg=TopAlg,
    EvtSel='NONE',
    EvtMax=-1,
    ExtSvc=ExtSvc,
    StopOnSignal=True,
)

# for algo in applicationMgr.TopAlg:
#     algo.AuditExecute = True
