from Gaudi.Configuration import *

from Configurables import ApplicationMgr, k4DataSvc, PodioOutput

podioevent = k4DataSvc("EventDataSvc", input="output_ecalSim_e50GeV_1events.root")

# reads HepMC text file and write the HepMC::GenEvent to the data service
from Configurables import PodioInput
podioinput = PodioInput("PodioReader", collections = ["ECalHits", "ECalPositionedHits"], OutputLevel = DEBUG)

from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc", detectors = [  'file:Detector/DetFCChhBaseline1/compact/FCChh_DectEmptyMaster.xml',
                                             'file:Detector/DetFCChhECalSimple/compact/FCChh_ECalBarrel_Mockup.xml'],
                    OutputLevel = DEBUG)

# common ECAL specific information
# readout name
ecalReadoutName = "ECalHitsPhiEta"
# active material identifier name
ecalIdentifierName = "active_layer"
# active material volume name
ecalVolumeName = "LAr_sensitive"
# number of active layers to be merged to create cells
ecalNumberOfLayersToMerge = [19,71,9]
# number of ECAL layers
ecalNumberOfLayers = len(ecalNumberOfLayersToMerge)
# ECAL bitfield names & values
ecalFieldNames = ["system","ECAL_Cryo","bath","EM_barrel"]
ecalFieldValues = [5,1,1,1]

from Configurables import MergeLayers
mergelayers = MergeLayers("MergeLayers",
                   # take the bitfield description from the geometry service
                   readout = ecalReadoutName,
                   # cells in which field should be merged
                   identifier = ecalIdentifierName,
                   volumeName = ecalVolumeName,
                   # how many cells to merge
                   # merge first 19 into new cell (id=0), next 71 into second cell (id=1), ...
                   merge = ecalNumberOfLayersToMerge,
                   OutputLevel = INFO)
mergelayers.inhits.Path = "ECalHits"
mergelayers.outhits.Path = "mergedECalHits"

#Configure tools for calo cell positions                    
from Configurables import CellPositionsDummyTool
dummyCells = CellPositionsDummyTool("CellPositions",
                                         OutputLevel = INFO)
#Configure tools for calo reconstruction
from Configurables import CalibrateCaloHitsTool, NoiseCaloCellsFromFileTool, TubeLayerPhiEtaCaloTool
calibcells = CalibrateCaloHitsTool("CalibrateCaloHitsTool", invSamplingFraction="5.4")
noise = NoiseCaloCellsFromFileTool("NoiseCaloCellsFromFileTool",
                                   cellPositionsTool=dummyCells,
    noiseFileName="http://fccsw.web.cern.ch/fccsw/testsamples/fcc_ecalCellNoise_mu1000_3radialLayers_v1.root")

ecalgeo = TubeLayerPhiEtaCaloTool("EcalGeo",
                                  readoutName = ecalReadoutName,
                                  activeVolumeName = ecalVolumeName,
                                  activeFieldName = ecalIdentifierName,
                                  fieldNames = ecalFieldNames,
                                  fieldValues = ecalFieldValues,
                                  # to make it working with MergeLayers algorithm
                                  activeVolumesNumber = ecalNumberOfLayers,
                                  )

from Configurables import CreateCaloCells
createcells = CreateCaloCells("CreateCaloCells",
                              geometryTool = ecalgeo,
                              doCellCalibration = True,
                              calibTool = calibcells,
                              addCellNoise = True, filterCellNoise = False,
                              noiseTool = noise,
                              OutputLevel=DEBUG,
                              )
createcells.hits.Path="mergedECalHits"
createcells.cells.Path="ecalBarrelCells"

from Configurables import CreateEmptyCaloCellsCollection
createemptycells = CreateEmptyCaloCellsCollection("CreateEmptyCaloCells")
createemptycells.cells.Path = "emptyCaloCells"

#Create calo clusters
from Configurables import CreateCaloClustersSlidingWindow, CaloTowerTool
from GaudiKernel.PhysicalConstants import pi
towers = CaloTowerTool("towers",
                               deltaEtaTower = 0.01, deltaPhiTower = 2*pi/629.,
                               ecalBarrelReadoutName = ecalReadoutName,
                               ecalEndcapReadoutName = "",
                               ecalFwdReadoutName = "",
                               hcalBarrelReadoutName = "",
                               hcalExtBarrelReadoutName = "",
                               hcalEndcapReadoutName = "",
                               hcalFwdReadoutName = "",
                               )
towers.ecalBarrelCells.Path = "ecalBarrelCells"
towers.ecalEndcapCells.Path = "emptyCaloCells"
towers.ecalFwdCells.Path = "emptyCaloCells"
towers.hcalBarrelCells.Path = "emptyCaloCells"
towers.hcalExtBarrelCells.Path = "emptyCaloCells"
towers.hcalEndcapCells.Path = "emptyCaloCells"
towers.hcalFwdCells.Path = "emptyCaloCells"

createclusters = CreateCaloClustersSlidingWindow("CreateCaloClusters",
                                                 towerTool = towers,
                                                 nEtaWindow = 7, nPhiWindow = 15,
                                                 nEtaPosition = 5, nPhiPosition = 11,
                                                 nEtaDuplicates = 5, nPhiDuplicates = 11,
                                                 nEtaFinal = 7, nPhiFinal = 15,
                                                 energyThreshold = 8,
                                                 )
createclusters.clusters.Path = "caloClusters"

out = PodioOutput("output", filename = "output_ecalReco_noiseFromFile_test.root",
                   OutputLevel = DEBUG)
out.outputCommands = ["keep *"]

#CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
podioinput.AuditExecute = True
createclusters.AuditExecute = True
createcells.AuditExecute = True
mergelayers.AuditExecute = True
out.AuditExecute = True

ApplicationMgr(
    TopAlg = [podioinput,
              mergelayers,
              createcells,
              createemptycells,
              createclusters,
              out
              ],
    EvtSel = 'NONE',
    EvtMax = 1,
    ExtSvc = [podioevent, geoservice],
 )
