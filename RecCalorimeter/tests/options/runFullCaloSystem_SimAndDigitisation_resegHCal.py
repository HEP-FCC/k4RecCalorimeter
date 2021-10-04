import os

from GaudiKernel.SystemOfUnits import MeV,GeV

# simulations setup
energy=50*GeV
num_events=3
magnetic_field=1
particleType="pi-"

from Gaudi.Configuration import *

from Configurables import ApplicationMgr, FCCDataSvc, PodioOutput

podioevent   = FCCDataSvc("EventDataSvc")

from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc", detectors=[  'file:Detector/DetFCChhBaseline1/compact/FCChh_DectEmptyMaster.xml',
                                           # Tracker disabled to save cpu time
                                           #'file:Detector/DetFCChhTrackerTkLayout/compact/Tracker.xml',
                                           'file:Detector/DetFCChhECalInclined/compact/FCChh_ECalBarrel_withCryostat.xml',
                                           'file:Detector/DetFCChhHCalTile/compact/FCChh_HCalBarrel_TileCal.xml',
                                           'file:Detector/DetFCChhHCalTile/compact/FCChh_HCalExtendedBarrel_TileCal.xml',
                                           'file:Detector/DetFCChhCalDiscs/compact/Endcaps_coneCryo.xml',
                                           'file:Detector/DetFCChhCalDiscs/compact/Forward_coneCryo.xml'
                                           ],
                    OutputLevel = INFO)

# Geant4 service
# Configures the Geant simulation: geometry, physics list and user actions
from Configurables import SimG4Svc
geantservice = SimG4Svc("SimG4Svc", detector='SimG4DD4hepDetector', physicslist="SimG4FtfpBert", actions="SimG4FullSimActions")

# range cut
geantservice.g4PreInitCommands += ["/run/setCut 0.1 mm"]

# Magnetic field
from Configurables import SimG4ConstantMagneticFieldTool
if magnetic_field==1:
    field = SimG4ConstantMagneticFieldTool("SimG4ConstantMagneticFieldTool",FieldOn=True,IntegratorStepper="ClassicalRK4")
else:
    field = SimG4ConstantMagneticFieldTool("SimG4ConstantMagneticFieldTool",FieldOn=False)

# Geant4 algorithm
# Translates EDM to G4Event, passes the event to G4, writes out outputs via tools
# and a tool that saves the calorimeter hits

# ECAL readouts
ecalBarrelReadoutName = "ECalBarrelEta"
ecalBarrelReadoutNamePhiEta = "ECalBarrelPhiEta"
ecalEndcapReadoutName = "EMECPhiEta"
ecalFwdReadoutName = "EMFwdPhiEta"
# HCAL readouts
hcalBarrelReadoutName = "HCalBarrelReadout"
hcalBarrelReadoutNamePhiEta = "BarHCal_Readout_phieta"
extHcalReadoutName = "HCalExtBarrelReadout"
hcalEndcapReadoutName = "HECPhiEta"
hcalFwdReadoutName = "HFwdPhiEta"
# layers to be merged in endcaps & forward calo
ecalEndcapNumberOfLayersToMerge = [26]*5+[27]
ecalFwdNumberOfLayersToMerge = [7]*5+[8]
hcalEndcapNumberOfLayersToMerge = [13]+[14]*5
hcalFwdNumberOfLayersToMerge = [8]+[9]*5
identifierName = "layer"
volumeName = "layer"

from Configurables import SimG4Alg, SimG4SaveCalHits
saveecalbarreltool = SimG4SaveCalHits("saveECalBarrelHits", readoutNames = [ecalBarrelReadoutName])
saveecalbarreltool.positionedCaloHits.Path = "ECalBarrelPositionedHits"
saveecalbarreltool.caloHits.Path = "ECalBarrelHits"
saveecalendcaptool = SimG4SaveCalHits("saveECalEndcapHits", readoutNames = [ecalEndcapReadoutName])
saveecalendcaptool.positionedCaloHits.Path = "ECalEndcapPositionedHits"
saveecalendcaptool.caloHits.Path = "ECalEndcapHits"
saveecalfwdtool = SimG4SaveCalHits("saveECalFwdHits", readoutNames = [ecalFwdReadoutName])
saveecalfwdtool.positionedCaloHits.Path = "ECalFwdPositionedHits"
saveecalfwdtool.caloHits.Path = "ECalFwdHits"
savehcaltool = SimG4SaveCalHits("saveHCalHits",readoutNames = [hcalBarrelReadoutName])
savehcaltool.positionedCaloHits.Path = "HCalPositionedHits"
savehcaltool.caloHits.Path = "HCalHits"
saveexthcaltool = SimG4SaveCalHits("saveExtHCalHits",readoutNames = [extHcalReadoutName])
saveexthcaltool.positionedCaloHits.Path = "ExtHCalPositionedHits"
saveexthcaltool.caloHits.Path = "ExtHCalHits"
savehcalendcaptool = SimG4SaveCalHits("saveHCalEndcapHits", readoutNames = [hcalEndcapReadoutName])
savehcalendcaptool.positionedCaloHits.Path = "HCalEndcapPositionedHits"
savehcalendcaptool.caloHits.Path = "HCalEndcapHits"
savehcalfwdtool = SimG4SaveCalHits("saveHCalFwdHits", readoutNames = [hcalFwdReadoutName])
savehcalfwdtool.positionedCaloHits.Path = "HCalFwdPositionedHits"
savehcalfwdtool.caloHits.Path = "HCalFwdHits"

# next, create the G4 algorithm, giving the list of names of tools ("XX/YY")
from Configurables import SimG4SingleParticleGeneratorTool
pgun = SimG4SingleParticleGeneratorTool("SimG4SingleParticleGeneratorTool",saveEdm=True,
                particleName=particleType,energyMin=energy,energyMax=energy,etaMin=-.0,etaMax=.0,
                OutputLevel = DEBUG)

geantsim = SimG4Alg("SimG4Alg",
                       outputs= ["SimG4SaveCalHits/saveECalBarrelHits", "SimG4SaveCalHits/saveECalEndcapHits",
                                 "SimG4SaveCalHits/saveECalFwdHits", "SimG4SaveCalHits/saveHCalHits",
                                 "SimG4SaveCalHits/saveExtHCalHits", "SimG4SaveCalHits/saveHCalEndcapHits",
                                 "SimG4SaveCalHits/saveHCalFwdHits"],
                       eventProvider=pgun,
                       OutputLevel=INFO)

#Configure tools for calo reconstruction
# EM scale calibration
from Configurables import CalibrateCaloHitsTool
calibHcells = CalibrateCaloHitsTool("CalibrateHCal", invSamplingFraction="41.66")

from Configurables import CalibrateInLayersTool
calibEcalBarrel = CalibrateInLayersTool("CalibrateECalBarrel",
                                   # sampling fraction obtained using SamplingFractionInLayers from DetStudies package
                                   samplingFraction = [0.36571381189697705] * 1 + [0.09779064189677973] * 1 + [0.12564152224404024] * 1 + [0.14350599973146283] * 1 + [0.1557126972314961] * 1 + [0.16444759076233928] * 1 + [0.17097165096847836] * 1 + [0.17684775359805122] * 1 + [0.18181154293837265] * 1 + [0.18544247938196395] * 1 + [0.18922747431624687] * 1 + [0.21187001375505543] * 1,
                                   readoutName = ecalBarrelReadoutName,
                                   layerFieldName = "layer")

calibEcalEndcap = CalibrateCaloHitsTool("CalibrateECalEndcap", invSamplingFraction="13.89")
calibEcalFwd = CalibrateCaloHitsTool("CalibrateECalFwd", invSamplingFraction="303.03")
calibHcalEndcap = CalibrateCaloHitsTool("CalibrateHCalEndcap", invSamplingFraction="33.62")
calibHcalFwd = CalibrateCaloHitsTool("CalibrateHCalFwd", invSamplingFraction="1207.7")

# Create cells in ECal barrel
# 1. step - merge hits into cells with default Eta segmentation
# 2. step - rewrite the cellId using the Phi-Eta segmentation
from Configurables import CreateCaloCells
createEcalBarrelCellsStep1 = CreateCaloCells("CreateECalBarrelCellsStep1",
                               doCellCalibration=True,
                               calibTool = calibEcalBarrel,
                               addCellNoise=False, filterCellNoise=False,
                               OutputLevel=INFO,
                               hits="ECalBarrelHits",
                               cells="ECalBarrelCellsStep1")
# Ecal barrel cell positions
from Configurables import CreateVolumeCaloPositions
positionsEcalBarrel = CreateVolumeCaloPositions("positionsBarrelEcal", OutputLevel = INFO)
positionsEcalBarrel.hits.Path = "ECalBarrelCellsStep1"
positionsEcalBarrel.positionedHits.Path = "ECalBarrelPositions"
# Use Phi-Eta segmentation in ECal barrel
from Configurables import RedoSegmentation
resegmentEcalBarrel = RedoSegmentation("ReSegmentationEcal",
                             # old bitfield (readout)
                             oldReadoutName = ecalBarrelReadoutName,
                             # specify which fields are going to be altered (deleted/rewritten)
                             oldSegmentationIds = ["module"],
                             # new bitfield (readout), with new segmentation
                             newReadoutName = ecalBarrelReadoutNamePhiEta,
                             OutputLevel = INFO,
                             inhits = "ECalBarrelPositions",
                             outhits = "ECalBarrelCellsStep2")
createEcalBarrelCells = CreateCaloCells("CreateECalBarrelCells",
                               doCellCalibration=False,
                               addCellNoise=False, filterCellNoise=False,
                               OutputLevel=INFO,
                               hits="ECalBarrelCellsStep2",
                               cells="ECalBarrelCells")


# Create Ecal cells in endcaps
# 1. step - merge layer IDs
# 2. step - create cells
from Configurables import MergeLayers
mergelayersEcalEndcap = MergeLayers("MergeLayersEcalEndcap",
                   # take the bitfield description from the geometry service
                   readout = ecalEndcapReadoutName,
                   # cells in which field should be merged
                   identifier = identifierName,
                   volumeName = volumeName,
                   # how many cells to merge
                   merge = ecalEndcapNumberOfLayersToMerge,
                   OutputLevel = INFO)
mergelayersEcalEndcap.inhits.Path = "ECalEndcapHits"
mergelayersEcalEndcap.outhits.Path = "mergedECalEndcapHits"

createEcalEndcapCells = CreateCaloCells("CreateEcalEndcapCaloCells",
                                    doCellCalibration=True,
                                    calibTool=calibEcalEndcap,
                                    addCellNoise=False, filterCellNoise=False,
                                    OutputLevel=INFO)
createEcalEndcapCells.hits.Path="mergedECalEndcapHits"
createEcalEndcapCells.cells.Path="ECalEndcapCells"

# Create Ecal cells in forward
mergelayersEcalFwd = MergeLayers("MergeLayersEcalFwd",
                   # take the bitfield description from the geometry service
                   readout = ecalFwdReadoutName,
                   # cells in which field should be merged
                   identifier = identifierName,
                   volumeName = volumeName,
                   # how many cells to merge
                   merge = ecalFwdNumberOfLayersToMerge,
                   OutputLevel = INFO)
mergelayersEcalFwd.inhits.Path = "ECalFwdHits"
mergelayersEcalFwd.outhits.Path = "mergedECalFwdHits"

createEcalFwdCells = CreateCaloCells("CreateEcalFwdCaloCells",
                                 doCellCalibration=True,
                                 calibTool=calibEcalFwd,
                                 addCellNoise=False, filterCellNoise=False,
                                 OutputLevel=INFO)
createEcalFwdCells.hits.Path="mergedECalFwdHits"
createEcalFwdCells.cells.Path="ECalFwdCells"

# Create cells in HCal
# 1. step - merge hits into cells with the default readout
# 2. step - rewrite the cellId using the Phi-Eta segmentation
# 3. step - merge new cells corresponding to eta-phi segmentation
createHcalCells = CreateCaloCells("CreateHCaloCells",
                               doCellCalibration=True,
                               calibTool=calibHcells,
                               addCellNoise = False, filterCellNoise = False,
                               OutputLevel = INFO,
                               hits="HCalHits",
                               cells="HCalBarrelCellsStep1")

# Hcal barrel cell positions                                                                                                                                                                                                              
from Configurables import CreateVolumeCaloPositions
positionsHcalBarrel = CreateVolumeCaloPositions("positionsBarrelHcal", OutputLevel = INFO)
positionsHcalBarrel.hits.Path = "HCalBarrelCellsStep1"
positionsHcalBarrel.positionedHits.Path = "HCalBarrelPositions"
# Use Phi-Eta segmentation in Hcal barrel                                                                                                                                                                                                  
from Configurables import RedoSegmentation
resegmentHcalBarrel = RedoSegmentation("ReSegmentationHcal",
                             # old bitfield (readout)                                                                                                                                                                                      
                             oldReadoutName = hcalBarrelReadoutName,
                             # specify which fields are going to be altered (deleted/rewritten)                                                                                                                                             
                             oldSegmentationIds = ["module","row","layer","eta","phi"],
                             # new bitfield (readout), with new segmentation                                                                                                                                                              
                             newReadoutName = hcalBarrelReadoutNamePhiEta,
                             OutputLevel = INFO,
                             inhits = "HCalBarrelPositions",
                             outhits = "HCalBarrelCellsStep2")
createHcalBarrelCells = CreateCaloCells("CreateHCalBarrelCells",
                               doCellCalibration=False,
                               addCellNoise=False, filterCellNoise=False,
                               OutputLevel=INFO,
                               hits="HCalBarrelCellsStep2",
                               cells="HCalBarrelCells")

# Hcal extended barrel cells
createExtHcalCells = CreateCaloCells("CreateExtHcalCaloCells",
                                  doCellCalibration=True,
                                  calibTool=calibHcells,
                                  addCellNoise = False, filterCellNoise = False,
                                  OutputLevel = INFO,
                                  hits="ExtHCalHits",
                                  cells="HCalExtBarrelCells")

# Create Hcal cells in endcaps
mergelayersHcalEndcap = MergeLayers("MergeLayersHcalEndcap",
                   # take the bitfield description from the geometry service
                   readout = hcalEndcapReadoutName,
                   # cells in which field should be merged
                   identifier = identifierName,
                   volumeName = volumeName,
                   # how many cells to merge
                   merge = hcalEndcapNumberOfLayersToMerge,
                   OutputLevel = INFO)
mergelayersHcalEndcap.inhits.Path = "HCalEndcapHits"
mergelayersHcalEndcap.outhits.Path = "mergedHCalEndcapHits"

createHcalEndcapCells = CreateCaloCells("CreateHcalEndcapCaloCells",
                                    doCellCalibration=True,
                                    calibTool=calibHcalEndcap,
                                    addCellNoise=False, filterCellNoise=False,
                                    OutputLevel=INFO)
createHcalEndcapCells.hits.Path="mergedHCalEndcapHits"
createHcalEndcapCells.cells.Path="HCalEndcapCells"

# Create Hcal cells in forward
mergelayersHcalFwd = MergeLayers("MergeLayersHcalFwd",
                   # take the bitfield description from the geometry service
                   readout = hcalFwdReadoutName,
                   # cells in which field should be merged
                   identifier = identifierName,
                   volumeName = volumeName,
                   # how many cells to merge
                   merge = hcalFwdNumberOfLayersToMerge,
                   OutputLevel = INFO)
mergelayersHcalFwd.inhits.Path = "HCalFwdHits"
mergelayersHcalFwd.outhits.Path = "mergedHCalFwdHits"

createHcalFwdCells = CreateCaloCells("CreateHcalFwdCaloCells",
                                 doCellCalibration=True,
                                 calibTool=calibHcalFwd,
                                 addCellNoise=False, filterCellNoise=False,
                                 OutputLevel=INFO)
createHcalFwdCells.hits.Path="mergedHCalFwdHits"
createHcalFwdCells.cells.Path="HCalFwdCells"

out = PodioOutput("out",
                  OutputLevel=INFO)
out.outputCommands = ["drop *", "keep ECalBarrelCells", "keep ECalEndcapCells", "keep ECalFwdCells", "keep HCalBarrelCells", "keep HCalBarrelPositions", "keep HCalExtBarrelCells", "keep HCalEndcapCells", "keep HCalFwdCells", "keep GenParticles","keep GenVertices"]
out.filename = "output_fullCalo_SimAndDigi_pi50GeV_"+str(num_events)+"events.root"

#CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
geantsim.AuditExecute = True
createEcalBarrelCellsStep1.AuditExecute = True
positionsEcalBarrel.AuditExecute = True
resegmentEcalBarrel.AuditExecute = True
createEcalBarrelCells.AuditExecute = True
createEcalEndcapCells.AuditExecute = True
createEcalFwdCells.AuditExecute = True
createHcalCells.AuditExecute = True
createExtHcalCells.AuditExecute = True
createHcalEndcapCells.AuditExecute = True
createHcalFwdCells.AuditExecute = True
out.AuditExecute = True

ApplicationMgr(
    TopAlg = [geantsim,
              createEcalBarrelCellsStep1,
              positionsEcalBarrel,
              resegmentEcalBarrel,
              createEcalBarrelCells,
              mergelayersEcalEndcap,
              createEcalEndcapCells,
              mergelayersEcalFwd,
              createEcalFwdCells,
              createHcalCells,
              positionsHcalBarrel,
              resegmentHcalBarrel,
              createHcalBarrelCells,
              createExtHcalCells,
              mergelayersHcalEndcap,
              createHcalEndcapCells,
              mergelayersHcalFwd,
              createHcalFwdCells,
              out
              ],
    EvtSel = 'NONE',
    EvtMax   = int(num_events),
    ExtSvc = [geoservice, podioevent, geantservice, audsvc],
 )
