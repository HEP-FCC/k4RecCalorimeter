import os

from GaudiKernel.SystemOfUnits import GeV
from Gaudi.Configuration import INFO, WARNING

# Application manager
from Configurables import ApplicationMgr
ApplicationMgr().EvtSel = 'NONE'
ApplicationMgr().EvtMax = 2
ApplicationMgr().OutputLevel = INFO


# Data service
from Configurables import k4DataSvc
DATASERVICE = k4DataSvc("EventDataSvc")
ApplicationMgr().ExtSvc += [DATASERVICE]


# Detector geometry
from Configurables import GeoSvc
GEOSERVICE = GeoSvc("GeoSvc")
# if FCC_DETECTORS is empty, this should use relative path to working directory
DETECTORPATH = os.environ.get("FCCDETECTORS", "")
DETECTORS = [
    'Detector/DetFCCeeIDEA-LAr/compact/FCCee_DectMaster.xml',
]
GEOSERVICE.detectors = [os.path.join(DETECTORPATH, d) for d in DETECTORS]
GEOSERVICE.OutputLevel = INFO
ApplicationMgr().ExtSvc += [GEOSERVICE]


# Empty collections of calorimeter cells
from Configurables import CreateEmptyCaloCellsCollection
CREATEEMPTYECALCELLS = CreateEmptyCaloCellsCollection("CreateEmptyECalCells")
CREATEEMPTYECALCELLS.cells.Path = "EcalBarrelCells"
ApplicationMgr().TopAlg += [CREATEEMPTYECALCELLS]

CREATEEMPTYHCALCELLS = CreateEmptyCaloCellsCollection("CreateEmptyHCalCells")
CREATEEMPTYHCALCELLS.cells.Path = "HcalBarrelCells"
ApplicationMgr().TopAlg += [CREATEEMPTYHCALCELLS]


# Benchmark method calibration
from Configurables import CalibrateBenchmarkMethod
BENCHMARKCALIB = CalibrateBenchmarkMethod("CalibrateBenchmarkMethod",
                                          readoutNames=["ECalBarrelEta",
                                                        "HCalBarrelReadout"],
                                          energy=50*GeV,
                                          ECalSystemID=4,
                                          HCalSystemID=8,
                                          numLayersECal=12,
                                          firstLayerHCal=0,
                                          OutputLevel=WARNING)
BENCHMARKCALIB.ecalBarrelCells.Path = "EcalBarrelCells"
BENCHMARKCALIB.hcalBarrelCells.Path = "HcalBarrelCells"
ApplicationMgr().TopAlg += [BENCHMARKCALIB]
