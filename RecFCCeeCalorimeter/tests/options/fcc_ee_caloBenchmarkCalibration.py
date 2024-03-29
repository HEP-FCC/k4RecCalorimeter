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
DETECTORPATH = os.environ.get("K4GEO", "")
DETECTORS = [
        'FCCee/ALLEGRO/compact/ALLEGRO_o1_v01/ALLEGRO_o1_v01.xml'
]
for det in DETECTORS:
    GEOSERVICE.detectors += [os.path.join(DETECTORPATH, det)]
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
