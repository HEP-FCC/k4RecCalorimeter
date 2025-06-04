from Configurables import GeoSvc
from Configurables import ApplicationMgr
from Configurables import CreateFCCeeCaloXTalkNeighbours
import os
from Gaudi.Configuration import INFO, DEBUG

# Detector geometry
geoservice = GeoSvc("GeoSvc")
# if K4GEO is empty, this should use relative path to working directory
path_to_detector = os.environ.get("K4GEO", "")
print(path_to_detector)
detectors_to_use = [
    'FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml'
]

# prefix all xmls with path_to_detector
geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]
geoservice.OutputLevel = INFO

# create the crosstalk neighbour file for ECAL barrel cells
neighbours = CreateFCCeeCaloXTalkNeighbours("xtalk_neighbours",
                                       outputFileName="xtalk_neighbours_map_ecalB_thetamodulemerged.root",
                                       readoutNames=["ECalBarrelModuleThetaMerged"],
                                       systemNames=["system"],
                                       systemValues=[4],
                                       activeFieldNames=["layer"],
                                       activeVolumesNumbers=[11],
                                       xtalkCoefRadial=0.7e-2,
                                       xtalkCoefTheta=0.2e-2,
                                       xtalkCoefDiagonal=0.04e-2,
                                       xtalkCoefTower=0.1e-2,
                                       OutputLevel=DEBUG)

# ApplicationMgr
ApplicationMgr(TopAlg=[],
               EvtSel='NONE',
               EvtMax=1,
               # order is important, as GeoSvc is needed by G4SimSvc
               ExtSvc=[geoservice, neighbours],
               OutputLevel=INFO
               )
