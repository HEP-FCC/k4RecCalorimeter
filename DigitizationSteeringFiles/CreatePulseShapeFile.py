from Configurables import GeoSvc
from Configurables import ApplicationMgr
from Configurables import CreateFCCeeCaloSignalShapes
import os
from Gaudi.Configuration import INFO, DEBUG

# Detector geometry
geoservice = GeoSvc("GeoSvc")
# if K4GEO is empty, this should use relative path to working directory
path_to_detector = os.environ.get("K4GEO", "")
detectors_to_use = [
    '%s/FCCee/ALLEGRO/compact/ALLEGRO_o1_v02/ALLEGRO_o1_v02.xml' % path_to_detector
]

# prefix all xmls with path_to_detector
geoservice.detectors = detectors_to_use
geoservice.OutputLevel = INFO

ecalBarrelReadoutName = "ECalBarrelModuleThetaMerged"  # name of the ECal barrel readout - found in GEANT4 xml

# create the signal shape file for ECAL barrel cells - based on xtalk neighbours code
neighbours = CreateFCCeeCaloSignalShapes("SignalShapes",
                                       outputFileName="SignalPulseShapes.root",
                                       readoutNames=[ecalBarrelReadoutName], # Can be a list
                                       systemNames=["system"],
                                       systemValues=[4],
                                       activeFieldNames=["layer"],
                                       activeVolumesNumbers=[11],
                                       pulseSampling = 10,
                                       pulseIndex = 2,
                                       pulseAmplitude = 1,
                                       OutputLevel=DEBUG)

# ApplicationMgr
ApplicationMgr(TopAlg=[],
               EvtSel='NONE',
               EvtMax=1,
               # order is important, as GeoSvc is needed by G4SimSvc
               ExtSvc=[geoservice, neighbours],
               OutputLevel=INFO
               )
