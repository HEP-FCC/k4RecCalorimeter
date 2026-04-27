#
# File: RecFCCeeCalorimeter/tests/options/TubeLayerModuleThetaCaloTool_test.py
# Author: scott snyder <snyder@bnl.gov>
# Date: Jan, 2026
# Purpose: Test for TubeLayerModuleThetaCaloTool
#

import os
import Configurables as C

compactFile = 'ALLEGRO_o1_v03.xml'
pathToDetector = os.environ.get('K4GEO','') + '/FCCee/ALLEGRO/compact/' + os.path.splitext(compactFile)[0]

geoSvc = C.GeoSvc('GeoSvc',
                  detectors = [os.path.join (pathToDetector, compactFile)])

caloTool = C.TubeLayerModuleThetaCaloTool \
    ('ecalBarrelGeometryTool',
     readoutName = 'ECalBarrelModuleThetaMerged',
     activeVolumeName = "LAr_sensitive",
     activeFieldName = "layer",
     activeVolumesNumber = 11,
     fieldNames = ["system"],
     fieldValues = [4])  # ECAL_Barrel


appmgr = C.ApplicationMgr \
    (TopAlg = [C.k4__recCalo__TubeLayerModuleThetaCaloToolTestAlg(CalorimeterTool = caloTool)],
     ExtSvc = [geoSvc])
