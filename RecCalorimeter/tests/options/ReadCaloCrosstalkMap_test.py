#
# File: RecCalorimeter/tests/options/ReadCaloCrosstalkMap_test.py
# Author: scott snyder <snyder@bnl.gov>
# Date: Jun, 2026
# Purpose: Test for ReadCaloCrosstalkMap
#

import os
import Configurables as C

fileName = "crosstalkTest.root"
ECAL_Barrel = 4

compactFile = "ALLEGRO_o1_v03.xml"
pathToDetector = (
    os.environ.get("K4GEO", "")
    + "/FCCee/ALLEGRO/compact/"
    + os.path.splitext(compactFile)[0]
)

geoSvc = C.GeoSvc("GeoSvc", detectors=[os.path.join(pathToDetector, compactFile)])

ecalBarrelTool = C.TubeLayerModuleThetaCaloTool(
    "ecalBarrelGeometryTool",
    readoutName="ECalBarrelModuleThetaMerged",
    activeVolumeName="LAr_sensitive",
    activeFieldName="layer",
    activeVolumesNumber=11,
    fieldNames=["system"],
    fieldValues=[ECAL_Barrel],
)

indexerSvc = C.k4__recCalo__CaloCellIndexerSvc(GeoTools=[ecalBarrelTool])

crosstalkTool = C.ReadCaloCrosstalkMap(
    "ReadCaloCrosstalkMap", fileName=fileName, detID=ECAL_Barrel
)

appmgr = C.ApplicationMgr(
    TopAlg=[
        C.k4__recCalo__ReadCaloCrosstalkMapTestAlg(
            Tool=crosstalkTool, fileName=fileName, detID=ECAL_Barrel
        )
    ],
    ExtSvc=[geoSvc, indexerSvc],
)
