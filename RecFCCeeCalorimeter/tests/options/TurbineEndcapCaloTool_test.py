#
# File: RecFCCeeCalorimeter/tests/options/TurbineEndcapCaloTool_test.py
# Author: Giovanni Marchiori <giovanni.marchiori@cern.ch>
# Date: June, 2026
# Purpose: Test for TurbineEndcapCaloTool
#

import os
import Configurables as C

compactFile = "ALLEGRO_o1_v03.xml"
pathToDetector = (
    os.environ.get("K4GEO", "")
    + "/FCCee/ALLEGRO/compact/"
    + os.path.splitext(compactFile)[0]
)

geoSvc = C.GeoSvc("GeoSvc", detectors=[os.path.join(pathToDetector, compactFile)])

caloTool = C.TurbineEndcapCaloTool(
    "ecalEndcapGeometryTool",
    readoutName="ECalEndcapTurbine",
)


appmgr = C.ApplicationMgr(
    TopAlg=[C.k4__recCalo__TurbineEndcapCaloToolTestAlg(CalorimeterTool=caloTool)],
    ExtSvc=[geoSvc],
)
