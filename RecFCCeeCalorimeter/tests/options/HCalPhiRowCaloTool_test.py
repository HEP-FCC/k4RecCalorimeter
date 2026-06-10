#
# File: RecFCCeeCalorimeter/tests/options/HCalPhiRowCaloTool_test.py
# Author: scott snyder <snyder@bnl.gov>
# Date: May, 2026
# Purpose: Test for HCalPhiTowCaloTool
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

hcalBarrelTool = C.HCalPhiRowCaloTool(
    "hcalBarrelGeometryTool", readoutName="HCalBarrelReadoutPhiRow"
)
hcalEndcapTool = C.HCalPhiRowCaloTool(
    "hcalEndcapGeometryTool", readoutName="HCalEndcapReadoutPhiRow"
)


# Sorry about the lack of whitespace and ugly formatting, but ruff-format
# complains if this is written legibly.
appmgr = C.ApplicationMgr(
    TopAlg=[
        C.k4__recCalo__HCalPhiCaloToolTestAlg(
            HCalBarrelTool=hcalBarrelTool,
            HCalEndcapTool=hcalEndcapTool,
            ExpectedBarrelReadout="HCalBarrelReadoutPhiRow",
            ExpectedEndcapReadout="HCalEndcapReadoutPhiRow",
            ExpectedBarrelCells=1031680,
            ExpectedEndcapCells=1164800,
        )
    ],
    ExtSvc=[geoSvc],
)
