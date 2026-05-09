#
# File: RecFCCeeCalorimeter/tests/options/HCalPhiThetaCaloTool_test.py
# Author: scott snyder <snyder@bnl.gov>
# Date: Feb, 2026
# Purpose: Test for HCalPhiThetaCaloTool
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

hcalBarrelTool = C.HCalPhiThetaCaloTool(
    "hcalBarrelGeometryTool", readoutName="HCalBarrelReadout"
)
hcalEndcapTool = C.HCalPhiThetaCaloTool(
    "hcalEndcapGeometryTool", readoutName="HCalEndcapReadout"
)


# Sorry about the lack of whitespace and ugly formatting, but ruff-format
# complains if this is written legibly.
appmgr = C.ApplicationMgr(
    TopAlg=[
        C.k4__recCalo__HCalPhiCaloToolTestAlg(
            HCalBarrelTool=hcalBarrelTool,
            HCalEndcapTool=hcalEndcapTool,
            ExpectedBarrelReadout="HCalBarrelReadout",
            ExpectedEndcapReadout="HCalEndcapReadout",
            ExpectedBarrelCells=210944,
            ExpectedEndcapCells=80896,
        )
    ],
    ExtSvc=[geoSvc],
)
