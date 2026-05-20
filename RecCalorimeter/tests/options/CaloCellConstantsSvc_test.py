#
# File: RecCalorimeter/tests/options/CaloCellConstantsSvc_test.py
# Author: scott snyder <snyder@bnl.gov>
# Date: Jan, 2026
# Purpose: Test for CaloCellConstantsSvc
#

import Configurables as C

appmgr = C.ApplicationMgr()
appmgr.TopAlg += [C.k4__recCalo__CaloCellConstantsSvcTestAlg()]
