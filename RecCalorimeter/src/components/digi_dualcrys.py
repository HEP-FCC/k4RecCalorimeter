from Configurables import DualCrysCalDigi

digi = DualCrysCalDigi("DualCrystalDigitizer")
digi.calCollections = "SimCaloHits"
digi.outputCalCollection = "DigitizedCaloHits"
digi.outputRelCollection = "CaloHitLinks"
digi.CalThreshold = 0.03  # MeV
digi.calibrationCoeffcal = 120000.0
digi.maxCalHitEnergy = 2.0
digi.detectornameEcal = "DRCrystal"
digi.detectornameHcal = "DRFtubeFiber"

# Then include this in your ApplicationMgr.TopAlg or ComponentAccumulator if you're using GaudiConfig2
