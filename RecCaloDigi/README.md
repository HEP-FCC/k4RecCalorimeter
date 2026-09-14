# RecCaloDigi

Calorimeter digitisation and cell energy calibration for sampling calorimeters read out with
silicon or scintillator + PPD (SiPM/MPPC) sensors.

These algorithms are ported from the `CaloDigi/Realistic` family in
[MarlinReco](https://github.com/iLCSoft/MarlinReco/tree/master/CaloDigi/Realistic),
originally by D. Jeans (02/2016) with the charge-integration and timing rewrite by R. Ete
(11/2020).

## Algorithms

The package is organised as two technology-blind base classes with technology-specific
implementations inheriting from them. The suffix of each concrete algorithm names the energy
scale it produces.

| Algorithm | Purpose |
|---|---|
| `BaseDigitiser` | Base digitiser: timing cuts and charge integration, miscalibration, dead cells, electronics noise and dynamic range. `SimCalorimeterHit` → `CalorimeterHit`. |
| `DigitiserSiliconMip` | Silicon-specific digitisation; deposited GeV → MIP scale. |
| `DigitiserScinPpdNpe` | Scintillator + PPD digitisation including SiPM saturation and pixel statistics; deposited GeV → saturated photo-electrons (NPE). |
| `BaseCellEnergyCalibrator` | Base calibration: applies layer-group calibration constants to go from the digitised scale to shower GeV. |
| `CellEnergyCalibratorSiliconMip` | MIP → shower GeV. |
| `CellEnergyCalibratorScinPpdNpe` | PPD de-saturation, NPE → MIP → shower GeV. |

## The two-stage chain

Digitisation and calibration are separate algorithms, so the digitised hits remain available
on their own scale (MIP or NPE) in between. Both stages emit a `CaloHitSimCaloHitLink`
collection relating their output hits back to the originating simulated hits, and the
calibrator consumes the digitiser's *link* collection rather than its hit collection:

```
SimCalorimeterHits ─┐
                    ├─▶ DigitiserSiliconMip ─┬─▶ CalorimeterHits      (MIP scale)
EventHeader ────────┘                        └─▶ CaloHitLinks ─┐
                                                               │
                     CellEnergyCalibratorSiliconMip ◀──────────┘
                                    │
                                    ├─▶ CalorimeterHitsRec   (shower GeV)
                                    └─▶ CaloHitLinksRec
```

Both stages need `GeoSvc` for the cellID encoding, and the digitiser additionally needs
`UniqueIDGenSvc`, from which it draws a per-event random seed so that results depend only on
(run number, event number, algorithm name) and not on how events are scheduled across threads.

A minimal configuration of the silicon chain:

```python
from Gaudi.Configuration import INFO
from Configurables import GeoSvc, UniqueIDGenSvc
from Configurables import DigitiserSiliconMip, CellEnergyCalibratorSiliconMip
from k4FWCore import ApplicationMgr, IOSvc

geo_svc = GeoSvc("GeoSvc", detectors=["path/to/compact/detector.xml"])

io_svc = IOSvc()
io_svc.Input = "sim.root"
io_svc.Output = "digi.root"

digitiser = DigitiserSiliconMip(
    "EcalBarrelDigitiser",
    inputHitCollection="ECalBarrelCollection",
    inputHeaderCollection="EventHeader",
    outputHitCollection="EcalBarrelCalorimeterHits",
    outputRelationCollection="EcalBarrelCaloHitLinks",
    CaloType="em",
    CaloID="ecal",
    CaloLayout="barrel",
    calibration_mip=1.0e-4,  # average energy deposited by a MIP, in GeV
    threshold=0.5,
    thresholdUnit="MIP",
)

calibrator = CellEnergyCalibratorSiliconMip(
    "EcalBarrelCalibrator",
    inputLinkCollection="EcalBarrelCaloHitLinks",
    outputHitCollection="EcalBarrelCalorimeterHitsRec",
    outputRelationCollection="EcalBarrelCaloHitLinksRec",
    calibration_layergroups=[20, 20],            # group sizes, in layers
    calibration_factorsMipGev=[1.0e-4, 2.0e-4],  # one constant per group
)

ApplicationMgr(
    TopAlg=[digitiser, calibrator],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[geo_svc, UniqueIDGenSvc("UniqueIDGenSvc")],
    OutputLevel=INFO,
)
```

`calibration_layergroups` and `calibration_factorsMipGev` must have the same length and every
entry must be positive; the groups are expanded into a per-layer lookup at initialisation, so
together they should cover every layer the detector produces hits in. A hit in a layer beyond
the last group is calibrated to zero energy and warns.

Both algorithms decode the layer number out of the cellID, so `EncodingStringParameterName`
must name a DD4hep constant whose encoding contains a `layer` field. It defaults to
`GlobalCalorimeterReadoutID`; pointing it at a tracker encoding by mistake decodes `layer`
from the wrong bits rather than failing, which stays invisible for as long as only one
calibration group is in use.

## Why the calibrator takes a link collection as input

`BaseCellEnergyCalibrator` takes a `CaloHitSimCaloHitLinkCollection` as its only input, where
a `CalorimeterHitCollection` would be the natural expectation. This is deliberate.

The calibrator produces a *new* hit collection — the energies change, so the hits cannot be
modified in place — and it must also carry the association to the simulated hits forward, so
that the calibrated hits stay usable by anything needing MC truth. Taking the link collection
supplies both halves at once, and podio's typed accessors let the code say which half it wants
without anyone having to remember which way the link points:

```cpp
const auto hit = link.get<edm4hep::CalorimeterHit>();     // the digitised hit to calibrate
const auto sim = link.get<edm4hep::SimCalorimeterHit>();  // what the new link must point back at
```

Nothing has to be matched up by hand.

The alternative — taking the hit and link collections as two separate inputs — requires a
hit-to-link map to be built each event, and admits a configuration in which the two inputs do
not correspond to one another.

The cost of this choice is worth stating plainly: the digitised `CalorimeterHit` collection is
reached *through* the links and is therefore **not** a declared data dependency of the
calibrator. The scheduler sees it consume only the link collection.

## `CalorimeterHitType`

`CalorimeterHitType` (`CHT`) encodes the calorimeter type, ID, layout and layer number into
the `CalorimeterHit` type field. It originates in MarlinUtil and is also carried, in
near-identical form, by [k4GaudiPandora](https://github.com/key4hep/k4GaudiPandora). That
package compiles it privately into its Gaudi plugin module rather than exposing an installed
target, so there is nothing here to link against and the sources are duplicated instead.

This copy has been made a little safer to use: the three enums are scoped (`enum class`), and
`CaloType` gained the `unknown` fallback that `CaloID` and `Layout` already had, so a
calorimeter type matching nothing is no longer silently reported as electromagnetic. The
encoded integer is unchanged for every pre-existing value.

The duplication should be resolved by promoting these sources to a single shared, installed
location that both packages depend on.

## Documentation

The upstream package carries more detailed documentation of the digitisation model at
[MarlinReco `CaloDigi/Realistic/doc`](https://github.com/iLCSoft/MarlinReco/tree/master/CaloDigi/Realistic/doc).
Porting that material is deliberately left to a follow-up pull request.
