from Gaudi.Configuration import *
import os

# CI mode -> barrel wedge (IDEA_o2_v01_CI.xml); otherwise the full detector.
CI = bool(os.environ.get("IDEA_O2_CI"))

# ---- I/O ----
from k4FWCore import IOSvc, ApplicationMgr
from Configurables import EventDataSvc

io_svc = IOSvc("IOSvc")
io_svc.Input = "testIDEA_o2_v01_sim.root"
io_svc.Output = "testIDEA_o2_v01_reco.root"
io_svc.outputCommands = ["keep *"]

# ---- Geometry ----
from Configurables import GeoSvc

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [
    os.path.join(
        os.environ.get("K4GEO", ""),
        "FCCee/IDEA/compact/IDEA_o2_v01",
        "IDEA_o2_v01_CI.xml" if CI else "IDEA_o2_v01.xml",
    )
]
geoservice.OutputLevel = INFO

# ---- Tracks from gen particles ----
from Configurables import TracksFromGenParticles

tracksFromGenParticles = TracksFromGenParticles(
    "CreateTracksFromGenParticles",
    InputGenParticles=["MCParticles"],
    InputSimTrackerHits=[
        "VertexBarrelCollection",
        "VertexEndcapCollection",
        "DCHCollection",
        "SiWrBCollection",
        "SiWrDCollection",
    ],
    OutputTracks=["TracksFromGenParticles"],
    OutputMCRecoTrackParticleAssociation=["TracksFromGenParticlesAssociation"],
    ExtrapolateToECal=True,
    OnlyCaloReachingParticles=True,
    OutputLevel=INFO,
)

# ---- Optical calo digitization ----
from Configurables import CreateOpticalCaloCells


def optical_digi(name, optical, truth, out, link, calib, mask=False):
    return CreateOpticalCaloCells(
        name,
        OpticalHits=optical,
        TruthHits=truth,
        OutputCollection=out,
        OutputLinks=link,
        calibrationConstant=calib,
        maskCherenkovForTruthLink=mask,
        OutputLevel=INFO,
    )


# ECAL crystals: scint & Cherenkov share one deposit -> mask cherenkov in the link
ecalDigis = [
    optical_digi("createEcalScintCells", "SCEPCal_MainScounts", "SCEPCal_MainEdep",
                 "SCEPCal_digi_scint", "SCEPCal_scint_link", 1.0 / 1965.0, mask=True),
    optical_digi("createEcalCherenCells", "SCEPCal_MainCcounts", "SCEPCal_MainEdep",
                 "SCEPCal_digi_cheren", "SCEPCal_cheren_link", 1.0 / 97.75, mask=True),
]

# HCAL tubes: scint & Cherenkov are distinct cells -> full-cellID link (default)
hcalDigis = [
    optical_digi("createHcalBarrelScintCells", "DRBTScin", "DRTubeEdep",
                 "DRBTScin_digi", "DRBTScin_link", 1.0 / 206.25),
    optical_digi("createHcalBarrelCherenCells", "DRBTCher", "DRTubeEdep",
                 "DRBTCher_digi", "DRBTCher_link", 1.0 / 68.25),
]
if not CI:  # DR endcap exists only in the full geometry
    hcalDigis += [
        optical_digi("createHcalEndcapRScintCells", "DRETScinRight", "DRTubeEdep",
                     "DRETScinRight_digi", "DRETScinRight_link", 1.0 / 206.25),
        optical_digi("createHcalEndcapRCherenCells", "DRETCherRight", "DRTubeEdep",
                     "DRETCherRight_digi", "DRETCherRight_link", 1.0 / 68.25),
        optical_digi("createHcalEndcapLScintCells", "DRETScinLeft", "DRTubeEdep",
                     "DRETScinLeft_digi", "DRETScinLeft_link", 1.0 / 206.25),
        optical_digi("createHcalEndcapLCherenCells", "DRETCherLeft", "DRTubeEdep",
                     "DRETCherLeft_digi", "DRETCherLeft_link", 1.0 / 68.25),
    ]

# ---- Seed -> merge -> grow clustering (SCEPCal) ----
from Configurables import (
    CaloDrivenClusterSeeding,
    TrackDrivenClusterSeeding,
    ClusterSeedMerging,
    ClusterSeedGrower,
)

caloDrivenClusterSeed = CaloDrivenClusterSeeding(
    "CaloDrivenClusterSeed",
    InputCaloHitCollections=["SCEPCal_digi_scint"],
    OutputSeedsA="CaloDrivenSeedsA",
    OutputSeedsB="CaloDrivenSeedsB",
    ReadoutName="SCEPCal_Main",
    FieldStringsToFilter=["depth"],
    FieldValuesToFilter=[0],
    SeedEnergyThresholdA=0.04,
    SeedEnergyThresholdB=0.02,
    MinAboveThresholdNeighbours=2,
    VonNeumannDistance=2,
    OutputLevel=INFO,
)

trackDrivenClusterSeed = TrackDrivenClusterSeeding(
    "TrackDrivenClusterSeed",
    InputTrackCollection="TracksFromGenParticles",
    InputCaloHitCollections=["SCEPCal_digi_scint"],
    OutputSeedsC="TrackDrivenSeedsC",
    ReadoutName="SCEPCal_Main",
    FieldStringsToFilter=["depth"],
    FieldValuesToFilter=[0],
    SeedEnergyThreshold=0.02,
    TrackWindow=0.05,
    OutputLevel=INFO,
)

clusterSeedMerge = ClusterSeedMerging(
    "ClusterSeedMerge",
    CaloDrivenSeeds=["CaloDrivenSeedsA", "CaloDrivenSeedsB"],
    TrackDrivenSeeds=["TrackDrivenSeedsC"],
    OutputMergedSeeds="MergedSeeds",
    MergeDistance=45.0,
    W0=4.6,
    OutputLevel=INFO,
)

topoClusterGrower = ClusterSeedGrower(
    "TopoClusterGrower",
    InputSeeds=["MergedSeeds"],
    InputHits=["SCEPCal_digi_cheren", "SCEPCal_digi_scint"],
    OutputClusters="TopoGrownClusters",
    ReadoutName="SCEPCal_Main",
    FieldStringsToFilter=[],
    FieldValuesToFilter=[],
    FieldStringsToInclude=[],
    FieldValuesToInclude=[],
    GrowThreshold=0.0,
    HardThreshold=0.0,
    UnseededThreshold=0.02,
    W0=4.6,
    VNDistSeeded=2,
    VNDistUnseeded=2,
    MinUnseededHits=2,
    MaxOpeningAngle=0.3,
    AttachIsolatedHits=False,
    OutputLevel=INFO,
)

# ---- Truth links ----
from Configurables import CreateTruthLinks

createTruthLinksECAL = CreateTruthLinks(
    "CreateTruthLinksECAL",
    cell_hit_links=["SCEPCal_scint_link"],
    clusters=["TopoGrownClusters"],
    mcparticles="MCParticles",
    cell_mcparticle_links="SCEPCal_CaloHitMCParticleLinks",
    cluster_mcparticle_links="EcalClusterMCParticleLinks",
    OutputLevel=INFO,
)

hcal_links = ["DRBTScin_link", "DRBTCher_link"]
if not CI:
    hcal_links += [
        "DRETScinRight_link", "DRETCherRight_link",
        "DRETScinLeft_link", "DRETCherLeft_link",
    ]

createTruthLinksHCAL = CreateTruthLinks(
    "CreateTruthLinksHCAL",
    cell_hit_links=hcal_links,
    clusters=[],
    mcparticles="MCParticles",
    cell_mcparticle_links="DRTube_CaloHitMCParticleLinks",
    cluster_mcparticle_links="empty_ClusterMCParticleLinks",
    OutputLevel=INFO,
)

# ---- Services + app ----
from Configurables import UniqueIDGenSvc

ApplicationMgr(
    TopAlg=[tracksFromGenParticles]
    + ecalDigis
    + hcalDigis
    + [
        caloDrivenClusterSeed,
        trackDrivenClusterSeed,
        clusterSeedMerge,
        topoClusterGrower,
        createTruthLinksECAL,
        createTruthLinksHCAL,
    ],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc"), geoservice, UniqueIDGenSvc("uidSvc")],
    StopOnSignal=True,
)
