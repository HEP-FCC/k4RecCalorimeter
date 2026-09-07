#!/bin/sh
# CI: IDEA_o2_v01 dual-readout calorimeter reco on a barrel wedge.
#   ddsim pions into the CI wedge -> optical digi + seed/grow clustering -> check clusters.
set -e

SOURCE_DIR=$(dirname "$0")
COMPACT="$K4GEO/FCCee/IDEA/compact/IDEA_o2_v01_CI/IDEA_o2_v01_CI.xml"

# Steering wires the DR/SCEPCal SDs.  It ships with k4geo alongside the CI geometry
# above, so both must come from the same k4geo the environment provides.
STEERING="$K4GEO/example/SteeringFile_IDEA_o2_v01.py"

for required in "$COMPACT" "$STEERING"; do
  if [ ! -f "$required" ]; then
    echo "$(basename "$0"): missing $required" >&2
    exit 1
  fi
done

# IDEA_O2_CI=1 selects the CI wedge in both the steering and the reco config.
export IDEA_O2_CI=1
ddsim --compactFile="$COMPACT" --steeringFile="$STEERING" \
      --runType batch -G --numberOfEvents 5 \
      --gun.particle pi- --gun.energy 10*GeV \
      --gun.direction "0.966,0.259,0.01" --random.seed 1988301045 \
      --outputFile testIDEA_o2_v01_sim.root

k4run "$SOURCE_DIR/IDEA_o2_v01_reco.py"

# Sanity: digitization, truth links and clustering.
python3 - <<'EOF'
import sys
import podio.reading

nClusters = nCells = nSignal = nUnlinked = nCherenLinks = 0
cherenMasked = True

for frame in podio.reading.get_reader("testIDEA_o2_v01_reco.root").get("events"):
    nClusters += len(frame.get("TopoGrownClusters"))

    # the digitizer emits exactly one cell per optical hit that carries signal
    nSignal += sum(1 for h in frame.get("SCEPCal_MainScounts") if h.getEnergy() > 0)
    cells = frame.get("SCEPCal_digi_scint")
    nCells += len(cells)

    # every digitized cell must carry a truth link to its energy deposit
    scint = frame.get("SCEPCal_scint_link")
    linked = {l.getFrom().getCellID() for l in scint}
    nUnlinked += sum(1 for c in cells if c.getCellID() not in linked)

    # maskCherenkovForTruthLink: Cherenkov cells resolve to the same deposits as the
    # scintillation cells.  Without the mask this collection would come out empty.
    cheren = frame.get("SCEPCal_cheren_link")
    nCherenLinks += len(cheren)
    cherenMasked &= {l.getTo().getCellID() for l in cheren} <= {l.getTo().getCellID() for l in scint}

print(f"clusters {nClusters}, cells {nCells}/{nSignal}, "
      f"unlinked cells {nUnlinked}, cherenkov links {nCherenLinks}")

checks = [
    (nClusters > 0, "no topo clusters were grown"),
    (nCells == nSignal, "one cell per optical hit with signal was not produced"),
    (nUnlinked == 0, f"{nUnlinked} digitized cells carry no truth link"),
    (nCherenLinks > 0 and cherenMasked, "cherenkov cells are not linked to the scintillation deposits"),
]
for ok, msg in checks:
    if not ok:
        print("FAIL:", msg)

sys.exit(0 if all(ok for ok, _ in checks) else 1)
EOF
