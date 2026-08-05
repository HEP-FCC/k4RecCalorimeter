#!/bin/sh
# CI: IDEA_o2_v01 dual-readout calorimeter reco on a barrel wedge.
#   ddsim pions into the CI wedge -> optical digi + seed/grow clustering -> check clusters.
set -e

SOURCE_DIR=$(dirname "$0")
COMPACT="$K4GEO/FCCee/IDEA/compact/IDEA_o2_v01/IDEA_o2_v01_CI.xml"

# Steering wires the DR/SCEPCal SDs; use $K4GEO's copy if present, else fetch it.
STEERING="$K4GEO/example/SteeringFile_IDEA_o2_v01.py"
if [ ! -f "$STEERING" ]; then
  [ -d k4geo ] || { git clone --no-checkout https://github.com/key4hep/k4geo && ( cd k4geo && git checkout HEAD example ); }
  STEERING="k4geo/example/SteeringFile_IDEA_o2_v01.py"
fi

# IDEA_O2_CI=1 selects the CI wedge in both the steering and the reco config.
export IDEA_O2_CI=1
ddsim --compactFile="$COMPACT" --steeringFile="$STEERING" \
      --runType batch -G --numberOfEvents 5 \
      --gun.particle pi- --gun.energy 10*GeV \
      --gun.direction "0.966,0.259,0.01" --random.seed 1988301045 \
      --outputFile testIDEA_o2_v01_sim.root

k4run "$SOURCE_DIR/IDEA_o2_v01_reco.py"

# Sanity: at least one topo cluster must have been grown.
python3 - <<'EOF'
import sys, podio.reading
r = podio.reading.get_reader("testIDEA_o2_v01_reco.root")
n = sum(len(f.get("TopoGrownClusters")) for f in r.get("events"))
print("TopoGrownClusters produced:", n)
sys.exit(0 if n > 0 else 1)
EOF
