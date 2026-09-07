#!/bin/sh
# CI: IDEA_o2_v01 dual-readout calorimeter reco on a barrel wedge.
#   ddsim pions into the CI wedge -> optical digi + seed/grow clustering -> check clusters.
set -e

SOURCE_DIR=$(dirname "$0")
COMPACT="$K4GEO/FCCee/IDEA/compact/IDEA_o2_v01_CI/IDEA_o2_v01_CI.xml"

# The steering file wires the DR/SCEPCal sensitive detectors.  k4geo ships it in the
# repository but does not install it, so it has to be fetched -- and from the same
# revision as the installed geometry, or the two can disagree.  The installed version
# is recorded in k4geoConfig.cmake; nightly builds track main.
K4GEO_PREFIX=$(dirname "$(dirname "$K4GEO")")
K4GEO_VERSION=$(sed -n 's/^set(k4geo_VERSION \([0-9.]*\)).*/\1/p' \
                "$K4GEO_PREFIX/lib/cmake/k4geo/k4geoConfig.cmake")
case "$K4GEO" in
  */nightlies/*|*/HEAD/*) K4GEO_REF=main ;;
  *) K4GEO_REF=$(echo "$K4GEO_VERSION" | awk -F. '{printf "v%02d-%02d", $1, $2}') ;;
esac
echo "Taking the steering file from k4geo $K4GEO_REF (installed version ${K4GEO_VERSION:-unknown})"

if [ ! -d k4geo ]; then
  git clone --no-checkout --depth 1 --branch "$K4GEO_REF" https://github.com/key4hep/k4geo
fi
( cd k4geo && git checkout "$K4GEO_REF" -- example )
STEERING="k4geo/example/SteeringFile_IDEA_o2_v01.py"

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
python3 "$SOURCE_DIR/check_IDEAo2_reco.py" testIDEA_o2_v01_reco.root
