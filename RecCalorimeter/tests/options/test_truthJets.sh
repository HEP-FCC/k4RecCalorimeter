#!/bin/bash

# set this to 1 to use ee->Z->qq Pythia events, or 0 to use particle guns for e and pi, barrel and endcap
usePythia=0

# Define the unction for downloading files
# Attempt to download the file using wget, exit the script if wget failed
download_file() {
  local url="$1"
  wget "$url" || { echo "Download failed"; exit 1; }
}

# Check that the Key4hep environment is set
if [[ -z "${KEY4HEP_STACK}" ]]; then
  echo "Error: Key4hep environment not set"
  exit 1
fi

# run the SIM step (for debug do not run it if files already present. Comment the if and fi lines for production)
if [ "$usePythia" -gt 0 ]; then
  # download the events to reconstruct
  if ! test -f ./pythia_ee_z_qq_10evt.hepmc; then
    echo "Downloading files needed for simulation"
    download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/gen/pythia_ee_z_qq_10evt.hepmc"
  fi
  # if ! test -f truthJets_sim.root; then
  echo "Performing the Geant4 simulation with ddsim"
  ddsim --inputFiles pythia_ee_z_qq_10evt.hepmc --numberOfEvents -1 --outputFile truthJets_sim.root --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml || { retcode=$? ; echo "Simulation failed" ; exit $retcode ; }
  # fi
else
  #if ! test -f truthJets_sim.root; then
  echo "Generating events and performing the Geant4 simulation with ddsim"
  ddsim --enableGun --gun.distribution uniform --gun.energy "10*GeV" --gun.thetaMin "45*deg" --gun.thetaMax "135*deg" --gun.particle e- --numberOfEvents 2 --outputFile truthJets_sim.root --random.enableEventSeed --random.seed 4242 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml || exit 1
  #fi
fi

# run the RECO step
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # workaround to have ctests working
k4run $SCRIPT_DIR/test_truthJets.py --IOSvc.Input=truthJets_sim.root --IOSvc.Output=truthJets_rec.root || exit 1
