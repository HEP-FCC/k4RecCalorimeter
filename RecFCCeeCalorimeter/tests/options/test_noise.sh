#!/bin/bash

# Define the function for downloading files
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

# run ddsim
if ! test -f ./noise_sim.root; then
    echo "Generating events and performing the Geant4 simulation with ddsim"
    ddsim --enableGun --gun.distribution uniform --gun.momentumMin "10*GeV" --gun.momentumMax "10*GeV" --gun.thetaMin "0.5*deg" --gun.thetaMax "1.5*deg" --gun.particle geantino --numberOfEvents 1 --outputFile noise_sim.root --random.enableEventSeed --random.seed 4255 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml || exit 1
fi

# get the files needed for noise
if ! test -f ./elecNoise_ecalBarrelFCCee_theta.root; then  # assumes that if the last file exists, all the other as well
  echo "Downloading files needed to add noise"
  download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/elecNoise_ecalendcap.root"
  download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/elecNoise_ecalBarrelFCCee_theta.root"
fi

# run the DIGI step to add noise
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # workaround to have ctests working
# debug
nm -C $SCRIPT_DIR/../../RecFCCeeCalorimeter/libk4RecFCCeeCalorimeterPlugins.so | grep FromFile
# end debug
if ! test -f noise_digi.root; then
    k4run $SCRIPT_DIR/test_noise.py --IOSvc.Input=noise_sim.root --IOSvc.Output=noise_digi.root || exit 1
fi
