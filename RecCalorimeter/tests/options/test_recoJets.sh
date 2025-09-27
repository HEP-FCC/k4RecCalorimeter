#!/bin/bash

# Set this to 1 to use ee->Z->qq Pythia events, or 0 to use particle guns.
# This only applies to the reconstruction:
# simulation will be run on both ee->Z->qq events and on particle guns,
# since the output files are used by other reco scripts in this package
usePythia=1

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

# run the SIM step (for debug do not run it if files already present. Comment the if and fi lines for production)
# Z->qq
# download the events to reconstruct
if ! test -f ./pythia_ee_z_qq_10evt.hepmc; then
    echo "Downloading files needed for simulation"
    download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/gen/pythia_ee_z_qq_10evt.hepmc"
fi
# run ddsim
# if ! test -f ALLEGRO_sim_ee_z_qq.root; then
echo "Performing the Geant4 simulation with ddsim"
ddsim --inputFiles pythia_ee_z_qq_10evt.hepmc --numberOfEvents 5 --outputFile ALLEGRO_sim_ee_z_qq.root --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml || { retcode=$? ; echo "Simulation failed" ; exit $retcode ; }
# fi

# particle gun (for debug do not run it if files already present. Comment the if and fi lines for production)
# if ! test -f ALLEGRO_sim_e.root; then
echo "Generating events and performing the Geant4 simulation with ddsim"
ddsim --enableGun --gun.distribution uniform --gun.energy "10*GeV" --gun.thetaMin "20*deg" --gun.thetaMax "160*deg" --gun.particle e- --numberOfEvents 10 --outputFile ALLEGRO_sim_e.root --random.enableEventSeed --random.seed 4255 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml || exit 1
# ddsim --enableGun --gun.distribution uniform --gun.energy "10*GeV" --gun.thetaMin "90*deg" --gun.thetaMax "90*deg" --gun.particle e- --numberOfEvents 10 --outputFile ALLEGRO_sim_e_barrel.root --random.enableEventSeed --random.seed 4255 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml || exit 1
# ddsim --enableGun --gun.distribution uniform --gun.energy "10*GeV" --gun.thetaMin "160*deg" --gun.thetaMax "160*deg" --gun.particle e- --numberOfEvents 10 --outputFile ALLEGRO_sim_e_endcap.root --random.enableEventSeed --random.seed 3426 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml || exit 1
ddsim --enableGun --gun.distribution uniform --gun.energy "10*GeV" --gun.thetaMin "20*deg" --gun.thetaMax "160*deg" --gun.particle pi- --numberOfEvents 10 --outputFile ALLEGRO_sim_pi.root --random.enableEventSeed --random.seed 4316 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml || exit 1
# ddsim --enableGun --gun.distribution uniform --gun.energy "10*GeV" --gun.thetaMin "90*deg" --gun.thetaMax "90*deg" --gun.particle pi- --numberOfEvents 10 --outputFile ALLEGRO_sim_pi_barrel.root --random.enableEventSeed --random.seed 4316 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml || exit 1
# ddsim --enableGun --gun.distribution uniform --gun.energy "10*GeV" --gun.thetaMin "160*deg" --gun.thetaMax "160*deg" --gun.particle pi- --numberOfEvents 10 --outputFile ALLEGRO_sim_pi_endcap.root --random.enableEventSeed --random.seed 63462 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml || exit 1
# fi

# run the RECO step
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # workaround to have ctests working
if [ "$usePythia" -gt 0 ]; then
    k4run $SCRIPT_DIR/test_recoJets.py --IOSvc.Input=ALLEGRO_sim_ee_z_qq.root --IOSvc.Output=recoJets_rec.root || exit 1
else
    k4run $SCRIPT_DIR/test_recoJets.py --IOSvc.Input=ALLEGRO_sim_e.root --IOSvc.Output=recoJets_rec.root || exit 1
fi
