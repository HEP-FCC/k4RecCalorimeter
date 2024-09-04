#!/bin/bash

# set-up the Key4hep environment if not already set
if [[ -z "${KEY4HEP_STACK}" ]]; then
  #  source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
  echo "Error: Key4hep environment not set"
  return 1
else
  echo "The Key4hep stack was already loaded in this environment."
fi

# download the events to reconstruct
if ! -test https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/gen/pythia_ee_z_qq_10evt.hepmc; then
  echo "Downloading files needed for simulation"
  wget https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/gen/pythia_ee_z_qq_10evt.hepmc
  retcode=$?
  if [ $retcode -ne 0 ]; then
    echo "Download failed"
    exit $retcode
  fi
fi
# run the SIM step (for debug do not run it if files already present. Comment the if and fi lines for production)
#if ! test -f ALLEGRO_sim_ee_z_qq.root; then
echo "Performing the Geant4 simulation with ddsim"
ddsim --inputFiles pythia_ee_z_qq_10evt.hepmc --numberOfEvents -1 --outputFile ALLEGRO_sim_ee_z_qq.root --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml
#fi
#if ! test -f ALLEGRO_sim_e_barrel.root; then
#echo "Generating events and performing the Geant4 simulation with ddsim"
#ddsim --enableGun --gun.distribution uniform --gun.energy "10*GeV" --gun.particle e- --numberOfEvents 10 --outputFile ALLEGRO_sim_e_barrel.root --random.enableEventSeed --random.seed 42 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml
#ddsim --enableGun --gun.distribution uniform --gun.energy "10*GeV" --gun.particle e- --numberOfEvents 10 --outputFile ALLEGRO_sim_e_endcap.root --random.enableEventSeed --random.seed 42 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml
#ddsim --enableGun --gun.distribution uniform --gun.energy "10*GeV" --gun.particle pi- --numberOfEvents 10 --outputFile ALLEGRO_sim_pi_barrel.root --random.enableEventSeed --random.seed 42 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml
#ddsim --enableGun --gun.distribution uniform --gun.energy "10*GeV" --gun.particle pi- --numberOfEvents 10 --outputFile ALLEGRO_sim_pi_endcap.root --random.enableEventSeed --random.seed 42 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml
#fi

retcode=$?
if [ $retcode -ne 0 ]; then
  echo "Simulation failed"
  exit $retcode
fi

# get the files needed for calibration, noise, neighbor finding, etc
if ! test -f ./neighbours_map_ecalB_thetamodulemerged_hcalB_thetaphi.root; then  # assumes that if the last file exists, all the other as well
  echo "Downloading files needed for reconstruction"
  wget https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/capacitances_ecalBarrelFCCee_theta.root
  wget https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/cellNoise_map_electronicsNoiseLevel_ecalB_thetamodulemerged.root
  wget https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/cellNoise_map_electronicsNoiseLevel_ecalB_thetamodulemerged_hcalB_thetaphi.root
  wget https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/elecNoise_ecalBarrelFCCee_theta.root
  wget https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/lgbm_calibration-CaloClusters.onnx
  wget https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/lgbm_calibration-CaloTopoClusters.onnx
  wget https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/neighbours_map_ecalB_thetamodulemerged.root
  wget https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/neighbours_map_ecalB_thetamodulemerged_hcalB_thetaphi.root

  # add here the lines to get the files for the photon ID
fi

# run the RECO step
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # workaround to have ctests working
k4run $SCRIPT_DIR/ALLEGRO_o1_v03_digi_reco.py
