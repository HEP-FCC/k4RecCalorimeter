#!/bin/bash

# set this to 1 to use ee->Z->qq Pythia events, or 0 to use particle guns for e and pi, barrel and endcap
usePythia=1

# Define the unction for downloading files
# Attempt to download the file using wget, exit the script if wget failed
download_file() {
  local url="$1"
  wget "$url" || { echo "Download failed"; exit 1; }
}

# Check that the Key4hep environment is set
if [[ -z "${KEY4HEP_STACK}" ]]; then
  echo "Error: Key4hep environment not set"
  return 1
fi

# run the SIM step (for debug do not run it if files already present. Comment the if and fi lines for production)
if [ "$usePythia" -gt 0 ]; then
  # download the events to reconstruct
  if ! test -f ./pythia_ee_z_qq_10evt.hepmc; then
    echo "Downloading files needed for simulation"
    download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/gen/pythia_ee_z_qq_10evt.hepmc"
  fi
  # if ! test -f ALLEGRO_sim_ee_z_qq.root; then
  echo "Performing the Geant4 simulation with ddsim"
  ddsim --inputFiles pythia_ee_z_qq_10evt.hepmc --numberOfEvents -1 --outputFile ALLEGRO_sim_ee_z_qq.root --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml || { retcode=$? ; echo "Simulation failed" ; exit $retcode ; }
  # fi
else
  #if ! test -f ALLEGRO_sim_e_barrel.root; then
  echo "Generating events and performing the Geant4 simulation with ddsim"
  ddsim --enableGun --gun.distribution uniform --gun.energy "10*GeV" --gun.thetaMin "90*deg" --gun.thetaMax "90*deg" --gun.particle e- --numberOfEvents 100 --outputFile ALLEGRO_sim_e_barrel.root --random.enableEventSeed --random.seed 4255 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml || exit 1
  ddsim --enableGun --gun.distribution uniform --gun.energy "10*GeV" --gun.thetaMin "160*deg" --gun.thetaMax "160*deg" --gun.particle e- --numberOfEvents 100 --outputFile ALLEGRO_sim_e_endcap.root --random.enableEventSeed --random.seed 3426 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml || exit 1
  ddsim --enableGun --gun.distribution uniform --gun.energy "50*GeV" --gun.thetaMin "90*deg" --gun.thetaMax "90*deg" --gun.particle pi- --numberOfEvents 100 --outputFile ALLEGRO_sim_pi_barrel.root --random.enableEventSeed --random.seed 4316 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml || exit 1
  ddsim --enableGun --gun.distribution uniform --gun.energy "50*GeV" --gun.thetaMin "160*deg" --gun.thetaMax "160*deg" --gun.particle pi- --numberOfEvents 100 --outputFile ALLEGRO_sim_pi_endcap.root --random.enableEventSeed --random.seed 63462 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml || exit 1
  #fi
fi

# get the files needed for calibration, noise, neighbor finding, etc
if ! test -f ./neighbours_map_ecalB_thetamodulemerged_hcalB_thetaphi.root; then  # assumes that if the last file exists, all the other as well
  echo "Downloading files needed for reconstruction"
  download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/capacitances_ecalBarrelFCCee_theta.root"
  download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/cellNoise_map_electronicsNoiseLevel_ecalB_thetamodulemerged.root"
  download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/cellNoise_map_electronicsNoiseLevel_ecalB_thetamodulemerged_hcalB_thetaphi.root"
  download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/elecNoise_ecalBarrelFCCee_theta.root"
  download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/lgbm_calibration-CaloClusters.onnx"
  download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/lgbm_calibration-CaloTopoClusters.onnx"
  download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/neighbours_map_ecalB_thetamodulemerged.root"
  download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/neighbours_map_ecalB_thetamodulemerged_hcalB_hcalEndcap_phitheta.root"
  download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/neighbours_map_ecalE_turbine.root" 
  download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/cellNoise_map_endcapTurbine_electronicsNoiseLevel.root"
  
  # add here the lines to get the files for the photon ID
fi

# run the RECO step
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # workaround to have ctests working
if [ "$usePythia" -gt 0 ]; then
  k4run $SCRIPT_DIR/ALLEGRO_o1_v03_digi_reco.py
else
  k4run $SCRIPT_DIR/ALLEGRO_o1_v03_digi_reco.py --EventDataSvc.input=ALLEGRO_sim_e_barrel.root --out.filename=ALLEGRO_sim_digi_reco_e_barrel.root || exit 1
  k4run $SCRIPT_DIR/ALLEGRO_o1_v03_digi_reco.py --EventDataSvc.input=ALLEGRO_sim_e_endcap.root --out.filename=ALLEGRO_sim_digi_reco_e_endcap.root || exit 1
  k4run $SCRIPT_DIR/ALLEGRO_o1_v03_digi_reco.py --EventDataSvc.input=ALLEGRO_sim_pi_barrel.root --out.filename=ALLEGRO_sim_digi_reco_pi_barrel.root || exit 1
  k4run $SCRIPT_DIR/ALLEGRO_o1_v03_digi_reco.py --EventDataSvc.input=ALLEGRO_sim_pi_endcap.root --out.filename=ALLEGRO_sim_digi_reco_pi_endcap.root || exit 1
fi
