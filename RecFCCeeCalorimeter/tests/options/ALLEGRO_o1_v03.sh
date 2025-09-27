#!/bin/bash

# set this to 1 to use ee->Z->qq Pythia events, or 0 to use particle guns for e and pi, barrel and endcap
# this only applies to the reconstruction. Simulation will be run on both ee->Z->qq events and on particle
# guns, since the output files are used by other reco scripts in this package
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
  exit 1
fi

# get the files needed for calibration, noise, neighbor finding, etc
if ! test -f ./cellNoise_map_endcapTurbine_electronicsNoiseLevel.root; then  # assumes that if the last file exists, all the other as well
  echo "Downloading files needed for reconstruction"
  download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/capacitances_ecalBarrelFCCee_theta.root"
  download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/cellNoise_map_electronicsNoiseLevel_ecalB_thetamodulemerged.root"
  download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/cellNoise_map_electronicsNoiseLevel_ecalB_thetamodulemerged_hcalB_thetaphi.root"
  download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/cellNoise_map_electronicsNoiseLevel_ecalB_ECalBarrelModuleThetaMerged_ecalE_ECalEndcapTurbine_hcalB_HCalBarrelReadout_hcalE_HCalEndcapReadout.root"
  download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/elecNoise_ecalBarrelFCCee_theta.root"
  download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/lgbm_calibration-CaloClusters.onnx"
  download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/lgbm_calibration-CaloTopoClusters.onnx"
  download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/neighbours_map_ecalB_thetamodulemerged.root"
  download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/neighbours_map_ecalE_turbine.root"
  download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/neighbours_map_ecalB_thetamodulemerged_ecalE_turbine_hcalB_hcalEndcap_phitheta.root"
  download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/cellNoise_map_endcapTurbine_electronicsNoiseLevel.root"

  # add here the lines to get the files for the photon ID
fi

# run the RECO step
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # workaround to have ctests working
if [ "$usePythia" -gt 0 ]; then
  k4run $SCRIPT_DIR/ALLEGRO_o1_v03_digi_reco.py
else
#  k4run $SCRIPT_DIR/ALLEGRO_o1_v03_digi_reco.py --IOSvc.Input=ALLEGRO_sim_e_barrel.root --IOSvc.Output=ALLEGRO_sim_digi_reco_e_barrel.root || exit 1
#  k4run $SCRIPT_DIR/ALLEGRO_o1_v03_digi_reco.py --IOSvc.Input=ALLEGRO_sim_e_endcap.root --IOSvc.Output=ALLEGRO_sim_digi_reco_e_endcap.root || exit 1
#  k4run $SCRIPT_DIR/ALLEGRO_o1_v03_digi_reco.py --IOSvc.Input=ALLEGRO_sim_pi_barrel.root --IOSvc.Output=ALLEGRO_sim_digi_reco_pi_barrel.root || exit 1
#  k4run $SCRIPT_DIR/ALLEGRO_o1_v03_digi_reco.py --IOSvc.Input=ALLEGRO_sim_pi_endcap.root --IOSvc.Output=ALLEGRO_sim_digi_reco_pi_endcap.root || exit 1
  k4run $SCRIPT_DIR/ALLEGRO_o1_v03_digi_reco.py --IOSvc.Input=ALLEGRO_sim_e.root --IOSvc.Output=ALLEGRO_sim_digi_reco_e.root || exit 1
#  k4run $SCRIPT_DIR/ALLEGRO_o1_v03_digi_reco.py --IOSvc.Input=ALLEGRO_sim_pi.root --IOSvc.Output=ALLEGRO_sim_digi_reco_pi.root || exit 1
fi
