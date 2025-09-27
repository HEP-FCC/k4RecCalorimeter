#!/bin/bash

# set this to 1 to use ee->Z->qq Pythia events, or 0 to use particle guns
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

# run the RECO step
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # workaround to have ctests working
if [ "$usePythia" -gt 0 ]; then
    k4run $SCRIPT_DIR/test_truthJets.py --IOSvc.Input=ALLEGRO_sim_ee_z_qq.root --IOSvc.Output=truthJets_rec.root || exit 1
else
    k4run $SCRIPT_DIR/test_truthJets.py --IOSvc.Input=ALLEGRO_sim_e.root --IOSvc.Output=truthJets_rec.root || exit 1
fi
