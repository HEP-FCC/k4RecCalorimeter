#!/bin/bash

# Check that the Key4hep environment is set
if [[ -z "${KEY4HEP_STACK}" ]]; then
  echo "Error: Key4hep environment not set"
  exit 1
fi

# run the script with the cell positioning tools
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # workaround to have ctests working
k4run $SCRIPT_DIR/recoPositions_fullCaloSystem.py --IOSvc.Input=FCChh_sim_digi_reco_e.root --IOSvc.Output=FCChh_cell_positions_e.root || exit 1
k4run $SCRIPT_DIR/recoPositions_fullCaloSystem.py --IOSvc.Input=FCChh_sim_digi_reco_pi.root --IOSvc.Output=FCChh_cell_positions_pi.root || exit 1
