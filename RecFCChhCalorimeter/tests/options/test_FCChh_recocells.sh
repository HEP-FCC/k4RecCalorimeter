#!/bin/bash


# Check that the Key4hep environment is set
if [[ -z "${KEY4HEP_STACK}" ]]; then
  echo "Error: Key4hep environment not set"
  exit 1
fi


compactFile=$K4GEO/FCChh/compact/FCChhBaseline/FCChh_DectMaster.xml

# particle gun
echo "Generating events and performing the Geant4 simulation with ddsim"
ddsim --enableGun --gun.distribution uniform --gun.energy "10*GeV" --gun.thetaMin "20*deg" --gun.thetaMax "160*deg" --gun.particle e- --numberOfEvents 10 --outputFile FCChh_sim_e.root --random.enableEventSeed --random.seed 4255 --compactFile=$compactFile  || exit 1
ddsim --enableGun --gun.distribution uniform --gun.energy "10*GeV" --gun.thetaMin "20*deg" --gun.thetaMax "160*deg" --gun.particle pi- --numberOfEvents 10 --outputFile FCChh_sim_pi.root --random.enableEventSeed --random.seed 4316 --compactFile=$compactFile || exit 1

# run the RECO step
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # workaround to have ctests working
k4run $SCRIPT_DIR/FCChh_recocells.py --IOSvc.Input=FCChh_sim_e.root --IOSvc.Output=FCChh_e_cells.root || exit 1
k4run $SCRIPT_DIR/FCChh_recocells.py --IOSvc.Input=FCChh_sim_pi.root --IOSvc.Output=FCChh_pi_cells.root || exit 1

