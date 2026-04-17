#!/bin/bash

# set this to 1 to use ee->Z->qq Pythia events, or 0 to use particle guns for e and pi, barrel and endcap
# this only applies to the reconstruction. Simulation will be run on both ee->Z->qq events and on particle
# guns, since the output files are used by other reco scripts in this package
usePythia=0

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
if [ "$usePythia" -gt 0 ]; then
    # Z->qq
    # download the events to reconstruct
    if ! test -f ./pythia_ee_z_qq_10evt.hepmc; then
	echo "Downloading files needed for simulation"
	download_file "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/gen/pythia_ee_z_qq_10evt.hepmc"
    fi
    # run ddsim
    if ! test -f ALLEGRO_sim_ee_z_qq.root; then
	echo "Performing the Geant4 simulation with ddsim"
	ddsim --inputFiles pythia_ee_z_qq_10evt.hepmc --numberOfEvents 5 --outputFile FCChh_sim_ee_z_qq.root --compactFile $K4GEO/FCChh/compact/FCChhBaseline/FCChh_DectMaster.xml || { retcode=$? ; echo "Simulation failed" ; exit $retcode ; }
    fi
else
    # particle gun (for debug do not run it if files already present. Comment the if and fi lines for production)
    if ! test -f FCChh_sim_e.root; then
	echo "Generating events and performing the Geant4 simulation with ddsim"
	ddsim --enableGun --gun.distribution uniform --gun.momentumMin "50*GeV" --gun.momentumMax "50*GeV" --gun.thetaMin "0.5*deg" --gun.thetaMax "1.5*deg" --gun.particle e- --numberOfEvents 2 --outputFile FCChh_sim_e.root --random.enableEventSeed --random.seed 4255 --compactFile $K4GEO/FCChh/compact/FCChhBaseline/FCChh_DectMaster.xml || exit 1
	ddsim --enableGun --gun.distribution uniform --gun.momentumMin "20*GeV" --gun.momentumMax "20*GeV" --gun.thetaMin "80*deg" --gun.thetaMax "100*deg" --gun.particle pi+ --numberOfEvents 2 --outputFile FCChh_sim_pi.root --random.enableEventSeed --random.seed 4255 --compactFile $K4GEO/FCChh/compact/FCChhBaseline/FCChh_DectMaster.xml || exit 1
    fi
fi

# run the DIGI step
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # workaround to have ctests working
if [ "$usePythia" -gt 0 ]; then
    k4run $SCRIPT_DIR/runFullCaloSystem_Digitisation.py || exit 1
else
    if ! test -f FCChh_sim_digi_e.root; then
	k4run $SCRIPT_DIR/runFullCaloSystem_Digitisation.py --IOSvc.Input=FCChh_sim_e.root --IOSvc.Output=FCChh_sim_digi_e.root || exit 1
	k4run $SCRIPT_DIR/runFullCaloSystem_Digitisation.py --IOSvc.Input=FCChh_sim_pi.root --IOSvc.Output=FCChh_sim_digi_pi.root || exit 1
    fi
fi

# run the RECO step
if [ "$usePythia" -gt 0 ]; then
    k4run $SCRIPT_DIR/runFullCaloSystem_ReconstructionSW_noNoise.py || exit 1
else
    if ! test -f FCChh_sim_digi_reco_e.root; then
	k4run $SCRIPT_DIR/runFullCaloSystem_ReconstructionSW_noNoise.py --IOSvc.Input=FCChh_sim_digi_e.root --IOSvc.Output=FCChh_sim_digi_reco_e.root || exit 1
	k4run $SCRIPT_DIR/runFullCaloSystem_ReconstructionSW_noNoise.py --IOSvc.Input=FCChh_sim_digi_pi.root --IOSvc.Output=FCChh_sim_digi_reco_pi.root || exit 1
    fi
fi
