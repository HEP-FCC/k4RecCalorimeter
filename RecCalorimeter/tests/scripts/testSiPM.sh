#!/bin/sh

workdir=$PWD # ${PROJECT_SOURCE_DIR}

git clone --no-checkout https://github.com/key4hep/k4geo
cd k4geo

git checkout HEAD FCCee/IDEA/compact/IDEA_o1_v03
git checkout HEAD example

cd $workdir

ddsim --steeringFile=k4geo/example/SteeringFile_IDEA_o1_v03.py --compactFile=k4geo/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml -G --gun.distribution uniform --gun.particle e- --random.seed 1234567

k4run RecCalorimeter/tests/options/runSiPMDualReadout.py

