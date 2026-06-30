#!/bin/sh

SOURCE_DIR=`dirname $0`
workdir=$PWD

git clone --no-checkout https://github.com/key4hep/k4geo
cd k4geo
git checkout HEAD example
cd $workdir

ddsim --steeringFile=k4geo/example/SteeringFile_IDEA_o1_v03.py --compactFile=$K4GEO/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml -G --gun.distribution uniform --gun.particle e- --random.seed 1234567 --numberOfEvents=5

k4run $SOURCE_DIR/../options/runSiPMDualReadout.py

