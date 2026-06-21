#!/bin/sh

SOURCE_DIR=`dirname $0`

ddsim --steeringFile=$K4GEO/example/SteeringFile_IDEA_o1_v03.py --compactFile=$K4GEO/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml -G --gun.distribution uniform --gun.particle e- --random.seed 1234567 --numberOfEvents=5

k4run $SOURCE_DIR/../options/runSiPMDualReadout.py

