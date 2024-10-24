# k4RecCalorimeter - Key4hep Framework Components for Calorimeter Reconstruction

![key4hep-build](https://github.com/HEP-FCC/k4RecCalorimeter/actions/workflows/key4hep-build.yaml/badge.svg)

The components are available from the Key4hep stack on machines with CVMFS.

```
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
```

## Dependencies

* Gaudi
* ROOT
* EDM4HEP
* k4FWCore
* DD4hep
* k4geo
* ONNXRuntime


## Building

After fetching the repository, do
```
source /cvmfs/sw.hsf.org/key4hep/setup.sh
mkdir build install
cd build
cmake -DCMAKE_INSTALL_PREFIX=install ..
make -j4
make install
cd ..
k4_local_repo  # will update all necessary environment variables
```
