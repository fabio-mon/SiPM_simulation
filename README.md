# SiPM_simulation
#Optical simulation of SiPM reading scintillating crystal, particularly focused on LYSO

#Installing instructions on lxplus
   mkdir <work_dir>
   cd <work_dir>
   git clone git@github.com:fabio-mon/SiPM_simulation.git
   cd SiPM_simulation
   mkdir <build_dir>
   cd <build_dir>
   source /afs/cern.ch/sw/lcg/contrib/gcc/4.9/x86_64-slc6/setup.sh
   source /afs/cern.ch/sw/lcg/external/geant4/10.1.p02/x86_64-slc6-gcc49-opt-MT/CMake-setup.sh
   export CXX=/afs/cern.ch/sw/lcg/contrib/gcc/4.9/x86_64-slc6/bin/g++
   export CC=/afs/cern.ch/sw/lcg/contrib/gcc/4.9/x86_64-slc6/bin/gcc
   cmake -DGeant4_DIR=/afs/cern.ch/sw/lcg/external/geant4/10.1.p02/x86_64-slc6-gcc49-opt-MT/lib64/Geant4-10.1.2 ..

#to compile and get the executable (executable name = "exampleB4c") in the <build_dir>
   make

#when opening a new shell BEFORE running the executable 
   source /afs/cern.ch/sw/lcg/contrib/gcc/4.9/x86_64-slc6/setup.sh
   source /afs/cern.ch/sw/lcg/external/geant4/10.1.p02/x86_64-slc6-gcc49-opt-MT/CMake-setup.sh
   export CXX=/afs/cern.ch/sw/lcg/contrib/gcc/4.9/x86_64-slc6/bin/g++
   export CC=/afs/cern.ch/sw/lcg/contrib/gcc/4.9/x86_64-slc6/bin/gcc

#to run the executable  
   ./exampleB4c <configfile>

#NOTE: some configfile can be found in python folder

