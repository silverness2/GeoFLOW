
## Differences from Main Branch

### Download GitHub Repository
Your first step is to get the code from the GitHub reposity
```console
git clone https://github.com/NOAA-GSD/GeoFLOW.git
```

### Branch feature/hackathon
For the hackathon the code has been reduced to a small kernel executing 
the portion of code we plan to concentrate on. This code (and the respective 
inputs) are only contained within the **feature/hackathon** branch of the GitHub 
repository so make sure you are inside the correct branch before compiling.

```console
cd GeoFLOW
git checkout feature/hackathon
```

### Compiling
Compiling the code is identical to the instructions provided on the main 
[README](../../README.md). A pre-compiled Boost Library will/has be/been 
provided on the machines we will work on to eliminate that hurdle and it's 
recommended we all use the matching compiler those are built with.  

More info will be provided as we gain access to the Hackathon cluster.

### Provided Boost Library (on Hera)
```console
export BOOST_ROOT=/scratch2/BMC/gsd-hpcs/Bryan.Flynt/opt/boost/intel-2020.2
```

## Complete Example

### Clone the repository and create a **build** directory.

```console
User@machine: > git clone https://github.com/NOAA-GSD/GeoFLOW.git
User@machine: > cd GeoFLOW
User@machine: GeoFLOW> git checkout feature/hackathon
User@machine: GeoFLOW> mkdir build
User@machine: GeoFLOW> ls
analysis  build  cmake  CMakeLists.txt  DISCLAIMER  docs  extern  LICENCE.md  README.md  scripts  src  test
```

### Build the code using the proper environment setup.  

To make it easier, I typically use the following script which I run from within 
the top level GeoFLOW directory. For your first attempt, I'd recommend to copy 
and past each line to see messages if any errors result.
```bash
#!/bin/bash  

#                                                                                                                             
# This script should be placed inside the top level GeoFLOW directory.                                                        
# The following commands will attempt to ensure a build directory exists                                                      
# within that level directory                                                                                                 
#                                                                                                                             
if [ ! -d "build" ]; then
    echo "build/ directory is not found at current level"
    exit 1
fi

# Load our Modules                                                       
module purge
module load cmake
module load intel/2020.2
module load impi/2020.2

# Set environment variables for CMake to pick up
export CC=mpiicc
export CXX=mpiicpc
export FC=mpiifort
export BOOST_ROOT=/scratch2/BMC/gsd-hpcs/Bryan.Flynt/opt/boost/intel-2020.2

# Change into the "out of source" build directory
cd build

# Remove any cached files from previous attempts
rm -r cmake*
rm -r CMake*

# Configure using CMake                                                                                                       
cmake --log-level=VERBOSE ../.

# Make it                                                                                                                     
make install
```   

### Running the test stub

```console
User@machine: GeoFLOW> cd build/bin
User@machine: GeoFLOW/build/bin> ls
hackathon_driver  hack_input.jsn
User@machine: GeoFLOW/build/bin> module load intel/2020.2 impi/2020.2
User@machine: GeoFLOW/build/bin> ./hackathon_driver
main: ---------------------------derivative OK: 2.01001e-14 : direction=1 method: old
main: ---------------------------derivative OK: 2.01001e-14 : direction=1 method: new
```

