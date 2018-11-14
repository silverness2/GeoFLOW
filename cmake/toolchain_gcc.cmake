#
# Toolchain file for GCC Compiler already in the path
# 
# To Use:
# > cmake -DCMAKE_TOOLCHAIN_FILE=<filename>
#

set(CMAKE_SYSTEM_NAME Linux)

set(CMAKE_C_COMPILER       "gcc")
set(CMAKE_CXX_COMPILER     "g++")
set(CMAKE_Fortran_COMPILER "gfortran")

set(MPI_C_COMPILER         "mpicc")
set(MPI_CXX_COMPILER       "mpic++")
set(MPI_Fortran_COMPILER   "mpifort")
