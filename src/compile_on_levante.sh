#!/bin/bash

# Load compiler and netCDF
module load gcc/11.2.0-gcc-11.2.0 netcdf-c/4.8.1-gcc-11.2.0  netcdf-fortran/4.5.3-gcc-11.2.0  cdo/2.0.5-gcc-11.2.0

# Set LD_LIBRARY_PATH
# (Why is this not done automatically?)
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/sw/spack-levante/netcdf-fortran-4.5.3-l2ulgp/lib/

# Make
make

# Copy
cp hydrocode.x ../hydrocode.x
