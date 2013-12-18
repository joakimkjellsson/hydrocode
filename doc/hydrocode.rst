=======================
HydroCode documentation
=======================

--------------------------------------------------------------
A Fortran 95 code for calculating atmospheric stream functions
--------------------------------------------------------------

:Author: Joakim Kjellsson, Department of Meterology & Bolin Centre for Climate Research, Stockholm University
:Contact: joakimkjellsson@gmail.com

Compiling and running
=====================

A Makefile and runscript are provided with the code and should be modified to fit the user's needs. 

.. code:: 
	
	make
	./psi.x


The number of OpenMP threads is set by

.. code::
	
	export OMP_NUM_THREADS=8
	

An analysis run on several threads can consume a lot of memory, so it might be necessary to set 

.. code::
	
	ulimit -s unlimited
	export OMP_STACKSIZE=1g

or similar for other shells. 


Data
=====

Any dataset with data on terrain-following coordinates can be used. 
Some models use sigma, :math:`p_k = \sigma_k p_s`, and some use hybrid, :math:`p_k = a_k p_0 + b_k p_s`, or :math:`p_k = a_k + b_k p_s`. 

The following data has been used

* Reanalysis
	- ERA-Interim
	- MERRA
* CMIP5 models (atmosphere component only)
	- CanESM2
	- CCSM4
	- CNRM CM5
	- CSIRO Mk3-6-0
	- GFDL CM3
	- IPSL CM5A
	- MIROC5
	- NorESM1



Set up
======

Analysis can be controlled both via a few pre-processing flags at compile time and some values in a namelist at run time. 


Pre-processing
--------------

The pre-processing flags control the kinds of data that will be saved to file. 

========  =============================================================
Flag      Calculates & Outputs
========  =============================================================
-Dpsixy   Barotropic stream function in lon-lat coordinates
-Dpsixz   Zonal overturning stream function in model coordinates
-Dpsixr   Zonal overturining stream function in tracer coordinates
-Dpsiyz   Meridional overturning stream function in model coordinates
-Dpsiyr   Meridional overturning stream function in tracer coordinates
-Dpsirr   Stream function in tracer-tracer coordinates
========  =============================================================


Namelist
--------

The namelist controls a variety of settings.

==============  =======================================================
TIME                                                                              
==============  =======================================================
yearstart       Starting year                                                   
monstart        Starting month                                                   
daystart        Starting day
hourstart       Starting hour
hourstep        Time (in hours) between model output
intsend         Number of time steps to analyse
==============  =======================================================


==============  =======================================================
DIR                                                                              
==============  =======================================================
inDataDir       Folder for data to be analysed                                                   
outDataDir      Folder for storing analysis output                                                  
tmpDataDir      Folder for storing temporary files, e.g. scratch disk
topoDir         Folder for topography data, e.g. grid/orography
==============  =======================================================


==============  =======================================================
FILE
==============  =======================================================
project         Project to run. See list further down
prefix          Name of analysis run. Unique for each run
==============  =======================================================


==============  =======================================================
GRID
==============  =======================================================
imt             Number of zonal points of input
jmt             Number of meridional points of input
km              Number of vertical levels of input
mr              Number of tracer coordinate levels of output
mr2             Number of pressure levels for some of output
nrst            Number of tracer coordinates
lbas            Number of flags. Should be set to 1.
==============  =======================================================


==============  =======================================================
COORD
==============  =======================================================
rmin            Minimum pressure [hPa]
rmax            Maximum pressure [hPa]
tmin            Minimum temperature [K]
tmax            Maximum temperature [K]
smin            Minimum specific humidity [g/kg]
smax            Maximum specific humidity [g/kg]
dmin            Minimum dry static energy [kJ/kg]
dmax            Maximum dry static energy [kJ/kg]
mmin            Minimum moist static energy [kJ/kg]
mmax            Maximum moist static energy [kJ/kg]
amin            Minimum specific volume [m3/kg]
amax            Maximum specific volume [m3/kg]
==============  =======================================================


==============  =======================================================
MISC
==============  =======================================================
lverbose        Verbose mode. Show extra information
logp            Use log(pressure) instead of pressure
tweak_tmean     Use a temporal mean e.g. daily/monthly/annual mean
tweak_zmean     Use spatial average fields, e.g. re-grid or zonal mean.
tweak_tend      Output local tendency results
tweak_freq      If not 0 it controls the frequency of output
==============  =======================================================



Source code
===========

psi_main.f90
------------

The main program allocates much of the memory, sets up the OpenMP parallellisation, and contains the main algorithms. 




History
========

Parts of the original code was written by Kristofer Doos to calculate the thermohaline stream function for the NEMO ocean model. 
It was then adapted for the atmosphere (ERA-Interim) by Joakim Kjellsson. 

- January 2012: Original atmospheric code for ERA-Interim
- February 2013: Added OpenMP parallellisation and GRIB->netCDF conversion with CDO
- June 2013: Added options to use time-mean and/or zonal-mean variables. Also support for EC-Earth GFDL CM3. 
- July 2013: Restructured code into different modules, subroutines, etc. Unified different versions. 
- August 2013: Added preprocessing flags to control what stream functions are outputed. 
- September 2013: Added support from CanESM2, CCSM4, IPSL-CM5A, NorESM1
- December 2013: Added support for CSIRO-Mk-3-6-0, CNRM-CM5, MIROC5



