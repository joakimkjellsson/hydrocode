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

.. code:: 
	
	make
	./psi.x




Any dataset with data on terrain-following coordinates can be used. 
Some models use sigma, :math:`p_k = \sigma_k * p_s`, and some use hybrid, :math:`p_k = a_k * p_0 + b_k * p_s`, or :math:`p_k = a_k + b_k * p_s`. 

The following data has been used

* Reanalysis
	- ERA-Interim
	- MERRA
* CMIP5 models (atmosphere component only)
	- CanESM2
	- CCSM4
	- CSIRO-Mk3-6-0
	- GFDL CM3
	- IPSL CM5A
	- NorESM1



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
- December 2013: Added support for CSIRO-Mk-3-6-0



