#!/bin/bash

#
# Get ERA-5 data from Levante and prepare for use in hydrocode
# * Remap to 1x1 deg grid (see era5_1deg_grid file)
# * Make 6-hr means
#

# Usage: ./prepare_era5.sh 1979 U
# will prepare zonal wind for year 1979

#SBATCH -N 1 
#SBATCH -n 8
#SBATCH -t 2-00:00:00 
#SBATCH -p shared 
#SBATCH -A bb0519
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out

module load cdo

e5dir=/work/bk1099/data/ml00_1H/
e5dir_sf=/work/bk1099/data/sf00_1H/
tmpdir=/work/bb0519/from_Mistral/bb0519/b350090/hydrothermal/tmp/

max_jobs=2

year=$1
var=$2 

# Loop over all months
for (( month=1 ; month<=12 ; month++ )) ; do

    # Find last day in month
    mm=$( printf "%02d" ${month} )
    daymax=$( cal $(date +"${mm} ${year}") | awk 'NF {DAYS = $NF}; END {print DAYS}' )

    # Loop over all days
    for (( day=1 ; day<=${daymax} ; day++ )) ; do
        
        dd=$( printf "%02d" ${day} )

	# Loop over every 6 hours
        for (( hour=0 ; hour<=18 ; hour=hour+6 )) ; do
            
            hh1=$( printf "%02d" ${hour} )
            hh2=$( printf "%02d" $(( ${hour} + 5 )) )            
            startdate=${year}-${mm}-${dd}T${hh1}:00:00
            enddate=${year}-${mm}-${dd}T${hh2}:59:59
            #echo "$startdate $enddate"

	    # Use cdo to
	    # * Select a 6-hour interval (-select command)
	    # * Compute time mean over the interval (-timmean)
	    # * Remap to regular grid (-setgridtype or -sp2gp)
	    # * Remap to 1x1 grid (-remapbil)
	    # * Convert from GRIB to netCDF (-f nc -t ecmwf)
	    
	    # For u,v,t we need to use sp2gp to go from spectral to grid-point space, then remap
	    # For q, p, z we need to use setgridtype,regular to go from reduced Gaussian to regular, then remap
	    
	    if [[ "$var" == "U" ]] ; then
               uargs="-O -f nc -t ecmwf -remapbil,era5_1deg_grid -sp2gp,cubic -timmean -select,startdate=${startdate},enddate=${enddate}"
               ugrb=${e5dir}/${year}/E5ml00_1H_${year}-${mm}-${dd}_131
               unc=${tmpdir}/E5ml00_1H_${year}-${mm}-${dd}_${hh1}_ERA5_U.nc
               ulog=${tmpdir}/E5ml00_1H_${year}-${mm}-${dd}_${hh1}_ERA5_U.log
	    elif [[ "$var" == "V" ]] ; then
	       uargs="-O -f nc -t ecmwf -remapbil,era5_1deg_grid -sp2gp,cubic -timmean -select,startdate=${startdate},enddate=${enddate}"
               ugrb=${e5dir}/${year}/E5ml00_1H_${year}-${mm}-${dd}_132
               unc=${tmpdir}/E5ml00_1H_${year}-${mm}-${dd}_${hh1}_ERA5_V.nc
               ulog=${tmpdir}/E5ml00_1H_${year}-${mm}-${dd}_${hh1}_ERA5_V.log
	    elif [[ "$var" == "T" ]] ; then
	       uargs="-O -f nc -t ecmwf -remapbil,era5_1deg_grid -sp2gp,cubic -timmean -select,startdate=${startdate},enddate=${enddate}"
               ugrb=${e5dir}/${year}/E5ml00_1H_${year}-${mm}-${dd}_130
               unc=${tmpdir}/E5ml00_1H_${year}-${mm}-${dd}_${hh1}_ERA5_T.nc
               ulog=${tmpdir}/E5ml00_1H_${year}-${mm}-${dd}_${hh1}_ERA5_T.log
	    elif [[ "$var" == "Q" ]] ; then
	       uargs="-O -f nc -t ecmwf -remapbil,era5_1deg_grid -setgridtype,regular -timmean -select,startdate=${startdate},enddate=${enddate}"
               ugrb=${e5dir}/${year}/E5ml00_1H_${year}-${mm}-${dd}_133
               unc=${tmpdir}/E5ml00_1H_${year}-${mm}-${dd}_${hh1}_ERA5_Q.nc
               ulog=${tmpdir}/E5ml00_1H_${year}-${mm}-${dd}_${hh1}_ERA5_Q.log
	    elif [[ "$var" == "P" ]] ; then
	       uargs="-O -f nc -t ecmwf -remapbil,era5_1deg_grid -setgridtype,regular -timmean -select,startdate=${startdate},enddate=${enddate}"
               ugrb=${e5dir_sf}/${year}/E5sf00_1H_${year}-${mm}_134
	       unc=${tmpdir}/E5ml00_1H_${year}-${mm}-${dd}_${hh1}_ERA5_P.nc
               ulog=${tmpdir}/E5ml00_1H_${year}-${mm}-${dd}_${hh1}_ERA5_P.log
            elif [[ "$var" == "Z" ]] ; then
	       uargs="-O -f nc -t ecmwf -remapbil,era5_1deg_grid -sp2gp,cubic -timmean -select,startdate=${startdate},enddate=${enddate}"
               ugrb=${e5dir}/${year}/E5ml00_1H_${year}-${mm}-${dd}_129
	       unc=${tmpdir}/E5ml00_1H_${year}-${mm}-${dd}_${hh1}_ERA5_Z.nc
               ulog=${tmpdir}/E5ml00_1H_${year}-${mm}-${dd}_${hh1}_ERA5_Z.log
	    fi

	    # This will allow max_jobs to run simultaneously
            while (( $(jobs -p | wc -l) >=  max_jobs )); do sleep 5s; done
            (                
                echo " Work on ${ugrb} " 
                cdo ${uargs} ${ugrb} ${unc} > ${ulog}
                
                if [[ $? -eq 0 ]]; then
                    echo "Preparation OK"
                else
                    echo "ERROR"
                fi
            ) &

        done
    done
done

# Wait for all jobs to finish before finishing script
wait
