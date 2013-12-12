MODULE mod_supp

USE netcdf

CONTAINS  
   
   !!
   !! Subroutine for checking netCDF errors
   !!
   SUBROUTINE err(ierr2)
   
      INTEGER :: ierr2
   
      IF( ierr2 /= 0) THEN
         PRINT*, NF90_STRERROR(ierr2)
         STOP
      END IF
   
   END SUBROUTINE   
   
   !!
   !! Subroutine for sending a process to the system
   !!
   subroutine child(string)
      
      character(len=*) :: string   
   
      print*, string
      call system( trim(string) )
   
   end subroutine
   
   !!
   !! Subroutine for interpolating using CDO
   !!
   SUBROUTINE cdo_interp(infile,outfile,topofile,tmpfile,string,var,inc_lev,inc_ab,inc_p0)
      
      CHARACTER(LEN=*), INTENT(IN) :: infile, outfile, topofile, tmpfile, string, var
      CHARACTER(LEN=40) :: varstring
      
      INTEGER, INTENT(IN)   :: inc_ab,inc_p0,inc_lev
      
      if (inc_lev == 1) then
         varstring = TRIM(var)//',lon_bnds,lat_bnds,lev_bnds'
      else
         varstring = TRIM(var)//',lon_bnds,lat_bnds'
      end if
      
      PRINT*,' Isolate '//TRIM(var)
      PRINT*,'ncks -O -v '//TRIM(varstring)//' '//TRIM(infile)//' '//TRIM(tmpfile)
      CALL SYSTEM('ncks -O -v '//TRIM(varstring)//' '//TRIM(infile)//' '//TRIM(tmpfile))
      
      PRINT*,' Include grid cell area '
      PRINT*,'ncks -A -v areacella '//TRIM(topofile)//' '//TRIM(tmpfile)
      CALL SYSTEM('ncks -A -v areacella '//TRIM(topofile)//' '//TRIM(tmpfile))
      
      PRINT*,' Interpolate to '//TRIM(string)//' grid'
      PRINT*,'cdo -O remapbil,'//TRIM(string)//' '//TRIM(tmpfile)//' '//TRIM(outfile)
      CALL SYSTEM('cdo -O remapbil,'//TRIM(string)//' '//TRIM(tmpfile)//' '//TRIM(outfile)) 
      
      CALL SYSTEM('rm '//TRIM(tmpFile))
         
      IF (inc_ab == 1) THEN
         
         if (inc_p0 == 1) then
            varstring = 'a_bnds,b_bnds,p0'
         else if (inc_p0 == 0) then
            varstring = 'a_bnds,b_bnds'
         else if (inc_p0 == 2) then
            varstring = 'ap_bnds,b_bnds'
         end if
         PRINT*,'ncks -A -v '//TRIM(varstring)//' '//TRIM(infile)//' '//TRIM(outfile)
         CALL system('ncks -A -v '//TRIM(varstring)//' '//TRIM(infile)//' '//TRIM(outfile))
         
      END IF
   
   END SUBROUTINE

   
   
END MODULE