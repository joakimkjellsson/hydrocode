MODULE mod_data
!!------------------------------------------------------------------------------
!!
!!  This is the module for reading data
!!
!!  To include another dataset, create another subroutine in this module 
!!  and call it in the main program.
!!
!!  There are subroutines a2cgrid_t, a2cgrid_u, and a2cgrid_v that interpolates
!!  data from a A-grid to C-grid so that T-point i=1,j=1 is surrounded by 
!!  u-points i=0,i=1 and v-points j=0,j=1
!!
!!
!!  Joakim Kjellsson
!!  February 2013
!!
!! 
!!------------------------------------------------------------------------------

   USE netcdf
   USE omp_lib
   
   USE mod_grid
   USE mod_vars
   USE mod_supp
   
   
   IMPLICIT NONE
   
!!------------------------------------------------------------------------------
   
   !! Id:s for netCDF, loop indices, and temporary variables 
   !! are kept private for this module
   INTEGER, PRIVATE                     ::  id_nc, id_u, id_v, id_t, id_q,     &
   &                                        id_x, id_y, id_ml, id_hl, id_a,    &
   &                                        id_w, id_p0, id_area, id_lon2,     &
   &                                        id_lat2,                           &
   &                                        id_ncu, id_ncv, id_ncw, id_nct,    &
   &                                        id_ncq, id_ncz, id_ncp, id_nca,    &
   &                                        id_lon, id_lat, id_lev,            &
   &                                        id_b, id_z, id_p, ierr,            &
   &                                        IMT2, JMT2, KM2, KM3
   INTEGER, PRIVATE                     ::  ji, jj, jk, ip, im, jp, jm,        &
   &                                        istep, istep2, ik, il, ku, kb
   
   REAL*4, PRIVATE                      ::  da, db, dp, dlon, dlat,            &
   &                                        pp, pc, pm, tv, undef2, p0,&
   &                                        zz, dz, tt, ss
   
   
   !! Matrices to fill
   INTEGER, ALLOCATABLE, DIMENSION(:,:,:), PRIVATE   ::  steps, steps2
   REAL*4, ALLOCATABLE, DIMENSION(:), PRIVATE        ::  aa, bb
   
   REAL*4, ALLOCATABLE, DIMENSION(:,:), PRIVATE      ::  pxy, zxy, pt, pu, pv
   REAL*4, ALLOCATABLE, DIMENSION(:,:), PRIVATE      ::  zh, th, qh, uh, vh, ph
   
   REAL*4, ALLOCATABLE, DIMENSION(:,:,:), PRIVATE    ::  txyz, qxyz, uxyz, vxyz,&
   &                                                 wxyz, zxyz, pxyz
   
   
   CHARACTER, PRIVATE                                ::  dataprefix*100,       &
   &                                                 grbFile*150,              &
   &                                                 gzFile*100,               &
   &                                                 tmpFile*150,              &
   &                                                 tmpFile2*150,             &
   &                                                 scaFile*150,              &
   &                                                 vecFile*150,              &
   &                                                 scaFile2*150,             &
   &                                                 vecFile2*150,             &
   &                                                 tFile*150,                &
   &                                                 qFile*150,                &
   &                                                 zFile*150,                &
   &                                                 pFile*150,                &
   &                                                 uFile*150,                &
   &                                                 vFile*150,                &
   &                                                 wFile*150,                &
   &                                                 tFile2*150,               &
   &                                                 qFile2*150,               &
   &                                                 zFile2*150,               &
   &                                                 pFile2*150,               &
   &                                                 uFile2*150,               &
   &                                                 vFile2*150,               &
   &                                                 wFile2*150,               &
   &                                                 topoFile*150,             &
   &                                                 topoFile2*150,            &
   &                                                 string*1000
   
   LOGICAL, PRIVATE                              ::  lfirst

!!------------------------------------------------------------------------------

CONTAINS   
   
   subroutine get_data_era()
   !!---------------------------------------------------------------------------
   !!
   !! This subroutine reads the input data and fills the arrays
   !! uflux, vflux, wflux, tem, sal, rho, geo, vol
   !! all of dimension IMT x JMT x KM
   !!
   !! Also fills vlon, vlat, vlev of dimensions IMT, JMT, KM respectively
   !!
   !! This subroutine is run only by the master thread
   !! Consider parallelizing in the future...
   !!
   !!---------------------------------------------------------------------------
   
!$OMP MASTER
   
   dtstep = FLOAT(hourstep) * 3600. 
   
   !!
   !! Allocate dzt, dxdy and dx the first time step
   !! and reset them
   !!
   IF (iyear == yearstart .AND. &
   &    imon == monstart  .AND. &
   &    iday == daystart  .AND. &
   &   ihour == hourstart ) THEN
      
      lfirst = .true.
      
      dzt(:,:,:,:) = 0.
      dxdy(:,:) = 0.
      dx(:,:) = 0.
      dy(:,:) = 0.
      vlat2(:) = 0.
   
   ELSE
      
      lfirst = .false.
      
   END IF
   
   !!
   !! Temporary arrays
   !!
   
   ALLOCATE( txyz(IMT,0:JMT,KM), &
   &         qxyz(IMT,0:JMT,KM), &
   &         uxyz(IMT,0:JMT,KM), &
   &         vxyz(IMT,0:JMT,KM), &
   &         pxy(IMT,0:JMT),     &
   &         pu(imt,jmt),        &
   &         pv(imt,0:jmt),      &
   &         pt(imt,jmt),        &
   &         zxy(IMT,0:JMT),     &
   &         aa(0:KM), bb(0:KM) )
   
   
   !!
   !! Store previous time step and reset current
   !!
   dzt(:,:,:,1) = dzt(:,:,:,2)
   dzt(:,:,:,2) = 0.
      
   !!
   !! File names
   !!
   dataprefix='0000/uvtqzp_00000000.0000'
   
   write (dataprefix(1:4),   '(i4)') iyear
   write (dataprefix(13:16), '(i4)') iyear
   write (dataprefix(17:18), '(i2.2)') imon
   write (dataprefix(19:20), '(i2.2)') iday
   write (dataprefix(22:23), '(i2.2)') ihour
   
   grbFile = TRIM(inDataDir)//TRIM(dataprefix)//'.grb'
   tmpFile = TRIM(inDataDir)//TRIM(dataprefix(6:25))//TRIM(prefix)//'.nc' 
   
   if (lverbose) then
      print*,grbFile
   end if
   
   !!
   !! Use CDO to copy the GRIB file to a netCDF file
   !!
   call system('cdo -t ecmwf -f nc copy '//TRIM(grbFile)//' '//TRIM(tmpFile)//&
   &' > '//TRIM(dataprefix(6:25))//'.log')
   
   if (tweak_zmean > 0) then
      
      if (lverbose) then
         print*,' Performing spatial averaging with CDO '
      end if
      
      scaFile = TRIM(inDataDir)//TRIM(dataprefix(6:25))//TRIM(prefix)//'_sca.nc'
      vecFile = TRIM(inDataDir)//TRIM(dataprefix(6:25))//TRIM(prefix)//'_vec.nc'
      scaFile2 = TRIM(inDataDir)//TRIM(dataprefix(6:25))//TRIM(prefix)//'_sca2.nc'
      vecFile2 = TRIM(inDataDir)//TRIM(dataprefix(6:25))//TRIM(prefix)//'_vec2.nc'
      
      CALL system('cdo -O selvar,T,Q,Z,LNSP '//TRIM(tmpFile)//' '//TRIM(scaFile))
      CALL system('cdo -O remapbil,'//TRIM(rstring)//' '//TRIM(scaFile)//' '//TRIM(scaFile2)) 
      CALL system('rm '//TRIM(scaFile))
      CALL system('cdo -O selvar,U,V '//TRIM(tmpFile)//' '//TRIM(vecFile))
      CALL system('cdo -O remapbic,'//TRIM(rstring)//' '//TRIM(vecFile)//' '//TRIM(vecFile2))
      CALL system('rm '//TRIM(vecFile))
      CALL system('cdo -O merge '//TRIM(scaFile2)//' '//TRIM(vecFile2)//' '//TRIM(tmpFile))
      CALL system('rm '//TRIM(scaFile2))
      CALL system('rm '//TRIM(vecFile2))
      
   end if
   
   
   !!
   !! Read data
   !!   
   ierr = NF90_OPEN( tmpFile, NF90_SHARE, id_nc)
   CALL err(ierr)
   
   ierr = NF90_INQ_DIMID (id_nc, 'lon', id_x)
   CALL err(ierr)
   ierr = NF90_INQ_DIMID (id_nc, 'lat', id_y)
   CALL err(ierr)
   ierr = NF90_INQ_DIMID (id_nc, 'lev', id_ml)
   CALL err(ierr)
   ierr = NF90_INQ_DIMID (id_nc, 'nhyi', id_hl)
   CALL err(ierr)
   
   ierr = NF90_INQUIRE_DIMENSION (id_nc, id_x, LEN=IMT2)
   CALL err(ierr)
   ierr = NF90_INQUIRE_DIMENSION (id_nc, id_y, LEN=JMT2)
   CALL err(ierr)
   ierr = NF90_INQUIRE_DIMENSION (id_nc, id_ml, LEN=KM2)
   CALL err(ierr)
   ierr = NF90_INQUIRE_DIMENSION (id_nc, id_hl, LEN=KM3)
   CALL err(ierr)
   
   !! Test so that IMT, JMT, KM that has been set agrees with the data
   IF (IMT2 /= IMT) THEN
      PRINT*,' Error: IMT in the data is not as you have set it! '
      PRINT*,IMT,IMT2
      STOP
   END IF
   IF (JMT2 /= JMT+1) THEN
      PRINT*,' Error: JMT in the data is not as you have set it! '
      PRINT*,JMT,JMT2
      STOP
   END IF
   IF (KM2 /= KM) THEN
      PRINT*,' Error: KM in the data is not as you have set it! '
      PRINT*,KM,KM2
      STOP
   END IF
   IF (KM3 /= KM+1) THEN
      PRINT*,' Error: KM+1 in the data is not as you have set it! '
      PRINT*,KM,KM+1,KM3
      STOP
   END IF
   
   IF (lfirst) THEN
      ierr = NF90_INQ_VARID (id_nc, 'lon', id_lon)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_nc, 'lat', id_lat)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_nc, 'lev', id_lev)
      CALL err(ierr)
   END IF
   
   ierr = NF90_INQ_VARID (id_nc, 'U', id_u) 
   CALL err(ierr)
   ierr = NF90_INQ_VARID (id_nc, 'V', id_v) 
   CALL err(ierr)
   ierr = NF90_INQ_VARID (id_nc, 'T', id_t) 
   CALL err(ierr)
   ierr = NF90_INQ_VARID (id_nc, 'Q', id_q) 
   CALL err(ierr)
   ierr = NF90_INQ_VARID (id_nc, 'Z', id_z) 
   CALL err(ierr)
   ierr = NF90_INQ_VARID (id_nc, 'LNSP', id_p) 
   CALL err(ierr)
   ierr = NF90_INQ_VARID (id_nc, 'hyai', id_a)
   CALL err(ierr)
   ierr = NF90_INQ_VARID (id_nc, 'hybi', id_b)
   CALL err(ierr)
   
   IF (lfirst) THEN
      
      ierr = NF90_GET_VAR (id_nc, id_lon, vlon, start=[1], count=[IMT])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_nc, id_lat, vlat2(0:JMT), start=[1], count=[JMT+1])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_nc, id_lev, vlev, start=[1], count=[KM])
      CALL err(ierr)
      
      if (tweak_zmean <= 0) then
         vlat2(0:jmt) = vlat2(jmt:0:-1)
      end if
      
      up = .true.
      
      call calc_dxdy()
      
   END IF
   
   
   
   ierr = NF90_GET_VAR (id_nc, id_u, uxyz(:,0:jmt,:), start=[1,1,1,1], count=[IMT,JMT+1,KM,1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_nc, id_v, vxyz(:,0:jmt,:), start=[1,1,1,1], count=[IMT,JMT+1,KM,1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_nc, id_t, txyz(:,0:jmt,:), start=[1,1,1,1], count=[IMT,JMT+1,KM,1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_nc, id_q, qxyz(:,0:jmt,:), start=[1,1,1,1], count=[IMT,JMT+1,KM,1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_nc, id_z, zxy(:,0:jmt), start=[1,1,1,1], count=[IMT,JMT+1,1,1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_nc, id_p, pxy(:,0:jmt), start=[1,1,1,1], count=[IMT,JMT+1,1,1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_nc, id_a, aa(0:km), start=[1], count=[KM+1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_nc, id_b, bb(0:km), start=[1], count=[KM+1])
   CALL err(ierr)
   
   pxy(:,:) = exp(pxy(:,:))
   
   if (tweak_zmean <= 0) then
      
      uxyz(:,0:jmt,:) = uxyz(:,jmt:0:-1,:)
      vxyz(:,0:jmt,:) = vxyz(:,jmt:0:-1,:)
      txyz(:,0:jmt,:) = txyz(:,jmt:0:-1,:)
      qxyz(:,0:jmt,:) = qxyz(:,jmt:0:-1,:)
      zxy(:,0:jmt) = zxy(:,jmt:0:-1)
      pxy(:,0:jmt) = pxy(:,jmt:0:-1)
   
   end if
         
   ierr = NF90_CLOSE(id_nc)
   CALL err(ierr)
   
   CALL system('rm '//TRIM(tmpFile)) !Remove the temporary file
   CALL system('rm '//TRIM(dataprefix(6:25))//'.log')  
   
   !!
   !! Interpolate surface pressure and topography to C-grid
   !!
   call a2cgrid_u(pxy(:,:),pu(:,:))
   call a2cgrid_v(pxy(:,:),pv(:,:))
   call a2cgrid_t(pxy(:,:),pt(:,:))
         
   call a2cgrid_t(zxy(:,:),geo_i(:,1:jmt,km))
   
   !!
   !! Set top pressure explicitly
   !!
   rho_i(:,:,0) = 10.
   
   do jk=1,KM
   
      !!
      !! A-grid -> C-grid
      !!
      call a2cgrid_t(txyz(:,:,jk),tem(:,:,jk))
      call a2cgrid_t(qxyz(:,:,jk),sal(:,:,jk))
      
      !! Pressure
      rho_i(:,:,jk) = aa(jk) + bb(jk)*pt(:,:)
      
      !! Grid box depth
      da = aa(jk)-aa(jk-1)
      db = bb(jk)-bb(jk-1)
      dzt(:,:,jk,2) = da + db * pt(:,:)
      
      !! Mass
      vol(:,:,jk) = dzt(:,:,jk,2) * dxdy(:,:) / dg
      
      !!
      !! Mass fluxes
      !!
      call a2cgrid_u(uxyz(:,:,jk),uxyz(:,:,jk))
      uflux(:,:,jk) = (da + db * pu(:,:)) * uxyz(:,:,jk) * dy(:,:) / dg
       
      call a2cgrid_v(vxyz(:,:,jk),vxyz(:,:,jk))
      vflux(:,:,jk) = (da + db * pv(:,:)) * vxyz(:,:,jk) * dx(:,:) / dg
      
   end do
   
   !!
   !! Calculate geopotential
   !!
   call int_geo()
   
   !! Geopotential and pressure at level mid-points
   geo(:,:,1:km) = 0.5 * (geo_i(:,:,0:km-1) + geo_i(:,:,1:km))
   rho(:,:,1:km) = 0.5 * (rho_i(:,:,0:km-1) + rho_i(:,:,1:km))
   
   deallocate ( txyz, qxyz, uxyz, vxyz, pxy, pt, pu, pv, zxy, aa, bb )
   
   !!
   !! Calculate vertical mass flux
   !!
   call calc_wflux()
      
!$OMP END MASTER

   return
   
   end subroutine
   
   
   
   
   SUBROUTINE get_data_ifs(run)
   !!---------------------------------------------------------------------------
   !!
   !! Subroutine to get data from IFS, the atmospheric component of EC-Earth.
   !! Data is stored every 6h with one file every month. 
   !!
   !!---------------------------------------------------------------------------
   
   CHARACTER, INTENT(IN) :: run*9
   
!$OMP MASTER
   
   dtstep = FLOAT(hourstep) * 3600.  
   
   !!
   !! Allocate dzt, dxdy and dx the first time step
   !! and reset them
   !!
   IF (iyear == yearstart .AND. &
   &    imon == monstart  .AND. &
   &    iday == daystart  .AND. &
   &   ihour == hourstart ) THEN
      
      lfirst = .true.
      
      !!
      !! Private arrays that will be kept throughout the integrations
      !!
      ALLOCATE( aa(0:KM), bb(0:KM), zxy(IMT,0:JMT) ) 
      
      dzt(:,:,:,:) = 0.
      dxdy(:,:) = 0.
      dx(:,:) = 0.
      dy(:,:) = 0.
      vlat2(:) = 0.
      aa(:) = 0.
      bb(:) = 0.
      
   ELSE
      
      lfirst = .false.
      
   END IF
   
   IF (iday == 1 .AND. ihour == 0) THEN
      istep = 0
   END IF
   istep = istep + 1
   
   !!
   !! Temporary arrays
   !!
   ALLOCATE( txyz(IMT,0:JMT,KM), &
   &         qxyz(IMT,0:JMT,KM), &
   &         uxyz(IMT,0:JMT,KM), &
   &         vxyz(IMT,0:JMT,KM), &
   &         pxy(IMT,0:JMT),     &
   &         pt(imt,jmt),        &
   &         pu(imt,jmt),        &
   &         pv(imt,0:jmt) )
   
   geo_i(:,:,:) = 0.
   txyz(:,:,:) = 0.
   qxyz(:,:,:) = 0.
   uxyz(:,:,:) = 0.
   vxyz(:,:,:) = 0.
   pxy(:,:) = 0.
   pt(:,:) = 0.
   pu(:,:) = 0.
   pv(:,:) = 0.
   
   !!
   !! Store previous time step and reset current
   !!
   dzt(:,:,:,1) = dzt(:,:,:,2)
   dzt(:,:,:,2) = 0.
   
   IF (lverbose) THEN
      PRINT*,' Day, month, step ',iday, imon, istep
   END IF
      
   IF (istep == 1 .or. lfirst) THEN
      
      !!
      !! File names
      !!
      SELECT CASE (TRIM(run))
         CASE ('IFS-HISTR')
            dataprefix='0000/ICMSHSHC1+000000'
         CASE ('IFS-RCP45')
            dataprefix='0000/ICMSHSS41+000000'
         CASE ('IFS-RCP85')
            dataprefix='0000/ICMSHSS81+000000'
      END SELECT
   
      write (dataprefix(1:4),   '(i4.4)') iyear
      write (dataprefix(16:19), '(i4.4)') iyear
      write (dataprefix(20:21), '(i2.2)') imon
      
      grbFile = TRIM(inDataDir)//TRIM(dataprefix)
      tmpFile = TRIM(tmpDataDir)//TRIM(dataprefix(6:21))//TRIM(prefix)
      
      qFile = TRIM(tmpFile)//'_Q.nc' 
      uFile = TRIM(tmpFile)//'_U.nc' 
      vecFile = TRIM(tmpFile)//'_U_vectors.nc'
      scaFile = TRIM(tmpFile)//'_U_scalars.nc' 
      vecFile2 = TRIM(tmpFile)//'_U_vectors2.nc'
      scaFile2 = TRIM(tmpFile)//'_U_scalars2.nc' 
      
      zFile = '/nobackup/vagn2/x_joakj/data/ifs/topo/ec-earth_t159-topo.nc'
      
      IF (lverbose) THEN
         PRINT*,grbFile
      END IF
      
      !!
      !! Use CDO to copy the GRIB file to netCDF files
      !!
      
      !!
      !! Take all spectral variables
      !! 1. Copy all model level variables to grib file 
      !!
      string = 'grib_copy -w "typeOfLevel=hybrid" '//&
      &           TRIM(grbFile)//' '//TRIM(tmpFile)//'.grb'
      IF (lverbose) THEN
         PRINT*,' Copying spectral variables (div, vort, w, T, Z0, LNSP) '
         PRINT*,trim(string)
      END IF
      CALL system( TRIM(string) )
      
      !! 2. Convert div,vort to u and v and save the file as netCDF
      string = 'cdo -O -f nc -t ecmwf dv2uvl '//TRIM(tmpFile)//'.grb '//TRIM(tmpFile)//'.nc'
      IF (lverbose) THEN
         PRINT*,' Transforming div,vort -> u,v '
         PRINT*,trim(string)
      END IF
      CALL system( TRIM(string) ) 
      
      !! 3. Transform from spectral to grid point coordinates
      string = 'cdo -O -f nc -t ecmwf sp2gpl '//TRIM(tmpFile)//'.nc '//TRIM(uFile)
      IF (lverbose) THEN
         PRINT*,' Transforming from spectral to grid point '
         PRINT*,trim(string)
      END IF
      CALL system( TRIM(string) ) 
      
      !! 4. Remove temporary files
      CALL system('rm -f '//TRIM(tmpFile)//'.grb')
      CALL system('rm -f '//TRIM(tmpFile)//'.nc')
      
      string = 'cdo -O selvar,U,V '//TRIM(uFile)//' '//TRIM(vecFile)
      if (lverbose) then
         print*,' Separating u,v from other variables '
         print*,trim(string)
      end if
      call system(string)
      
      string = 'cdo -O remapbil,'//trim(rstring)//' '//TRIM(vecFile)//' '//TRIM(vecFile2)
      if (lverbose) then
         print*,' Re-mapping to '//trim(rstring)
         print*,trim(string)
      end if
      call system(string)
      
      string = 'cdo -O selvar,T,LNSP '//TRIM(uFile)//' '//TRIM(scaFile)
      if (lverbose) then
         print*,' Separating t and lnsp '
         print*,trim(string)
      end if
      call system(string)
      
      string = 'cdo -O remapbic,'//trim(rstring)//' '//TRIM(scaFile)//' '//TRIM(scaFile2)
      if (lverbose) then
         print*,' Re-mapping to '//trim(string)
         print*,trim(string)
      end if
      call system(string)
      
      string = 'cdo -O merge '//TRIM(scaFile2)//' '//TRIM(vecFile2)//' '//TRIM(uFile)
      if (lverbose) then
         print*,' Merging u,v with t,lnsp '
         print*,trim(string)
      end if
      call system(string)
      
      string = 'rm -f '//trim(scaFile)//' ; rm -f '//trim(vecFile)//&
      &        ' ; rm -f '//trim(scaFile2)//' ; rm -f '//trim(vecFile2)
      if (lverbose) then
         print*,' Cleaning up '
         print*,trim(string)
      end if
      call system(string)
      
      
      !! 5. Open the file and read the identifiers for the variables
      ierr = NF90_OPEN( uFile, NF90_SHARE, id_ncu)
      CALL err(ierr)
      
      ierr = NF90_INQ_VARID (id_ncu, 'U', id_u) 
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'V', id_v) 
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'T', id_t) 
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'LNSP', id_p) 
      CALL err(ierr)
      
      !!
      !! Q file
      !! 1. Copy q on model levels to grib file
      !!
      dataprefix(9:10) = 'GG'
      grbFile = TRIM(inDataDir)//TRIM(dataprefix)
      tmpFile = TRIM(tmpDataDir)//TRIM(dataprefix(6:21))//TRIM(prefix)
      string = 'grib_copy -w "typeOfLevel=hybrid,shortName=q" '//&
      &           TRIM(grbFile)//' '//TRIM(tmpFile)//'.grb'
      IF (lverbose) THEN
         PRINT*,' Copying specific humidity (q) '
         PRINT*,trim(string)
      END IF
      CALL system( TRIM(string) )
      
      !! 2. Copy the contents to a netCDF file and transform from reduced 
      !!    gaussian grid to regular grid
      string = 'cdo -O -R -f nc -t ecmwf -r copy '//TRIM(tmpFile)//'.grb '//TRIM(qFile)
      IF (lverbose) THEN
         PRINT*,' Transforming from reduced gaussian to regular grid point '
         PRINT*,trim(string)
      END IF
      CALL system( TRIM(string) ) 
      CALL system('rm -f '//TRIM(tmpFile)//'.grb')
      
      string = 'cdo -O remapbil,'//trim(rstring)//' '//TRIM(qFile)//' '//TRIM(scaFile)
      if (lverbose) then
         print*,' Re-mapping q to '//trim(rstring)
         print*,trim(string)
      end if
      call system(string)
      call system('rm -f '//trim(qFile))
      qFile = TRIM(scaFile)
      
      !! 3. Open the file
      ierr = NF90_OPEN( qFile, NF90_SHARE, id_ncq)
      CALL err(ierr)
      
      ierr = NF90_INQ_VARID (id_ncq, 'Q', id_q) 
      CALL err(ierr)
      
      
      !!
      !! Check dimensions
      !!
      ierr = NF90_INQ_DIMID (id_ncu, 'lon', id_x)
      CALL err(ierr)
      ierr = NF90_INQ_DIMID (id_ncu, 'lat', id_y)
      CALL err(ierr)
      ierr = NF90_INQ_DIMID (id_ncu, 'lev', id_ml)
      CALL err(ierr)
      ierr = NF90_INQ_DIMID (id_ncu, 'nhyi', id_hl)
      CALL err(ierr)
      
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_x, LEN=IMT2)
      CALL err(ierr)
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_y, LEN=JMT2)
      CALL err(ierr)
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_ml, LEN=KM2)
      CALL err(ierr)
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_hl, LEN=KM3)
      CALL err(ierr)
      
      !! Test so that IMT, JMT, KM that has been set agrees with the data
      IF (IMT2 /= IMT) THEN
         PRINT*,' Error: IMT in the data is not as you have set it! '
         PRINT*,IMT,IMT2
         STOP
      END IF
      IF (JMT2 /= JMT+1) THEN
         PRINT*,' Error: JMT in the data is not as you have set it! '
         PRINT*,JMT,JMT2
         STOP
      END IF
      IF (KM2 /= KM) THEN
         PRINT*,' Error: KM in the data is not as you have set it! '
         PRINT*,KM,KM2
         STOP
      END IF
      IF (KM3 /= KM+1) THEN
         PRINT*,' Error: KM+1 in the data is not as you have set it! '
         PRINT*,KM,KM+1,KM3
         STOP
      END IF
   
   END IF
   
   IF (lfirst) THEN
      
      IF (lverbose) THEN
         PRINT*,' Reading lon, lat, lev, A(k), B(k) '
      END IF
      
      ierr = NF90_INQ_VARID (id_ncu, 'lon', id_lon)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'lat', id_lat)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'lev', id_lev)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'hyai', id_a)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'hybi', id_b)
      CALL err(ierr)
      
      ierr = NF90_GET_VAR (id_ncu, id_lon, vlon, start=[1], count=[IMT])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_lat, vlat2(0:JMT), start=[1], count=[JMT+1])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_lev, vlev, start=[1], count=[KM])
      CALL err(ierr)
      
      ierr = NF90_GET_VAR (id_ncu, id_a, aa(0:KM), start=[1], count=[KM+1])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_b, bb(0:KM), start=[1], count=[KM+1])
      CALL err(ierr)
      
      !! 
      !! dx, dy, dxdy 
      !! dx at v-points, dy at u-points and dxdy at T-points
      !!
      call calc_dxdy()
      
      print*,vlat
      
      up = .true.
      
      !!
      !! Orography
      !!
      tmpFile = TRIM(tmpDataDir)//TRIM(dataprefix(6:21))//TRIM(prefix)//'_Z0.nc'
      string = 'cdo -O remapbic,'//trim(rstring)//' '//trim(zFile)//' '//trim(tmpFile)
      if (lverbose) then
         print*,' Re-mapping topography to '//trim(rstring)
         print*,trim(string)
      end if
      call system(trim(string))
      zFile = trim(tmpFile)
      
      !!
      !! Read 
      !!
      ierr = NF90_OPEN( zFile, NF90_SHARE, id_ncz)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncz, 'Z', id_z) 
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncz, id_z, zxy(:,0:JMT),  start=[  1,     1], &
      &                                                 count=[IMT, JMT+1])
      CALL err(ierr)
      
      ierr = NF90_CLOSE(id_ncz)
      CALL err(ierr)
      
      call system('rm -f '//trim(zFile))
      
      
   END IF
   
   if (lverbose) then
      print*,' reading u,v,t,p,q '
      print*,' imt, jmt, km ',imt,jmt,km
   end if
   ierr = NF90_GET_VAR (id_ncu, id_u, uxyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep],     &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncu, id_v, vxyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep], &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncu, id_t, txyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep], &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncu, id_p, pxy(:,0:JMT),  &
   &                                        start=[  1,     1,  1, istep], &
   &                                        count=[IMT, JMT+1,  1,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncq, id_q, qxyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep], &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   
         
   IF ( (iday == idmax(imon,iyear) .AND. ihour == 18) .OR. llast) THEN
      
      ierr = NF90_CLOSE(id_ncu)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncq)
      CALL err(ierr)
      
      !Remove the temporary files
      call system('rm '//TRIM(uFile)) 
      call system('rm '//TRIM(qFile)) 
      
   END IF
   
   
   !!
   !! Convert LNSP to ps
   !!
   pxy(:,:) = EXP( pxy(:,:) ) ![Pa]
   
   call a2cgrid_t(zxy(:,:),geo_i(:,:,km))
   call a2cgrid_t(pxy(:,:),pt(:,:))
   call a2cgrid_u(pxy(:,:),pu(:,:))
   call a2cgrid_v(pxy(:,:),pv(:,:))
   
   rho_i(:,:,0) = 5.
   
   DO jk=1,KM
   
      !!
      !! A-grid -> C-grid 
      !! u,v -> uflux,vflux
      !!
      call a2cgrid_t(txyz(:,:,jk),tem(:,:,jk))
      call a2cgrid_t(qxyz(:,:,jk),sal(:,:,jk))
      
      rho_i(:,:,jk)  =  aa(jk) + bb(jk) * pt(:,:)
      
      da = aa(jk)-aa(jk-1)
      db = bb(jk)-bb(jk-1)
      dzt  (:,:,jk,2)=  da + db * pt(:,:)
      vol  (:,:,jk)  =  dzt(:,:,jk,2) * dxdy(:,:) / dg 
      
      call a2cgrid_u(uxyz(:,:,jk),uxyz(:,:,jk))
      uflux(:,:,jk) = uxyz(:,:,jk) * (da + db * pu(:,:)) * dy(:,:) / dg 
            
      call a2cgrid_v(vxyz(:,:,jk),vxyz(:,:,jk))
      vflux(:,:,jk) =  vxyz(:,:,jk) * dx(:,:) * (da + db * pv(:,:)) / dg 
         
   END DO
      
   DEALLOCATE( txyz, qxyz, uxyz, vxyz, pxy, pt, pu, pv )
   
   
   !!
   !! Calculate geopotential
   !!
   call int_geo()
   
   
   !!
   !! Geopotential and pressue
   !!
   geo(:,:,1:km) = 0.5 * (geo_i(:,:,1:km) + geo_i(:,:,0:km-1))
   rho(:,:,1:km) = 0.5 * (rho_i(:,:,1:km) + rho_i(:,:,0:km-1))
   
   !!
   !! Calculate vertical mass flux
   !!
   call calc_wflux()
   
!$OMP END MASTER


   RETURN
   
   END SUBROUTINE

   
   
   
   
   SUBROUTINE get_data_merra()
   !!---------------------------------------------------------------------------
   !!
   !! Subroutine to get data from MERRA reanalysis from NASA.
   !! Data is stored every 3h with one file every day. 
   !!
   !!---------------------------------------------------------------------------
   
!$OMP MASTER
   
   dtstep = FLOAT(hourstep) * 3600. 
   
   !!
   !! Allocate dzt, dxdy and dx the first time step
   !! and reset them
   !!
   IF (iyear == yearstart .AND. &
   &    imon == monstart  .AND. &
   &    iday == daystart  .AND. &
   &   ihour == hourstart ) THEN
      
      lfirst = .true.
      
      !!
      !! Private arrays that will be kept throughout the integrations
      !!
      ALLOCATE( dxdy(IMT,JMT), dx(IMT,0:JMT), dy(IMT,JMT), &
      &         vlat2(0:JMT) ) 
      
      dxdy(:,:) = 0.
      dx(:,:) = 0.
      dy(:,:) = 0.
      vlat2(:) = 0.
      
   ELSE
      
      lfirst = .false.
      
   END IF
   
   IF (ihour == 0) THEN
      istep = 0
   END IF
   istep = istep + 1
   
   !!
   !! Temporary arrays
   !!
   ALLOCATE( txyz(IMT,JMT,KM), &
   &         qxyz(IMT,JMT,KM), &
   &         uxyz(IMT,JMT,KM), &
   &         vxyz(IMT,JMT,KM), &
   &         wxyz(IMT,JMT,KM), &
   &         zxyz(IMT,JMT,KM), &
   &         pxyz(IMT,JMT,0:KM), &
   &         pxy(IMT,JMT)      &
   &         )
   
   txyz(:,:,:) = 0.
   qxyz(:,:,:) = 0.
   uxyz(:,:,:) = 0.
   vxyz(:,:,:) = 0.
   wxyz(:,:,:) = 0.
   pxyz(:,:,:) = 0.
   zxyz(:,:,:) = 0.
   pxy(:,:) = 0.
   
   
   IF (lverbose) THEN
      PRINT*,' Day, month, step ',iday, imon, istep
   END IF
      
   IF (istep == 1) THEN
      
      !!
      !! File names
      !!
      dataprefix='0000/MERRA100.prod.assim.inst3_3d_asm_Cp.00000000.SUB.nc'
   
      WRITE (dataprefix(1:4),   '(i4)') iyear
      WRITE (dataprefix(42:45), '(i4)') iyear
      
      IF( imon < 10) THEN
         WRITE (dataprefix(47:47),'(i1)') imon
      ELSE
         WRITE (dataprefix(46:47),'(i2)') imon
      END IF   
      
      IF( iday < 10) THEN
         WRITE (dataprefix(49:49),'(i1)') iday
      ELSE
         WRITE (dataprefix(48:49),'(i2)') iday
      END IF   
      
      gzFile = TRIM(inDataDir)//TRIM(dataprefix)//'.gz'
      grbFile = TRIM(tmpDataDir)//TRIM(dataprefix(6:56))//'.gz'
      tmpFile = TRIM(tmpDataDir)//TRIM(dataprefix(6:56))
      
      IF (lverbose) THEN
         PRINT*,gzFile,grbFile,tmpFile
      END IF
      
      !!
      !! Copy compressed file to temporary location
      !!
      string = 'cp '//TRIM(gzFile)//' '//TRIM(grbFile)
      IF (lverbose) THEN
         PRINT*,' Copying .nc.gz file '
         PRINT*,string
      END IF
      CALL system( TRIM(string) )
      
      !!
      !! Unzip netCDF file
      !!
      string = 'gunzip --force '//TRIM(grbFile)
      IF (lverbose) THEN
         PRINT*,' Unzipping file '
         PRINT*,string
      END IF
      CALL system( TRIM(string) )
      
      !! Open the file and read the identifiers for the variables
      ierr = NF90_OPEN( tmpFile, NF90_SHARE, id_ncu)
      CALL err(ierr)
      
      ierr = NF90_INQ_VARID (id_ncu, 'u', id_u) 
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'v', id_v) 
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'omega', id_w) 
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 't', id_t) 
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'ps', id_p) 
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'qv', id_q) 
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'h', id_z) 
      CALL err(ierr)
      
      !!
      !! Check dimensions
      !!
      ierr = NF90_INQ_DIMID (id_ncu, 'longitude', id_x)
      CALL err(ierr)
      ierr = NF90_INQ_DIMID (id_ncu, 'latitude', id_y)
      CALL err(ierr)
      ierr = NF90_INQ_DIMID (id_ncu, 'levels', id_ml)
      CALL err(ierr)
      
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_x, LEN=IMT2)
      CALL err(ierr)
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_y, LEN=JMT2)
      CALL err(ierr)
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_ml, LEN=KM2)
      CALL err(ierr)
      
      !! Test so that IMT, JMT, KM that has been set agrees with the data
      IF (IMT2 /= IMT) THEN
         PRINT*,' Error: IMT in the data is not as you have set it! '
         PRINT*,IMT,IMT2
         STOP
      END IF
      IF (JMT2 /= JMT) THEN
         PRINT*,' Error: JMT in the data is not as you have set it! '
         PRINT*,JMT,JMT2
         STOP
      END IF
      IF (KM2 /= KM) THEN
         PRINT*,' Error: KM in the data is not as you have set it! '
         PRINT*,KM,KM2
         STOP
      END IF
      
   END IF
   
   IF (lfirst) THEN
      
      IF (lverbose) THEN
         PRINT*,' Reading longitudes, latitudes, levels '
      END IF
      
      ierr = NF90_INQ_VARID (id_ncu, 'longitude', id_lon)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'latitude', id_lat)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'levels', id_lev)
      CALL err(ierr)
      
      ierr = NF90_GET_VAR (id_ncu, id_lon, vlon, start=[1], count=[IMT])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_lat, vlat, start=[1], count=[JMT])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_lev, vlev(KM:1:-1), start=[1], count=[KM])
      CALL err(ierr)
      
      
      !! 
      !! dx, dy, dxdy 
      !! dx at v-points, dy at u-points and dxdy at T-points
      !!
      dlat = vlat(2) - vlat(1) ! Constant dy
      dlon = vlon(2) - vlon(1) ! Constant dx
      
      vlat2(1:JMT) = vlat(1:JMT) + dlat/2.  !Should be 88.75S - 90N
      vlat2(0) = vlat(1) - dlat/2.          !Should be 90S
      
      DO jj=0,JMT
      
         dx(:,jj) = ABS( deg * dlon * DCOS( vlat2(jj) * radian ) ) !dx at v-point
         
      END DO
      
      DO jj=1,JMT
         
         DO ji=1,IMT
            
            ip = ji+1
            IF (ip == IMT+1) THEN
               ip = 1
            END IF
            
            dy(ji,jj) = ABS( deg * dlat ) ! dy at u-point
            ! Take average of dx and multiply by dy
            ! Neglects curvature of latitude lines
            dxdy(ji,jj) = 0.5 * ( dx(ji,jj) + dx(ji,jj-1) ) * dy(ji,jj)
            
         END DO
         
      END DO
            
      IF (lverbose) THEN
         PRINT*,'dx(1,:)',dx(1,:)
         PRINT*,'dy(1,:)',dy(1,:)
         PRINT*,'dxdy(1,:)',dxdy(1,:)
      END IF
            
   END IF
   
   
   ierr = NF90_GET_VAR (id_ncu, id_u, uxyz(:,:,:), &
   &                                        start=[  1,   1,  1, istep],     &
   &                                        count=[IMT, JMT, KM,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncu, id_v, vxyz(:,:,:), &
   &                                        start=[  1,   1,  1, istep], &
   &                                        count=[IMT, JMT, KM,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncu, id_w, wxyz(:,:,:), &
   &                                        start=[  1,   1,  1, istep], &
   &                                        count=[IMT, JMT, KM,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncu, id_t, txyz(:,:,:), &
   &                                        start=[  1,   1,  1, istep], &
   &                                        count=[IMT, JMT, KM,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncu, id_p, pxy(:,:),  &
   &                                        start=[  1,   1,  1, istep], &
   &                                        count=[IMT, JMT,  1,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncu, id_q, qxyz(:,:,:), &
   &                                        start=[  1,   1,  1, istep], &
   &                                        count=[IMT, JMT, KM,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncu, id_z, zxyz(:,:,:), &
   &                                        start=[  1,   1,  1, istep], &
   &                                        count=[IMT, JMT, KM,     1])
   CALL err(ierr)      
   
   !!
   !! Get fill value
   !!
   ierr = NF90_GET_ATT (id_ncu, id_u, 'missing_value', undef2) 
   CALL err(ierr)
   
         
   IF ( ihour == (24-hourstep) .OR. llast) THEN
      
      ierr = NF90_CLOSE(id_ncu)
      CALL err(ierr)
      
      !Remove the temporary files
      CALL system('rm '//TRIM(tmpFile)) 
      
   END IF
   
   
   DO jk=1,KM
      
      !!
      !! Flip data vertically
      !!
      ik = KM - jk + 1
      tem(:,:,ik) = txyz(:,:,jk)
      sal(:,:,ik) = qxyz(:,:,jk)
      geo(:,:,ik) = dg * zxyz(:,:,jk)
      
      !!
      !! Data is on pressure levels, so 3D pressure is easily obtained 
      !!
      rho(:,:,ik) = vlev(ik) * 100. ![Pa]
      
   END DO
   
   
   DO jk=1,KM
      
      ik = KM - jk + 1
      
      !!
      !! u,v,omega -> uflux,vflux,wflux
      !!
      DO jj=1,JMT
         
         jm = jj-1
         jp = jj+1
         
         DO ji=1,IMT
            
            ip = ji+1
            IF (ip == IMT+1) THEN
               ip = 1
            END IF
            
            im = ji-1
            IF (im == 0) THEN
               im = IMT
            END IF
            
            IF (jk /= KM) THEN
               pxyz(ji,jj,jk) = 0.5 * (rho(ji,jj,jk) + rho(ji,jj,jk+1))
            ELSE
               pxyz(ji,jj,jk) = 120000. 
            END IF
            
            !!
            !! Calculate depth of grid box and whether it is active or not
            !!
            
            ! If surface pressure is higher than level interface
            IF ( pxy(ji,jj) >= pxyz(ji,jj,jk) ) THEN
               dp = pxyz(ji,jj,jk) - pxyz(ji,jj,jk-1)
            ! If surface pressure is lower than level interface, but higher than
            ! interface above (bottom box)
            ELSE IF ( pxy(ji,jj) <  pxyz(ji,jj,jk) .AND. &
            &         pxy(ji,jj) >= pxyz(ji,jj,jk-1) ) THEN
               dp = pxy(ji,jj) - pxyz(ji,jj,jk-1)
            ! Grid box inactive
            ELSE
               dp = 0.
            END IF
            
            IF ( dp < 0. ) THEN
               PRINT*,'depth of grid box is negative, ji,jj,jk',ji,jj,jk
               PRINT*,pxy(ji,jj),pxyz(ji,jj,jk-1),pxyz(ji,jj,jk)
               STOP
            END IF
            
            vol  (ji,jj,jk) =  dxdy(ji,jj) * dp / dg
            
            IF (txyz(ji,jj,ik) == undef2) THEN
               
               tem(ji,jj,jk) = undef
               sal(ji,jj,jk) = undef
               geo(ji,jj,jk) = undef
               vol(ji,jj,jk) = undef
               rho(ji,jj,jk) = undef
               
            END IF
            
            
            IF (txyz(ji,jj,ik) /= undef2 .AND. txyz(ip,jj,ik) /= undef2) THEN   
               uflux(ji,jj,jk) =  0.5 *( uxyz(ip,jj,ik) + uxyz(ji,jj,ik) ) / &
               &                  dg * dy(ji,jj) * dp 
            ELSE
               uflux(ji,jj,jk) = 0.
            END IF
            
            IF (jj /= JMT) THEN
               IF (txyz(ji,jj,ik) /= undef2 .AND. txyz(ji,jp,ik) /= undef2) THEN   
                  vflux(ji,jj,jk) =  0.5 *( vxyz(ji,jp,ik) + vxyz(ji,jj,ik) ) / &
                  &               dg * dx(ji,jj) * dp
               END IF
            ELSE
               vflux(ji,jj,jk) = 0.
            END IF
            
            !IF (ik /= 1) THEN
            !   IF (txyz(ji,jj,ik) /= undef2 .AND. txyz(ji,jj,ik-1) /= undef2) THEN
            !      wflux(ji,jj,jk) = 0.5 * ( wxyz(ji,jj,ik-1) + wxyz(ji,jj,ik) ) * &
            !      &                 dxdy(ji,jj) / dg
            !   END IF
            !ELSE
            !   wflux(ji,jj,jk) = 0.
            !END IF
            
            
         END DO
         
      END DO
      
      !! Zero mass flux at poles
      vflux(:,0,:)   = 0.
      vflux(:,JMT,:) = 0.
      
      !! Zero mass flux at TOA and surface
      wflux(:,:,0)   = 0.
      wflux(:,:,KM)  = 0.
      
      
   END DO
   
   DEALLOCATE( txyz, qxyz, uxyz, vxyz, wxyz, pxy, zxyz, pxyz )
   
   !!
   !! Calculate vertical mass flux
   !! Should we use this, or the data from MERRA?
   !!
   
   !wflux(:,:,0) = 0.
   !
   DO jj=1,JMT
      
      jm = jj-1
      
      DO jk=1,KM
         
         DO ji=1,IMT
            
            im = ji-1
            IF (im == 0) THEN
               im = IMT
            END IF
            
            wflux(ji,jj,jk) = wflux(ji,jj,jk-1) - &
            &                         ( uflux(ji,jj,jk) - uflux(im,jj,jk) + &
            &                           vflux(ji,jj,jk) - vflux(ji,jm,jk) )
            
         END DO
      END DO
   END DO
   
!$OMP END MASTER


   RETURN
   
   END SUBROUTINE
   
   
   
   
   
   SUBROUTINE get_data_gfdl(run)
   !!---------------------------------------------------------------------------
   !!
   !! Subroutine to get data from AM, the atmospheric model from GFDL.
   !! Data is stored every 6h with one file every year. 
   !!
   !!---------------------------------------------------------------------------
   
   REAL*4, ALLOCATABLE, DIMENSION(:,:) :: lat_bnds, lon_bnds
   CHARACTER, INTENT(IN) :: run*9
   
!$OMP MASTER
   
   dtstep = FLOAT(hourstep) * 3600.  
   
   !!
   !! Allocate dzt, dxdy and dx the first time step
   !! and reset them
   !!
   IF (iyear == yearstart .AND. &
   &    imon == monstart  .AND. &
   &    iday == daystart  .AND. &
   &   ihour == hourstart ) THEN
      
      lfirst = .true.
      
      !!
      !! Private arrays that will be kept throughout the integrations
      !!
      ALLOCATE( aa(0:KM), bb(0:KM), zxy(IMT,0:JMT) ) 
      
      dzt(:,:,:,:) = 0.
      dxdy(:,:) = 0.
      dx(:,:) = 0.
      dy(:,:) = 0.
      vlat2(:) = 0.
      aa(:) = 0.
      bb(:) = 0.
      
   ELSE
      
      lfirst = .false.
      
   END IF
   
   IF (imon == 1 .AND. iday == 1 .AND. ihour == 0) THEN
      istep = 0
   END IF
   istep = istep + 1
   
   !!
   !! Temporary arrays
   !!
   
   allocate ( txyz(IMT,0:JMT,KM), &
   &          qxyz(IMT,0:JMT,KM), &
   &          uxyz(IMT,0:JMT,KM), &
   &          vxyz(IMT,0:JMT,KM), &
   &          pxy (IMT,0:JMT),    &
   &          pu  (imt,jmt),      &
   &          pv  (imt,0:jmt),    &
   &          pt  (imt,jmt)  )
   
   geo_i(:,:,:) = 0.
   txyz(:,:,:) = 0.
   qxyz(:,:,:) = 0.
   uxyz(:,:,:) = 0.
   vxyz(:,:,:) = 0.
   pxy(:,:) = 0.
   
   !!
   !! Store previous time step and reset current
   !!
   dzt(:,:,:,1) = dzt(:,:,:,2)
   dzt(:,:,:,2) = 0.
   
   IF (lverbose) THEN
      PRINT*,' Day, month, step ',iday, imon, istep
   END IF
      
   IF (istep == 1 .OR. lfirst) THEN
      
      !!
      !! File names
      !! One file per year
      !!
      SELECT CASE (TRIM(run))
         
         CASE ('CM3-HISTR')
            string='0000010100-0000123123.nc'
            WRITE (string(1:4),   '(i4)') iyear
            WRITE (string(12:15), '(i4)') iyear
            dataprefix='6hrLev_GFDL-CM3_historical_r1i1p1_'//TRIM(string)
            string='/0000/'
            WRITE (string(2:5),       '(i4)') iyear
            uFile=TRIM(inDataDir)//TRIM(string)//'ua_'//TRIM(dataprefix)
            vFile=TRIM(inDataDir)//TRIM(string)//'va_'//TRIM(dataprefix)
            pFile=TRIM(inDataDir)//TRIM(string)//'ps_'//TRIM(dataprefix)
            tFile=TRIM(inDataDir)//TRIM(string)//'ta_'//TRIM(dataprefix)
            qFile=TRIM(inDataDir)//TRIM(string)//'hus_'//TRIM(dataprefix)
            zFile=TRIM(topoDir)//'orog_fx_GFDL-CM3_historical_r0i0p0.nc'
            topoFile=TRIM(topoDir)//'areacella_fx_GFDL-CM3_historical_r0i0p0.nc'
         
         CASE ('CM3-RCP85')
            string='0000010100-0000123123.nc'
            WRITE (string(1:4),   '(i4)') iyear
            WRITE (string(12:15), '(i4)') iyear
            dataprefix='6hrLev_GFDL-CM3_rcp85_r1i1p1_'//TRIM(string)
            string='/0000/'
            WRITE (string(2:5),       '(i4)') iyear
            uFile=TRIM(inDataDir)//TRIM(string)//'ua_'//TRIM(dataprefix)
            vFile=TRIM(inDataDir)//TRIM(string)//'va_'//TRIM(dataprefix)
            pFile=TRIM(inDataDir)//TRIM(string)//'ps_'//TRIM(dataprefix)
            tFile=TRIM(inDataDir)//TRIM(string)//'ta_'//TRIM(dataprefix)
            qFile=TRIM(inDataDir)//TRIM(string)//'hus_'//TRIM(dataprefix)
            zFile=TRIM(topoDir)//'orog_fx_GFDL-CM3_historical_r0i0p0.nc'
            topoFile=TRIM(topoDir)//'areacella_fx_GFDL-CM3_historical_r0i0p0.nc'
         
         CASE ('AMIP-CONT')
            string='0000010100-0000123123.nc'
            WRITE (string(1:4),   '(i4)') iyear
            WRITE (string(12:15), '(i4)') iyear
            dataprefix='6hrLev_GFDL-CM3_amip_r1i1p1_'//TRIM(string)
            string='0000'
            WRITE (string(1:4),       '(i4)') iyear
            uFile=TRIM(inDataDir)//TRIM(string)//'/ua_'//TRIM(dataprefix)
            vFile=TRIM(inDataDir)//TRIM(string)//'/va_'//TRIM(dataprefix)
            pFile=TRIM(inDataDir)//TRIM(string)//'/ps_'//TRIM(dataprefix)
            tFile=TRIM(inDataDir)//TRIM(string)//'/ta_'//TRIM(dataprefix)
            qFile=TRIM(inDataDir)//TRIM(string)//'/hus_'//TRIM(dataprefix)
            zFile=TRIM(topoDir)//'orog_fx_GFDL-CM3_historical_r0i0p0.nc'
            topoFile=TRIM(topoDir)//'areacella_fx_GFDL-CM3_historical_r0i0p0.nc'
      END SELECT
      
      
      !!
      !! Re-grid if necessary
      !!
      IF (tweak_zmean >= 1) THEN
         
         tmpFile = TRIM(tmpDataDir)//TRIM(prefix)//'tmp.nc'
         pFile2 = TRIM(tmpDataDir)//TRIM(prefix)//'ps.nc'
         uFile2 = TRIM(tmpDataDir)//TRIM(prefix)//'ua.nc'
         vFile2 = TRIM(tmpDataDir)//TRIM(prefix)//'va.nc'
         tFile2 = TRIM(tmpDataDir)//TRIM(prefix)//'ta.nc'
         qFile2 = TRIM(tmpDataDir)//TRIM(prefix)//'hus.nc'
         
         CALL cdo_interp(pFile,pFile2,topoFile,tmpFile,rstring,'ps',0,0,0)
         CALL cdo_interp(uFile,uFile2,topoFile,tmpFile,rstring,'ua',1,1,1)
         CALL cdo_interp(vFile,vFile2,topoFile,tmpFile,rstring,'va',1,0,0)
         CALL cdo_interp(tFile,tFile2,topoFile,tmpFile,rstring,'ta',1,0,0)
         CALL cdo_interp(qFile,qFile2,topoFile,tmpFile,rstring,'hus',1,0,0)
                  
         pFile = pFile2
         uFile = uFile2
         vFile = vFile2
         tFile = tFile2
         qFile = qFile2
      
      END IF
      
      !!
      !! Open the files and read the identifiers for the variables
      !!
      
      !! U file
      IF (lverbose) THEN
         PRINT*,uFile
      END IF
      ierr = NF90_OPEN( uFile, NF90_SHARE, id_ncu)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'ua', id_u) 
      CALL err(ierr)
      
      !! P file
      IF (lverbose) THEN
         PRINT*,pFile
      END IF
      ierr = NF90_OPEN( pFile, NF90_SHARE, id_ncp)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncp, 'ps', id_p) 
      CALL err(ierr)
      
      !! V file
      IF (lverbose) THEN
         PRINT*,vFile
      END IF
      ierr = NF90_OPEN( vFile, NF90_SHARE, id_ncv)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncv, 'va', id_v) 
      CALL err(ierr)
      
      !! T file
      IF (lverbose) THEN
         PRINT*,tFile
      END IF
      ierr = NF90_OPEN( tFile, NF90_SHARE, id_nct)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_nct, 'ta', id_t) 
      CALL err(ierr)
      
      !! Q file
      IF (lverbose) THEN
         PRINT*,qFile
      END IF
      ierr = NF90_OPEN( qFile, NF90_SHARE, id_ncq)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncq, 'hus', id_q) 
      CALL err(ierr)
      
      
      !!
      !! Check dimensions of uFile
      !!
      ierr = NF90_INQ_DIMID (id_ncu, 'lon', id_x)
      CALL err(ierr)
      ierr = NF90_INQ_DIMID (id_ncu, 'lat', id_y)
      CALL err(ierr)
      ierr = NF90_INQ_DIMID (id_ncu, 'lev', id_ml)
      CALL err(ierr)
      
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_x, LEN=IMT2)
      CALL err(ierr)
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_y, LEN=JMT2)
      CALL err(ierr)
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_ml, LEN=KM2)
      CALL err(ierr)
      
      !! Test so that IMT, JMT, KM that has been set agrees with the data
      IF (IMT2 /= IMT) THEN
         PRINT*,' Error: IMT in the data is not as you have set it! '
         PRINT*,IMT,IMT2
         STOP
      END IF
      IF (JMT2 /= JMT+1) THEN
         PRINT*,' Error: JMT in the data is not as you have set it! '
         PRINT*,JMT,JMT2
         STOP
      END IF
      IF (KM2 /= KM) THEN
         PRINT*,' Error: KM in the data is not as you have set it! '
         PRINT*,KM,KM2
         STOP
      END IF
      
   END IF
   
   IF (lfirst) THEN
      
      IF (lverbose) THEN
         PRINT*,' Reading lon, lat, lev, A(k), B(k) '
      END IF
      
      !!
      !! Read topography and grid info from separate files
      !!
      ierr = NF90_INQ_VARID (id_ncu, 'lon', id_lon)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'lat', id_lat)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'lev', id_lev)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'a_bnds', id_a)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'b_bnds', id_b)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'p0', id_p0)
      CALL err(ierr)
      
      ierr = NF90_GET_VAR (id_ncu, id_lon, vlon, start=[1], count=[IMT])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_lat, vlat2(0:jmt), start=[1], count=[JMT+1])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_lev, vlev, start=[1], count=[KM])
      CALL err(ierr)
      
      ierr = NF90_GET_VAR (id_ncu, id_a, aa(1:KM), start=[2,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_a, aa(0:KM-1), start=[1,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_b, bb(1:KM), start=[2,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_b, bb(0:KM-1), start=[1,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_p0, p0)
      CALL err(ierr)
      
      
      !!
      !! Define dy and calculate dx
      !!
      call calc_dxdy()
      
      
      !!
      !! Re-grid if necessary
      !!
      IF (tweak_zmean >= 1) THEN
         zFile2 = trim(tmpDataDir)//trim(prefix)//'z0.nc'
         CALL cdo_interp(zFile,zFile2,topoFile,tmpFile,rstring,'orog',0,0,0)
         zFile = zFile2
      END IF
      
      !!
      !! Read orography (Z file)
      !!
      ierr = NF90_OPEN( zFile, NF90_SHARE, id_ncz)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncz, 'orog', id_z) 
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncz, id_z, zxy(:,0:jmt),    start=[  1,     1], &
      &                                                   count=[IMT, JMT+1])
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncz)
      CALL err(ierr)
      IF (tweak_tend >=1) THEN
         CALL SYSTEM('rm '//TRIM(zFile))
      END IF
      
   END IF
   
   
   ierr = NF90_GET_VAR (id_ncu, id_u, uxyz(:,0:JMT,:), &
   &                                        start=[  1,   1,  1, istep],     &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncv, id_v, vxyz(:,0:JMT,:), &
   &                                        start=[  1,   1,  1, istep], &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_nct, id_t, txyz(:,0:JMT,:), &
   &                                        start=[  1,   1,  1, istep], &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncp, id_p, pxy(:,0:JMT),  &
   &                                        start=[  1,   1, istep], &
   &                                        count=[IMT, JMT+1,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncq, id_q, qxyz(:,0:JMT,:), &
   &                                        start=[  1,   1,  1, istep], &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   
         
   IF ( (imon == 12 .AND. iday == 31 .AND. ihour == 18) .OR. llast) THEN
      
      ierr = NF90_CLOSE(id_ncu)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncv)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncp)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_nct)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncq)
      CALL err(ierr)
      
      IF (tweak_zmean >= 1) THEN
         CALL SYSTEM('rm '//TRIM(pFile))
         CALL SYSTEM('rm '//TRIM(uFile))
         CALL SYSTEM('rm '//TRIM(vFile))
         CALL SYSTEM('rm '//TRIM(qFile))
         CALL SYSTEM('rm '//TRIM(tFile))
      END IF
      
   END IF
   
   
   !!
   !! A-grid -> C-grid
   !!
   call a2cgrid_u(pxy(:,:),pu(:,:))
   call a2cgrid_v(pxy(:,:),pv(:,:))
   call a2cgrid_t(pxy(:,:),pt(:,:))
   call a2cgrid_t(zxy(:,:),geo_i(:,:,km))
   
   rho_i(:,:,0) = aa(km) * p0 + bb(km) * pt(:,:)
   
   DO jk=1,KM
      
      !! Flip vertical axis
      il = km - jk + 1
      kb = km - jk
      ku = km - jk + 1
      
      da = aa(kb) - aa(ku)
      db = bb(kb) - bb(ku)
      rho_i(:,:,jk) = aa(kb) * p0 + bb(kb) * pt(:,:)
      dzt(:,:,jk,2) = da * p0 + db * pt(:,:)
      vol(:,:,jk)   = dzt(:,:,jk,2) * dxdy(:,:) / dg
      
      call a2cgrid_t(txyz(:,:,il),tem(:,:,jk))
      call a2cgrid_t(qxyz(:,:,il),sal(:,:,jk))
      
      call a2cgrid_u(uxyz(:,:,il),uxyz(:,:,il))
      uflux(:,:,jk) = uxyz(:,:,il) * dy(:,:) * ( da * p0 + db * pu(:,:)) / dg
      call a2cgrid_v(vxyz(:,:,il),vxyz(:,:,il))
      vflux(:,:,jk) = vxyz(:,:,il) * dx(:,:) * ( da * p0 + db * pv(:,:)) / dg
      
      
   END DO
   
   deallocate ( txyz, qxyz, uxyz, vxyz, pxy, pu, pv, pt )
   
   geo_i(:,:,km) = geo_i(:,:,km) * dg
   
   call int_geo()
   
   geo(:,:,1:km) = 0.5 * (geo_i(:,:,0:km-1) + geo_i(:,:,1:km))
   rho(:,:,1:km) = 0.5 * (rho_i(:,:,0:km-1) + rho_i(:,:,1:km))
   
   call calc_wflux()
   
   
!$OMP END MASTER


   RETURN
   
   END SUBROUTINE
   
   
   SUBROUTINE get_data_hadgem(run)
   !!---------------------------------------------------------------------------
   !!
   !! Subroutine to get data from HadGAM2, the atmospheric model from HadGEM2.
   !! Data is stored every 6h with one file covering Dec - Dec for one year. 
   !!
   !!---------------------------------------------------------------------------
   
   REAL*4, ALLOCATABLE, DIMENSION(:,:) :: lat_bnds, lon_bnds
   CHARACTER, INTENT(IN) :: run*9
   
!$OMP MASTER
   
   dtstep = FLOAT(hourstep) * 3600.  
   
   !!
   !! Allocate dzt, dxdy and dx the first time step
   !! and reset them
   !!
   IF (iyear == yearstart .AND. &
   &    imon == monstart  .AND. &
   &    iday == daystart  .AND. &
   &   ihour == hourstart ) THEN
      
      lfirst = .true.
      
      !!
      !! Private arrays that will be kept throughout the integrations
      !!
      ALLOCATE( dzt(IMT,JMT,KM,2), dxdy(IMT,JMT), dx(IMT,0:JMT), dy(IMT,JMT), &
      &         vlat2(0:JMT), aa(0:KM), bb(0:KM), zxy(IMT,0:JMT), steps(12,30,4) ) 
      
      dzt(:,:,:,:) = 0.
      dxdy(:,:) = 0.
      dx(:,:) = 0.
      dy(:,:) = 0.
      vlat2(:) = 0.
      aa(:) = 0.
      bb(:) = 0.
      
      istep = 119 !step at last day of year
      DO ji=1,12
         DO jj=1,30
            DO jk=1,4
               IF (ji == 12 .AND. jj == 1 .AND. jk == 2) THEN
                  istep = 0
               END IF
               istep = istep + 1
               steps(ji,jj,jk) = istep
            END DO
         END DO
      END DO
   
      
   ELSE
      
      lfirst = .false.
      
   END IF
   
   
   istep = steps(imon,iday,INT(ihour/6)+1)
   
   !!
   !! Temporary arrays
   !!
   
   ALLOCATE( txyz(IMT,0:JMT,KM), &
   &         qxyz(IMT,0:JMT,KM), &
   &         uxyz(IMT,0:JMT,KM), &
   &         vxyz(IMT,0:JMT,KM), &
   &         pxy(IMT,0:JMT) )
   
   txyz(:,:,:) = 0.
   qxyz(:,:,:) = 0.
   uxyz(:,:,:) = 0.
   vxyz(:,:,:) = 0.
   pxy(:,:) = 0.
   
   !!
   !! Store previous time step and reset current
   !!
   dzt(:,:,:,1) = dzt(:,:,:,2)
   dzt(:,:,:,2) = 0.
   
   IF (lverbose) THEN
      PRINT*,' Day, month, step ',iday, imon, istep
   END IF
      
   IF (istep == 1 .OR. lfirst) THEN
      
      !!
      !! File names
      !! One file per Dec-Dec year
      !!
      SELECT CASE (TRIM(run))
         CASE ('HAD-HISTR')
            string='0000120106-0000120100.nc'
            IF ( imon < 12 ) THEN
               WRITE (string(1:4),   '(i4)') iyear-1
               WRITE (string(12:15), '(i4)') iyear
            ELSE
               WRITE (string(1:4),   '(i4)') iyear
               WRITE (string(12:15), '(i4)') iyear+1
            END IF
            dataprefix='6hrLev_HadGEM2-ES_historical_r2i1p1_'//TRIM(string)
            string='/0000/'
            IF ( imon < 12 ) THEN
               WRITE (string(2:5),   '(i4)') iyear
            ELSE
               WRITE (string(2:5),   '(i4)') iyear+1
            END IF
            uFile=TRIM(inDataDir)//TRIM(string)//'ua_'//TRIM(dataprefix)
            vFile=TRIM(inDataDir)//TRIM(string)//'va_'//TRIM(dataprefix)
            pFile=TRIM(inDataDir)//TRIM(string)//'ps_'//TRIM(dataprefix)
            tFile=TRIM(inDataDir)//TRIM(string)//'ta_'//TRIM(dataprefix)
            qFile=TRIM(inDataDir)//TRIM(string)//'hus_'//TRIM(dataprefix)
      END SELECT
            
      IF (lverbose) THEN
         PRINT*,uFile
      END IF
      
      
      !!
      !! Open the files and read the identifiers for the variables
      !!
      
      !! U file
      ierr = NF90_OPEN( uFile, NF90_SHARE, id_ncu)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'ua', id_u) 
      CALL err(ierr)
      
      !! P file
      ierr = NF90_OPEN( pFile, NF90_SHARE, id_ncp)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncp, 'ps', id_p) 
      CALL err(ierr)
      
      !! V file
      ierr = NF90_OPEN( vFile, NF90_SHARE, id_ncv)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncv, 'va', id_v) 
      CALL err(ierr)
      
      !! T file
      ierr = NF90_OPEN( tFile, NF90_SHARE, id_nct)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_nct, 'ta', id_t) 
      CALL err(ierr)
      
      !! Q file
      ierr = NF90_OPEN( qFile, NF90_SHARE, id_ncq)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncq, 'hus', id_q) 
      CALL err(ierr)
      
      
      !!
      !! Check dimensions of tFile
      !!
      ierr = NF90_INQ_DIMID (id_nct, 'lon', id_x)
      CALL err(ierr)
      ierr = NF90_INQ_DIMID (id_nct, 'lat', id_y)
      CALL err(ierr)
      ierr = NF90_INQ_DIMID (id_nct, 'lev', id_ml)
      CALL err(ierr)
      
      ierr = NF90_INQUIRE_DIMENSION (id_nct, id_x, LEN=IMT2)
      CALL err(ierr)
      ierr = NF90_INQUIRE_DIMENSION (id_nct, id_y, LEN=JMT2)
      CALL err(ierr)
      ierr = NF90_INQUIRE_DIMENSION (id_nct, id_ml, LEN=KM2)
      CALL err(ierr)
      
      !! Test so that IMT, JMT, KM that has been set agrees with the data
      IF (IMT2 /= IMT) THEN
         PRINT*,' Error: IMT in the data is not as you have set it! '
         PRINT*,IMT,IMT2
         STOP
      END IF
      IF (JMT2 /= JMT+1) THEN
         PRINT*,' Error: JMT in the data is not as you have set it! '
         PRINT*,JMT,JMT2
         STOP
      END IF
      IF (KM2 /= KM) THEN
         PRINT*,' Error: KM in the data is not as you have set it! '
         PRINT*,KM,KM2
         STOP
      END IF
      
   END IF
   
   IF (lfirst) THEN
      
      IF (lverbose) THEN
         PRINT*,' Reading lon, lat, lev, A(k), B(k) '
      END IF
      
      !!
      !! Read orography
      !!
      ierr = NF90_INQ_VARID (id_nct, 'orog', id_z)
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_nct, id_z, zxy(1:IMT,0:JMT), start=[1,1], count=[IMT,JMT+1])
      CALL err(ierr)
      
      !!
      !! Read grid info
      !!
      ierr = NF90_INQ_VARID (id_nct, 'lon', id_lon)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_nct, 'lat', id_lat)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_nct, 'lev', id_lev)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_nct, 'lev_bnds', id_a)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_nct, 'b_bnds', id_b)
      CALL err(ierr)
      
      ierr = NF90_GET_VAR (id_nct, id_lon, vlon, start=[1], count=[IMT])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_nct, id_lat, vlat2(0:JMT), start=[1], count=[JMT+1])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_nct, id_lev, vlev, start=[1], count=[KM])
      CALL err(ierr)
      
      ierr = NF90_GET_VAR (id_nct, id_a, aa(1:KM), start=[2,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_nct, id_a, aa(0:KM-1), start=[1,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_nct, id_b, bb(1:KM), start=[2,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_nct, id_b, bb(0:KM-1), start=[1,1], count=[1,KM])
      CALL err(ierr)
      
      
      !! 
      !! dx, dy, dxdy 
      !! dx at v-points, dy at u-points and dxdy at T-points
      !!
      DO jj=0,JMT
         
         DO ji=1,IMT
            im = ji-1
            IF (im == 0) THEN
               im = IMT
            END IF
            
            dlon = vlon(ji) - vlon(im)
            IF (dlon < 0.) THEN
               dlon = dlon + 360.
            END IF
            IF (dlon > 360.) THEN
               dlon = dlon - 360.
            END IF
            dx(ji,jj) = ABS( deg * dlon * DCOS( vlat2(jj) * radian ) ) !dx at v-point
            
            IF(jj /= 0) THEN
               dlat = vlat2(jj) - vlat2(jj-1)
               dy(ji,jj) = ABS( deg * dlat ) ! dy at u-point
            END IF
            
         END DO
         
      END DO
            
      !!
      !! dxdy at T-points
      !!
      DO jj=1,JMT
         dxdy(:,jj) = 0.5 * ( dx(:,jj) + dx(:,jj-1) ) * dy(:,jj)
      END DO
      
      !!
      !! latitude at T-points
      !!      
      vlat(1:JMT) = 0.5 * (vlat2(1:JMT) + vlat2(0:JMT-1))
      
      IF (lverbose) THEN
         PRINT*,'dx(1,:)',dx(1,:)
         PRINT*,'dy(1,:)',dy(1,:)
         PRINT*,'dxdy(1,:)',dxdy(1,:)
      END IF
      
   
   END IF
   
   
   ierr = NF90_GET_VAR (id_ncu, id_u, uxyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep],     &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   
   ierr = NF90_GET_VAR (id_ncv, id_v, vxyz(:,0:JMT-1,:), &
   &                                        start=[  1,     1,  1, istep], &
   &                                        count=[IMT,   JMT, KM,     1])
   CALL err(ierr)
   
   ierr = NF90_GET_VAR (id_nct, id_t, txyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep], &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   
   ierr = NF90_GET_VAR (id_ncp, id_p, pxy(:,0:JMT),  &
   &                                        start=[  1,     1,  1, istep], &
   &                                        count=[IMT, JMT+1,  1,     1])
   CALL err(ierr)
   
   ierr = NF90_GET_VAR (id_ncq, id_q, qxyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep], &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
         
   IF ( (imon == 12 .AND. iday == 1 .AND. ihour == 0) .OR. llast) THEN
      
      ierr = NF90_CLOSE(id_ncu)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncv)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncp)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_nct)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncq)
      CALL err(ierr)
      
   END IF
      
   !!
   !! Calculate pressure
   !! jk is the level counted from the surface and up
   !! il is the level counted from TOA and down
   !! il-1 is then the grid box above and il+1 the one below. 
   !!
   
   rho_i(:,:,0) = pxy(:,:)
   
   DO jk = 1,KM
      
      DO jj = 0,JMT
         
         DO ji = 1,IMT
            
            !! layer height
            da = aa(jk) - aa(jk-1)
            db = bb(jk) - bb(jk-1)
            zz = zxy(ji,jj)
            dz = da + db * zz
            
            !! Virtual temperature at layer mid-point
            tt = txyz(ji,jj,jk)
            ss = qxyz(ji,jj,jk)
            tv = ( 1.0 + 0.61 * ss ) * tt
            
            !! pressure at upper interface
            rho_i(ji,jj,jk) = rho_i(ji,jj,jk-1) * exp( (-1.) * dg * dz / (Rd * tv) )
            
         END DO
      END DO
   END DO
   
   
   !!
   !! Calculate mass fluxes and re-grid data to C-grid
   !! 
   
   DO jk=1,KM
      
      il = KM-jk+1
         
      DO jj=1,JMT
         jm = jj-1
         
         DO ji=1,IMT
            
            im = ji - 1
            IF (im == 0) THEN
               im = IMT
            END IF         
            
            !! A-grid to C-grid
            tem  (ji,jj,jk) = 0.5 *  (txyz(ji,jj,il)   +  txyz(ji,jm,il))
            sal  (ji,jj,jk) = 0.5 *  (qxyz(ji,jj,il)   +  qxyz(ji,jm,il))
            rho  (ji,jj,jk) = 0.25* (rho_i(ji,jj,il)   + rho_i(ji,jm,il)   + &
                            &        rho_i(ji,jj,il-1) + rho_i(ji,jm,il-1))
            
            dzt  (ji,jj,jk,2)= 0.5* (rho_i(ji,jj,il-1) + rho_i(ji,jm,il-1))  - &
                            &  0.5* (rho_i(ji,jj,il)   + rho_i(ji,jm,il))
                              
            vol  (ji,jj,jk)  =  dzt(ji,jj,jk,2) * dxdy(ji,jj) / dg
            
            zz = 0.5 * (zxy(ji,jj) + zxy(ji,jm))
            da = 0.5 * (aa(il-1) + aa(il))
            db = 0.5 * (bb(il-1) + bb(il))
            geo  (ji,jj,jk) = (da + db * zz) * dg
                        
            !! mass fluxes
            dp = 0.5 * (rho_i(ji,jj,il-1) + rho_i(ji,jm,il-1)) - &
                 0.5 * (rho_i(ji,jj,il  ) + rho_i(ji,jm,il  ))
            
            uflux(ji,jj,jk)    = 0.5 * (uxyz(ji,jj,il) + uxyz(ji,jm,il)) * &
                               & dy(ji,jj) * dp / dg
            
            IF (jj /= JMT) THEN
               dp = rho_i(ji,jj,il-1) - rho_i(ji,jj,il)
               vflux(ji,jj,jk) = 0.5 * (vxyz(ji,jj,il) + vxyz(ji,jm,il)) * &
                               & dx(ji,jj) * dp / dg
            END IF
            
         END DO
      END DO
   END DO
   
   !!
   !! No fluxes at poles
   !!
   vflux(:,0,:) = 0.
   vflux(:,JMT,:) = 0.
   
   DEALLOCATE( txyz, qxyz, uxyz, vxyz )
   
   DEALLOCATE ( rho_i, pxy )
   
   !!
   !! Calculate vertical mass flux
   !!
   wflux(:,:,0) = 0.
   
   DO jj=1,JMT
      
      jm = jj-1
      
      DO jk=1,KM
         
         DO ji=1,IMT
            
            im = ji-1
            IF (im == 0) THEN
               im = IMT
            END IF
            
            wflux(ji,jj,jk) = wflux(ji,jj,jk-1) - &
            &                         ( uflux(ji,jj,jk) - uflux(im,jj,jk) + &
            &                           vflux(ji,jj,jk) - vflux(ji,jm,jk) + &
            &                          (dzt(ji,jj,jk,2) - dzt(ji,jj,jk,1))* &
            &                           dxdy(ji,jj) / dtstep / dg )
            
         END DO
      END DO
   END DO
!$OMP END MASTER


   RETURN
   
   END SUBROUTINE   
   
   
   SUBROUTINE get_data_ipsl(run)
   !!---------------------------------------------------------------------------
   !!
   !! Subroutine to get data from the IPSL CM5A model.
   !! Data is stored every 6h with one file every 10 years. 
   !!
   !!---------------------------------------------------------------------------
   
   CHARACTER, INTENT(IN) :: run*10
   
!$OMP MASTER
   
   dtstep = FLOAT(hourstep) * 3600.  
   
   !!
   !! Allocate dzt, dxdy and dx the first time step
   !! and reset them
   !!
   IF (iyear == yearstart .AND. &
   &    imon == monstart  .AND. &
   &    iday == daystart  .AND. &
   &   ihour == hourstart ) THEN
      
      lfirst = .true.
      
      !!
      !! Private arrays that will be kept throughout the integrations
      !!
      allocate ( aa(0:KM), bb(0:KM), zxy(IMT,0:JMT), steps(12,31,4) ) 
      
      dzt(:,:,:,:) = 0.
      dxdy(:,:) = 0.
      dx(:,:) = 0.
      dy(:,:) = 0.
      vlat2(:) = 0.
      aa(:) = 0.
      bb(:) = 0.
      
      istep = 0
      DO ji=1,12
         DO jj=1,idmax(ji,1850)
            DO jk=1,4
               istep = istep + 1
               steps(ji,jj,jk) = istep
            END DO
         END DO
      END DO
      
   ELSE
      
      lfirst = .false.
      
   END IF
   
   SELECT CASE (TRIM(run))
      
      CASE ('IPSL-HISTR')
         istep = INT( steps(imon,iday,INT(ihour/6)+1) + 1460 * MOD(iyear,10) )
      
      CASE ('IPSL-RCP85')
         IF (MOD(iyear,10) <= 5) THEN
            istep = INT( steps(imon,iday,INT(ihour/6)+1) + 1460 * & 
            & ( MOD(iyear,10) + 4 ) )
         ELSE IF (MOD(iyear,10) >= 6) THEN
            istep = INT( steps(imon,iday,INT(ihour/6)+1) + 1460 * & 
            & ( MOD(iyear,10) - 6 ) )
         END IF
         
   END SELECT
   
   !!
   !! Temporary arrays
   !!
   
   allocate ( txyz(IMT,0:JMT,KM), &
   &          qxyz(IMT,0:JMT,KM), &
   &          uxyz(IMT,0:JMT,KM), &
   &          vxyz(IMT,0:JMT,KM), &
   &          pxy(IMT,0:JMT),     &
   &          pt(imt,jmt),        &
   &          pu(imt,jmt),        &
   &          pv(imt,0:jmt)  )
   
   geo_i(:,:,:) = 0.
   txyz(:,:,:) = 0.
   qxyz(:,:,:) = 0.
   uxyz(:,:,:) = 0.
   vxyz(:,:,:) = 0.
   pxy(:,:) = 0.
   
   
   !!
   !! Store previous time step and reset current
   !!
   dzt(:,:,:,1) = dzt(:,:,:,2)
   dzt(:,:,:,2) = 0.
   
   !IF (lverbose) THEN
      PRINT*,' Day, month, step ',iday, imon, istep
   !END IF
      
   IF (istep == 1 .OR. lfirst) THEN
      
      !!
      !! File names
      !!
      SELECT CASE (TRIM(run))
         
         CASE ('IPSL-HISTR')
            string='000001010300-000012312100.nc'
            WRITE (string(1:4),   '(i4)') iyear - MOD(iyear,10)
            WRITE (string(14:17), '(i4)') iyear - MOD(iyear,10) + 9
            dataprefix='6hrLev_IPSL-CM5A-LR_historical_r1i1p1_'//TRIM(string)
            string = '/0000/'
            WRITE (string(2:5),   '(i4)') iyear - MOD(iyear,10)
            uFile=TRIM(inDataDir)//TRIM(string)//'/ua_'//TRIM(dataprefix)
            vFile=TRIM(inDataDir)//TRIM(string)//'/va_'//TRIM(dataprefix)
            tFile=TRIM(inDataDir)//TRIM(string)//'/ta_'//TRIM(dataprefix)
            qFile=TRIM(inDataDir)//TRIM(string)//'/hus_'//TRIM(dataprefix)
            zFile=TRIM(topoDir)//'orog_fx_IPSL-CM5A-LR_historical_r0i0p0.nc'
            topoFile=TRIM(topoDir)//'areacella_fx_IPSL-CM5A-LR_historical_r0i0p0.nc'
            
         CASE ('IPSL-RCP85')
            string='0000010103-0000123121.nc'
            IF (MOD(iyear,10) <= 5) THEN
               WRITE (string(1:4),   '(i4)') iyear - MOD(iyear,10) - 4
               WRITE (string(12:15), '(i4)') iyear - MOD(iyear,10) + 5
            ELSE IF (MOD(iyear,10) >= 6) THEN
               WRITE (string(1:4),   '(i4)') iyear - MOD(iyear,10) + 6
               WRITE (string(12:15), '(i4)') iyear - MOD(iyear,10) + 15
            END IF
            dataprefix='6hrLev_IPSL-CM5A-LR_rcp85_r1i1p1_'//TRIM(string)
            string = '/0000/'
            IF (MOD(iyear,10) <= 5) THEN
               WRITE (string(2:5),       '(i4)') iyear - MOD(iyear,10) - 4
            ELSE IF (MOD(iyear,10) >= 6) THEN
               WRITE (string(2:5),       '(i4)') iyear - MOD(iyear,10) + 6
            END IF
            uFile=TRIM(inDataDir)//TRIM(string)//'ua_'//TRIM(dataprefix)
            vFile=TRIM(inDataDir)//TRIM(string)//'va_'//TRIM(dataprefix)
            tFile=TRIM(inDataDir)//TRIM(string)//'ta_'//TRIM(dataprefix)
            qFile=TRIM(inDataDir)//TRIM(string)//'hus_'//TRIM(dataprefix)
            zFile=TRIM(topoDir)//'orog_fx_IPSL-CM5A-LR_historical_r0i0p0.nc'
            topoFile=TRIM(topoDir)//'areacella_fx_IPSL-CM5A-LR_historical_r0i0p0.nc'
      
      END SELECT
            
      
      !!
      !! Re-grid if necessary
      !! 
      IF (tweak_zmean >= 1) THEN
         
         tmpFile = TRIM(tmpDataDir)//TRIM(prefix)//'tmp.nc'
         uFile2 = TRIM(tmpDataDir)//TRIM(prefix)//'ua.nc'
         vFile2 = TRIM(tmpDataDir)//TRIM(prefix)//'va.nc'
         tFile2 = TRIM(tmpDataDir)//TRIM(prefix)//'ta.nc'
         qFile2 = TRIM(tmpDataDir)//TRIM(prefix)//'hus.nc'
         
         CALL cdo_interp( uFile, uFile2, topoFile, tmpFile, rstring, 'ua,ps',1,1,2)
         CALL cdo_interp( vFile, vFile2, topoFile, tmpFile, rstring, 'va',1,0,0)
         CALL cdo_interp( tFile, tFile2, topoFile, tmpFile, rstring, 'ta',1,0,0)
         CALL cdo_interp( qFile, qFile2, topoFile, tmpFile, rstring, 'hus',1,0,0)
                                           
         
         uFile = uFile2
         vFile = vFile2
         tFile = tFile2
         qFile = qFile2
      
      END IF
      
      
      !!
      !! Open the file and read the identifiers for the variables
      !!
      PRINT*,uFile
      ierr = NF90_OPEN( uFile, NF90_SHARE, id_ncu)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'ua', id_u) 
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'ps', id_p) 
      CALL err(ierr)
      
      IF (lverbose) THEN
         PRINT*,vFile
      END IF
      ierr = NF90_OPEN( vFile, NF90_SHARE, id_ncv)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncv, 'va', id_v) 
      CALL err(ierr)
      
      IF (lverbose) THEN
         PRINT*,tFile
      END IF
      ierr = NF90_OPEN( tFile, NF90_SHARE, id_nct)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_nct, 'ta', id_t) 
      CALL err(ierr)
      
      IF (lverbose) THEN
         PRINT*,qFile
      END IF
      ierr = NF90_OPEN( qFile, NF90_SHARE, id_ncq)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncq, 'hus', id_q) 
      CALL err(ierr)
      
      !!
      !! Check dimensions
      !!
      ierr = NF90_INQ_DIMID (id_ncu, 'lon', id_x)
      CALL err(ierr)
      ierr = NF90_INQ_DIMID (id_ncu, 'lat', id_y)
      CALL err(ierr)
      ierr = NF90_INQ_DIMID (id_ncu, 'lev', id_ml)
      CALL err(ierr)
      
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_x, LEN=IMT2)
      CALL err(ierr)
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_y, LEN=JMT2)
      CALL err(ierr)
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_ml, LEN=KM2)
      CALL err(ierr)
      
      !!
      !! Test so that IMT, JMT, KM that has been set agrees with the data
      !!
      IF (IMT2 /= IMT) THEN
         PRINT*,' Error: IMT in the data is not as you have set it! '
         PRINT*,IMT,IMT2
         STOP
      END IF
      IF (JMT2 /= JMT+1) THEN
         PRINT*,' Error: JMT in the data is not as you have set it! '
         PRINT*,JMT,JMT2
         STOP
      END IF
      IF (KM2 /= KM) THEN
         PRINT*,' Error: KM in the data is not as you have set it! '
         PRINT*,KM,KM2
         STOP
      END IF
      
   END IF
   
   IF (lfirst) THEN
      
      IF (lverbose) THEN
         PRINT*,' Reading lon, lat, lev, A(k), B(k) '
      END IF
      
      ierr = NF90_INQ_VARID (id_ncu, 'lon', id_lon)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'lat', id_lat)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'lev', id_lev)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'ap_bnds', id_a)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'b_bnds', id_b)
      CALL err(ierr)
      
      ierr = NF90_GET_VAR (id_ncu, id_lon, vlon, start=[1], count=[IMT])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_lat, vlat2(0:JMT), start=[1], count=[JMT+1])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_lev, vlev, start=[1], count=[KM])
      CALL err(ierr)
      
      ierr = NF90_GET_VAR (id_ncu, id_a, aa(1:KM), start=[2,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_a, aa(0:KM-1), start=[1,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_b, bb(1:KM), start=[2,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_b, bb(0:KM-1), start=[1,1], count=[1,KM])
      CALL err(ierr)
      
      
      !! 
      !! dx, dy, dxdy 
      !! dx at v-points, dy at u-points and dxdy at T-points
      !!
      call calc_dxdy()
      
      !!
      !! Read orography (Z file)
      !!
      if (tweak_zmean >= 1) then
         zFile2 = trim(tmpDataDir)//trim(prefix)//'z0.nc'
         call cdo_interp(zFile,zFile2,topoFile,tmpFile,rstring,'orog',0,0,0)
         zFile = zFile2
      end if
      
      IF (lverbose) THEN
         PRINT*,zFile
      END IF
      ierr = NF90_OPEN( zFile, NF90_SHARE, id_ncz)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncz, 'orog', id_z) 
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncz, id_z, zxy(:,0:JMT),  start=[  1,     1], &
      &                                                count=[IMT, JMT+1])
      CALL err(ierr)
      
      ierr = NF90_CLOSE(id_ncz)
      CALL err(ierr)
      
      if (tweak_zmean >= 1) then
         call system('rm '//trim(zFile))
      end if
      
   END IF
   
   ! u [m/s]
   ierr = NF90_GET_VAR (id_ncu, id_u, uxyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep],     &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   ! v [m/s]
   ierr = NF90_GET_VAR (id_ncv, id_v, vxyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep], &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   ! t [K]
   ierr = NF90_GET_VAR (id_nct, id_t, txyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep], &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   ! ps [Pa]
   ierr = NF90_GET_VAR (id_ncu, id_p, pxy(:,0:JMT),  &
   &                                        start=[  1,     1,  istep], &
   &                                        count=[IMT, JMT+1,      1])
   CALL err(ierr)
   ! q [kg/kg]
   ierr = NF90_GET_VAR (id_ncq, id_q, qxyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep], &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
         
   IF ( istep == 14600 .OR. llast) THEN
      
      ierr = NF90_CLOSE(id_ncu)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncv)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_nct)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncq)
      CALL err(ierr)
      
      if (tweak_zmean >= 1) then
         call system('rm '//trim(uFile))
         call system('rm '//trim(vFile))
         call system('rm '//trim(tFile))
         call system('rm '//trim(qFile))
      end if
      
   END IF
   
   !!
   !! A-grid -> C-grid 
   !!
   call a2cgrid_t(pxy(:,:),pt(:,:))
   call a2cgrid_t(zxy(:,:),geo_i(:,:,km))
   call a2cgrid_u(pxy(:,:),pu(:,:))
   call a2cgrid_v(pxy(:,:),pv(:,:))
   
   rho_i(:,:,0) = 1.
   
   DO jk=1,KM
      
      il = km - jk + 1
      kb = km - jk
      ku = km - jk + 1
      
      
      call a2cgrid_t(txyz(:,:,il),tem(:,:,jk))
      call a2cgrid_t(qxyz(:,:,il),sal(:,:,jk))
      
      da = aa(kb) - aa(ku)
      db = bb(kb) - bb(ku)
      
      rho_i(:,:,jk) = aa(kb) + bb(kb) * pt(:,:)
      dzt(:,:,jk,2) = da + db * pt(:,:)
      vol(:,:,jk) = dzt(:,:,jk,2) * dxdy(:,:) / dg
      
      !!
      !! u,v -> uflux,vflux
      !!
      call a2cgrid_u(uxyz(:,:,il),uxyz(:,:,il))
      uflux(:,:,jk) = uxyz(:,:,il) * (da + db * pu(:,:)) * dy(:,:) / dg
      call a2cgrid_v(vxyz(:,:,il),vxyz(:,:,il))
      vflux(:,:,jk) = vxyz(:,:,il) * (da + db * pv(:,:)) * dx(:,:) / dg
      
      
   END DO
   
   deallocate ( txyz, qxyz, uxyz, vxyz, pxy, pt, pu, pv )
   
   !!
   !! Calculate geopotential
   !!
   call int_geo()
   
   geo(:,:,1:km) = 0.5 * (geo_i(:,:,0:km-1) + geo_i(:,:,1:km))
   rho(:,:,1:km) = 0.5 * (rho_i(:,:,0:km-1) + rho_i(:,:,1:km))
   
   
   !!
   !! Calculate vertical mass flux
   !!
   call calc_wflux()
   
!$OMP END MASTER
   

   RETURN
   
   END SUBROUTINE
   
   
   SUBROUTINE get_data_ccsm(run)
   !!---------------------------------------------------------------------------
   !!
   !! Subroutine to get data from the CCSM4 CAM4 model.
   !! Data is stored every 6h with one file every 3 months. 
   !!
   !!---------------------------------------------------------------------------
   
   CHARACTER, INTENT(IN) :: run*5
   
!$OMP MASTER
   
   dtstep = FLOAT(hourstep) * 3600.  
   
   !!
   !! Allocate dzt, dxdy and dx the first time step
   !! and reset them
   !!
   IF (iyear == yearstart .AND. &
   &    imon == monstart  .AND. &
   &    iday == daystart  .AND. &
   &   ihour == hourstart ) THEN
      
      lfirst = .true.
      
      !!
      !! Private arrays that will be kept throughout the integrations
      !!
      allocate ( aa(0:KM), bb(0:KM), zxy(IMT,0:JMT), steps(12,31,4), steps2(12,31,4) ) 
      
      dzt(:,:,:,:) = 0.
      dxdy(:,:) = 0.
      dx(:,:) = 0.
      dy(:,:) = 0.
      vlat2(:) = 0.
      aa(:) = 0.
      bb(:) = 0.
      
      istep = 0
      istep2 = 0
      DO ji=1,12
         DO jj=1,idmax(ji,1850)
            DO jk=1,4
               IF (jk == 1 .AND. jj == 1 .AND. &
                    (ji == 1 .OR. ji == 4 .OR. ji == 7 .OR. ji == 10) ) THEN
                  istep = 0
               END IF
               istep = istep + 1
               istep2 = istep2 + 1
               steps(ji,jj,jk) = istep
               steps2(ji,jj,jk) = istep2
            END DO
         END DO
      END DO
      
   ELSE
      
      lfirst = .false.
      
   END IF
   
   istep  = steps(imon,iday,INT(ihour/6)+1)
   istep2 = steps2(imon,iday,INT(ihour/6)+1)
   
   PRINT*,iyear,imon,iday,ihour,istep,istep2
   
   !!
   !! Temporary arrays
   !!
   allocate ( txyz(IMT,0:JMT,KM), &
   &          qxyz(IMT,0:JMT,KM), &
   &          uxyz(IMT,0:JMT,KM), &
   &          vxyz(IMT,0:JMT,KM), &
   &          pxy(IMT,0:JMT),     &
   &          pt(imt,jmt),        & 
   &          pv(imt,0:jmt),      &
   &          pu(imt,jmt)  )
   
   geo_i(:,:,:) = 0.
   txyz(:,:,:) = 0.
   qxyz(:,:,:) = 0.
   uxyz(:,:,:) = 0.
   vxyz(:,:,:) = 0.
   pxy(:,:) = 0.
   
   
   !!
   !! Store previous time step and reset current
   !!
   dzt(:,:,:,1) = dzt(:,:,:,2)
   dzt(:,:,:,2) = 0.
   
   IF (lverbose) THEN
      PRINT*,' Day, month, u-step ps-step ',iday, imon, istep, istep2
   END IF
      
   IF (istep == 1 .OR. lfirst) THEN
      
      !!
      !! File names
      !!
      SELECT CASE (TRIM(run))
         CASE ('HISTR')
            
            string='0000000100-0000000018.nc'
            WRITE (string(1:4),   '(i4.4)') iyear
            WRITE (string(5:6),   '(i2.2)') imon
            WRITE (string(12:15), '(i4.4)') iyear
            WRITE (string(16:17), '(i2.2)') imon+2
            WRITE (string(7:8),   '(i2.2)') iday
            WRITE (string(18:19), '(i2.2)') idmax(imon+2,1850)
            dataprefix='6hrLev_CCSM4_historical_r6i1p1_'//TRIM(string)
            string='/0000/'
            WRITE (string(2:5),'(i4.4)') iyear
            uFile=TRIM(inDataDir)//TRIM(string)//'ua_'//TRIM(dataprefix)
            vFile=TRIM(inDataDir)//TRIM(string)//'va_'//TRIM(dataprefix)
            tFile=TRIM(inDataDir)//TRIM(string)//'ta_'//TRIM(dataprefix)
            qFile=TRIM(inDataDir)//TRIM(string)//'hus_'//TRIM(dataprefix)
            zFile=TRIM(topoDir)//'orog_fx_CCSM4_historical_r0i0p0.nc'
            topoFile=TRIM(topoDir)//'areacella_fx_CCSM4_historical_r0i0p0.nc'
            
            IF ( istep2 == 1 ) THEN
               string='0000010100-0000123118.nc'
               WRITE (string(1:4),   '(i4)') iyear
               WRITE (string(12:15), '(i4)') iyear   
               dataprefix='6hrLev_CCSM4_historical_r6i1p1_'//TRIM(string)
               string='/0000/'
               WRITE (string(2:5),'(i4.4)') iyear
               pFile=TRIM(inDataDir)//TRIM(string)//'ps_'//TRIM(dataprefix)
            END IF
            
         CASE ('RCP85')
            
            string='0000000100-0000000018.nc'
            WRITE (string(1:4),   '(i4.4)') iyear
            WRITE (string(5:6),   '(i2.2)') imon
            WRITE (string(12:15), '(i4.4)') iyear
            WRITE (string(16:17), '(i2.2)') imon+2
            WRITE (string(7:8),   '(i2.2)') iday
            WRITE (string(18:19), '(i2.2)') idmax(imon+2,1850)
            dataprefix='6hrLev_CCSM4_rcp85_r6i1p1_'//TRIM(string)
            string='/0000/'
            WRITE (string(2:5),'(i4.4)') iyear
            uFile=TRIM(inDataDir)//TRIM(string)//'ua_'//TRIM(dataprefix)
            vFile=TRIM(inDataDir)//TRIM(string)//'va_'//TRIM(dataprefix)
            tFile=TRIM(inDataDir)//TRIM(string)//'ta_'//TRIM(dataprefix)
            qFile=TRIM(inDataDir)//TRIM(string)//'hus_'//TRIM(dataprefix)
            zFile=TRIM(topoDir)//'orog_fx_CCSM4_historical_r0i0p0.nc'
            topoFile=TRIM(topoDir)//'areacella_fx_CCSM4_historical_r0i0p0.nc'
            
            IF ( istep2 == 1 .OR. lfirst ) THEN
               string='0000010100-0000123118.nc'
               WRITE (string(1:4),   '(i4)') iyear
               WRITE (string(12:15), '(i4)') iyear   
               dataprefix='6hrLev_CCSM4_rcp85_r6i1p1_'//TRIM(string)
               string='/0000/'
               WRITE (string(2:5),'(i4.4)') iyear
               pFile=TRIM(inDataDir)//TRIM(string)//'ps_'//TRIM(dataprefix)
            END IF
            
      END SELECT
            
      
      !!
      !! Re-grid data if necessary
      !!
      IF (tweak_zmean >= 1) THEN
         
         tmpFile = TRIM(tmpDataDir)//TRIM(prefix)//'tmp.nc'
         pFile2 = TRIM(tmpDataDir)//TRIM(prefix)//'ps.nc'
         uFile2 = TRIM(tmpDataDir)//TRIM(prefix)//'ua.nc'
         vFile2 = TRIM(tmpDataDir)//TRIM(prefix)//'va.nc'
         tFile2 = TRIM(tmpDataDir)//TRIM(prefix)//'ta.nc'
         qFile2 = TRIM(tmpDataDir)//TRIM(prefix)//'hus.nc'
         
         call cdo_interp( uFile, uFile2, topoFile, tmpFile, rstring, 'ua',1,1,1)
         call cdo_interp( vFile, vFile2, topoFile, tmpFile, rstring, 'va',1,0,0)
         call cdo_interp( tFile, tFile2, topoFile, tmpFile, rstring, 'ta',1,0,0)
         call cdo_interp( qFile, qFile2, topoFile, tmpFile, rstring, 'hus',1,0,0)
                                           
         IF ( istep2 == 1 .OR. lfirst) THEN
            
            call cdo_interp( pFile, pFile2, topoFile, tmpFile, rstring, 'ps',0,0,0)
            
         END IF
         
         pFile = pFile2
         uFile = uFile2
         vFile = vFile2
         tFile = tFile2
         qFile = qFile2
      
      END IF
      
      
      !!
      !! Open the file and read the identifiers for the variables
      !!
      IF (lverbose) THEN
         PRINT*,uFile
      END IF
      ierr = NF90_OPEN( uFile, NF90_SHARE, id_ncu)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'ua', id_u) 
      CALL err(ierr)
      
      IF (istep2 == 1) THEN
         IF (lverbose) THEN
            PRINT*,pFile
         END IF
         ierr = NF90_OPEN( pFile, NF90_SHARE, id_ncp)
         CALL err(ierr)
         ierr = NF90_INQ_VARID (id_ncp, 'ps', id_p) 
         CALL err(ierr)
      END IF
      
      IF (lverbose) THEN
         PRINT*,vFile
      END IF
      ierr = NF90_OPEN( vFile, NF90_SHARE, id_ncv)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncv, 'va', id_v) 
      CALL err(ierr)
      
      IF (lverbose) THEN
         PRINT*,tFile
      END IF
      ierr = NF90_OPEN( tFile, NF90_SHARE, id_nct)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_nct, 'ta', id_t) 
      CALL err(ierr)
      
      IF (lverbose) THEN
         PRINT*,pFile
      END IF
      ierr = NF90_OPEN( qFile, NF90_SHARE, id_ncq)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncq, 'hus', id_q) 
      CALL err(ierr)
      
      !!
      !! Check dimensions
      !!
      ierr = NF90_INQ_DIMID (id_ncu, 'lon', id_x)
      CALL err(ierr)
      ierr = NF90_INQ_DIMID (id_ncu, 'lat', id_y)
      CALL err(ierr)
      ierr = NF90_INQ_DIMID (id_ncu, 'lev', id_ml)
      CALL err(ierr)
      
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_x, LEN=IMT2)
      CALL err(ierr)
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_y, LEN=JMT2)
      CALL err(ierr)
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_ml, LEN=KM2)
      CALL err(ierr)
      
      !!
      !! Test so that IMT, JMT, KM that has been set agrees with the data
      !!
      IF (IMT2 /= IMT) THEN
         PRINT*,' Error: IMT in the data is not as you have set it! '
         PRINT*,IMT,IMT2
         STOP
      END IF
      IF (JMT2 /= JMT+1) THEN
         PRINT*,' Error: JMT in the data is not as you have set it! '
         PRINT*,JMT,JMT2
         STOP
      END IF
      IF (KM2 /= KM) THEN
         PRINT*,' Error: KM in the data is not as you have set it! '
         PRINT*,KM,KM2
         STOP
      END IF
      
   END IF
   
   IF (lfirst) THEN
      
      IF (lverbose) THEN
         PRINT*,' Reading lon, lat, lev, A(k), B(k) '
      END IF
      
      ierr = NF90_INQ_VARID (id_ncu, 'lon', id_lon)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'lat', id_lat)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'lev', id_lev)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'a_bnds', id_a)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'b_bnds', id_b)
      CALL err(ierr)
      
      ierr = NF90_GET_VAR (id_ncu, id_lon, vlon, start=[1], count=[IMT])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_lat, vlat2(0:JMT), start=[1], count=[JMT+1])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_lev, vlev, start=[1], count=[KM])
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'p0', id_p0)
      CALL err(ierr)
      
      ierr = NF90_GET_VAR (id_ncu, id_a, aa(1:KM), start=[2,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_a, aa(0:KM-1), start=[1,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_b, bb(1:KM), start=[2,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_b, bb(0:KM-1), start=[1,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_p0, p0)
      CALL err(ierr)
      
      p0 = p0 * 100.
      
      !! 
      !! dx, dy, dxdy 
      !! dx at v-points, dy at u-points and dxdy at T-points
      !!
      call calc_dxdy()
      
      !!
      !! Read orography (Z file)
      !!
      if (tweak_zmean >= 1) then
         zFile2 = trim(tmpDataDir)//trim(prefix)//'.z0.nc'
         call cdo_interp(zFile,zFile2,topoFile,tmpFile,rstring,'orog',0,0,0)
         zFile = zFile2
      end if
      ierr = NF90_OPEN( zFile, NF90_SHARE, id_ncz)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncz, 'orog', id_z) 
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncz, id_z, zxy(:,0:JMT), start=[  1,     1], &
      &                                                count=[IMT, JMT+1])
      CALL err(ierr)
      
      ierr = NF90_CLOSE(id_ncz)
      CALL err(ierr)
      
      if (tweak_zmean >= 1) then
         call system('rm '//trim(zFile))
      end if
      
   END IF
   
   ! u [m/s]
   ierr = NF90_GET_VAR (id_ncu, id_u, uxyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep],     &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   ! v [m/s]
   ierr = NF90_GET_VAR (id_ncv, id_v, vxyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep], &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   ! t [K]
   ierr = NF90_GET_VAR (id_nct, id_t, txyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep], &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   ! ps [Pa]
   ierr = NF90_GET_VAR (id_ncp, id_p, pxy(:,0:JMT),  &
   &                                        start=[  1,     1,  istep2], &
   &                                        count=[IMT, JMT+1,      1])
   CALL err(ierr)
   
   ! q [kg/kg]
   ierr = NF90_GET_VAR (id_ncq, id_q, qxyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep], &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
         
   IF ( (mod(imon,3) == 0 .AND. iday == idmax(imon,1850) .AND. ihour == 18) .OR. llast) THEN
      
      ierr = NF90_CLOSE(id_ncu)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncv)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_nct)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncq)
      CALL err(ierr)
      
      CALL system('rm -f '//TRIM(qFile))
      CALL system('rm -f '//TRIM(tFile))
      CALL system('rm -f '//TRIM(uFile))
      CALL system('rm -f '//TRIM(vFile))
      
   END IF
   
   IF ( (imon == 12 .AND. iday == 31 .AND. ihour == 18) .OR. llast ) THEN
      ierr = NF90_CLOSE(id_ncp)
      CALL err(ierr)
      CALL system('rm -f '//trim(pFile))
   END IF
   
   
   !!
   !! A grid to C-grid
   !!
   call a2cgrid_t(pxy(:,:),pt(:,:))
   call a2cgrid_u(pxy(:,:),pu(:,:))
   call a2cgrid_v(pxy(:,:),pv(:,:))
   call a2cgrid_t(zxy(:,:),geo_i(:,:,km))
   
   rho_i(:,:,0) = aa(km) * p0 + bb(km) * pt(:,:)
   
   DO jk=1,KM
      
      il = km - jk + 1
      ku = km - jk + 1
      kb = km - jk
      
      !!
      !! A-grid -> C-grid 
      !!
      call a2cgrid_t(txyz(:,:,il),tem(:,:,jk))
      call a2cgrid_t(qxyz(:,:,il),sal(:,:,jk))
      
      da = aa(kb) - aa(ku)
      db = bb(kb) - bb(ku)
      rho_i(:,:,jk) = aa(kb) * p0 + bb(kb) * pt(:,:)
      dzt(:,:,jk,2) = da * p0 + db * pt(:,:)
      vol(:,:,jk) = dzt(:,:,jk,2) * dxdy(:,:) / dg
      
      !!
      !! u,v -> uflux,vflux
      !!
      call a2cgrid_u(uxyz(:,:,il),uxyz(:,:,il))
      uflux(:,:,jk) = uxyz(:,:,il) * (da * p0 + db * pu(:,:)) * dy(:,:) / dg
      call a2cgrid_v(vxyz(:,:,il),vxyz(:,:,il))
      vflux(:,:,jk) = vxyz(:,:,il) * (da * p0 + db * pv(:,:)) * dx(:,:) / dg
      
      
   END DO
   
   deallocate ( txyz, qxyz, uxyz, vxyz, pxy, pt, pu, pv )
   
   
   !!
   !! Calculate geopotential
   !!
   call int_geo()
   
   geo(:,:,1:km) = 0.5 * (geo_i(:,:,1:km) + geo_i(:,:,0:km-1))
   rho(:,:,1:km) = 0.5 * (rho_i(:,:,1:km) + rho_i(:,:,0:km-1))
   
   !!
   !! Calculate vertical mass flux
   !!
   call calc_wflux()
   
!$OMP END MASTER
   

   RETURN
   
   END SUBROUTINE
   
   
   
   SUBROUTINE get_data_canesm(run)
   !!---------------------------------------------------------------------------
   !!
   !! Subroutine to get data from CanAM, the atmospheric model from CanESM.
   !! Data is stored every 6h with one file every year. 
   !!
   !!---------------------------------------------------------------------------
   
   CHARACTER, INTENT(IN) :: run*10
   
!$OMP MASTER
   
   dtstep = FLOAT(hourstep) * 3600.  
   
   !!
   !! Allocate dzt, dxdy and dx the first time step
   !! and reset them
   !!
   IF (iyear == yearstart .AND. &
   &    imon == monstart  .AND. &
   &    iday == daystart  .AND. &
   &   ihour == hourstart ) THEN
      
      lfirst = .true.
      
      !!
      !! Private arrays that will be kept throughout the integrations
      !!
      ALLOCATE( aa(0:KM), bb(0:KM), zxy(IMT,0:JMT), steps(12,31,4) ) 
      
      dzt(:,:,:,:) = 0.
      dxdy(:,:) = 0.
      dx(:,:) = 0.
      dy(:,:) = 0.
      vlat2(:) = 0.
      aa(:) = 0.
      bb(:) = 0.
      zxy(:,:) = 0.
      
      istep = 0
      DO ji=1,12
        DO jj=1,idmax(ji,1850)
           DO jk=1,4
              istep = istep + 1
              steps(ji,jj,jk) = istep
           END DO
        END DO
      END DO
      
   ELSE
      
      lfirst = .false.
      
   END IF
   
   istep = steps(imon,iday,INT(ihour/6)+1)
   
   !!
   !! Temporary arrays
   !!
   ALLOCATE( txyz(IMT,0:JMT,KM), &
   &         qxyz(IMT,0:JMT,KM), &
   &         uxyz(IMT,0:JMT,KM), &
   &         vxyz(IMT,0:JMT,KM), &
   &         pxy(IMT,0:JMT),     &
   &         pt(imt,jmt),      & 
   &         pu(imt,jmt),      & 
   &         pv(imt,0:jmt)  )
   
   geo_i(:,:,:) = 0.
   txyz(:,:,:) = 0.
   qxyz(:,:,:) = 0.
   uxyz(:,:,:) = 0.
   vxyz(:,:,:) = 0.
   pxy(:,:) = 0.
   
   !!
   !! Store previous time step and reset current
   !!
   dzt(:,:,:,1) = dzt(:,:,:,2)
   dzt(:,:,:,2) = 0.
   
   IF (lverbose) THEN
      PRINT*,' Day, month, step ',iday, imon, istep
   END IF
      
   IF (istep == 1 .OR. lfirst) THEN
      
      !!
      !! File names
      !! One file per year
      !!
      SELECT CASE (TRIM(run))
         
         CASE ('ESM2-HISTR')
            string='000001010000-000012311800.nc'
            WRITE (string(1:4),   '(i4)') iyear
            WRITE (string(14:17), '(i4)') iyear
            dataprefix='6hrLev_CanESM2_historical_r1i1p1_'//TRIM(string)
            string='/0000/'
            WRITE (string(2:5),       '(i4)') iyear
            uFile=TRIM(inDataDir)//TRIM(string)//'ua_'//TRIM(dataprefix)
            vFile=TRIM(inDataDir)//TRIM(string)//'va_'//TRIM(dataprefix)
            pFile=TRIM(inDataDir)//TRIM(string)//'ps_'//TRIM(dataprefix)
            tFile=TRIM(inDataDir)//TRIM(string)//'ta_'//TRIM(dataprefix)
            qFile=TRIM(inDataDir)//TRIM(string)//'hus_'//TRIM(dataprefix)
            zFile=TRIM(topoDir)//'orog_fx_CanESM2_historical_r0i0p0.nc'
            topoFile=TRIM(topoDir)//'areacella_fx_CanESM2_historical_r0i0p0.nc'
         
         CASE ('ESM2-RCP85')
            string='000001010000-000012311800.nc'
            WRITE (string(1:4),   '(i4)') iyear
            WRITE (string(14:17), '(i4)') iyear
            dataprefix='6hrLev_CanESM2_rcp85_r1i1p1_'//TRIM(string)
            string='/0000/'
            WRITE (string(2:5),       '(i4)') iyear
            uFile=TRIM(inDataDir)//TRIM(string)//'ua_'//TRIM(dataprefix)
            vFile=TRIM(inDataDir)//TRIM(string)//'va_'//TRIM(dataprefix)
            pFile=TRIM(inDataDir)//TRIM(string)//'ps_'//TRIM(dataprefix)
            tFile=TRIM(inDataDir)//TRIM(string)//'ta_'//TRIM(dataprefix)
            qFile=TRIM(inDataDir)//TRIM(string)//'hus_'//TRIM(dataprefix)
            zFile=TRIM(topoDir)//'orog_fx_CanESM2_rcp85_r0i0p0.nc'
            topoFile=TRIM(topoDir)//'areacella_fx_CanESM2_rcp85_r0i0p0.nc'
         
      END SELECT
      
      
      !!
      !! Re-grid
      !!
      if (tweak_zmean >= 1) then
         
         tmpFile = trim(tmpDataDir)//trim(prefix)//'tmp.nc'
         uFile2  = trim(tmpDataDir)//trim(prefix)//'ua.nc'
         vFile2  = trim(tmpDataDir)//trim(prefix)//'va.nc'
         tFile2  = trim(tmpDataDir)//trim(prefix)//'ta.nc'
         qFile2  = trim(tmpDataDir)//trim(prefix)//'hus.nc'
         pFile2  = trim(tmpDataDir)//trim(prefix)//'ps.nc'
         
         call cdo_interp(uFile,uFile2,topoFile,tmpFile,rstring,'ua',1,1,2)
         call cdo_interp(vFile,vFile2,topoFile,tmpFile,rstring,'va',1,0,0)
         call cdo_interp(qFile,qFile2,topoFile,tmpFile,rstring,'hus',1,0,0)
         call cdo_interp(tFile,tFile2,topoFile,tmpFile,rstring,'ta',1,0,0)
         call cdo_interp(pFile,pFile2,topoFile,tmpFile,rstring,'ps',0,0,0)
         
         uFile = uFile2
         vFile = vFile2
         tFile = tFile2
         qFile = qFile2
         pFile = pFile2
         
      end if
            
      !!
      !! Open the files and read the identifiers for the variables
      !!
      
      !! U file
      IF (lverbose) THEN
         PRINT*,uFile
      END IF
      ierr = NF90_OPEN( uFile, NF90_SHARE, id_ncu)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'ua', id_u) 
      CALL err(ierr)
      
      !! P file
      IF (lverbose) THEN
         PRINT*,pFile
      END IF
      ierr = NF90_OPEN( pFile, NF90_SHARE, id_ncp)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncp, 'ps', id_p) 
      CALL err(ierr)
      
      !! V file
      IF (lverbose) THEN
         PRINT*,vFile
      END IF
      ierr = NF90_OPEN( vFile, NF90_SHARE, id_ncv)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncv, 'va', id_v) 
      CALL err(ierr)
      
      !! T file
      IF (lverbose) THEN
         PRINT*,tFile
      END IF
      ierr = NF90_OPEN( tFile, NF90_SHARE, id_nct)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_nct, 'ta', id_t) 
      CALL err(ierr)
      
      !! Q file
      IF (lverbose) THEN
         PRINT*,qFile
      END IF
      ierr = NF90_OPEN( qFile, NF90_SHARE, id_ncq)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncq, 'hus', id_q) 
      CALL err(ierr)
      
      
      !!
      !! Check dimensions of uFile
      !!
      ierr = NF90_INQ_DIMID (id_ncu, 'lon', id_x)
      CALL err(ierr)
      ierr = NF90_INQ_DIMID (id_ncu, 'lat', id_y)
      CALL err(ierr)
      ierr = NF90_INQ_DIMID (id_ncu, 'lev', id_ml)
      CALL err(ierr)
      
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_x, LEN=IMT2)
      CALL err(ierr)
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_y, LEN=JMT2)
      CALL err(ierr)
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_ml, LEN=KM2)
      CALL err(ierr)
      
      !! Test so that IMT, JMT, KM that has been set agrees with the data
      IF (IMT2 /= IMT) THEN
         PRINT*,' Error: IMT in the data is not as you have set it! '
         PRINT*,IMT,IMT2
         STOP
      END IF
      IF (JMT2 /= JMT+1) THEN
         PRINT*,' Error: JMT in the data is not as you have set it! '
         PRINT*,JMT,JMT2
         STOP
      END IF
      IF (KM2 /= KM) THEN
         PRINT*,' Error: KM in the data is not as you have set it! '
         PRINT*,KM,KM2
         STOP
      END IF
      
   END IF
   
   IF (lfirst) THEN
      
      IF (lverbose) THEN
         PRINT*,' Reading lon, lat, lev, A(k), B(k) '
      END IF
      
      !!
      !! Read orography
      !! 
      if (tweak_zmean >= 1) then
         
         zFile2 = trim(tmpDataDir)//trim(prefix)//'z0.nc'
         call cdo_interp(zFile,zFile2,topoFile,tmpFile,rstring,'orog',0,0,0)
         zFile = zFile2
         
      end if
      
      ierr = NF90_OPEN( zFile, NF90_SHARE, id_ncz)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncz, 'orog', id_z) 
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncz, id_z, zxy(:,0:jmt), start=[  1,     1], &
      &                                            count=[IMT, JMT+1])
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncz)
      CALL err(ierr)
      
      if (tweak_zmean >= 1) then
         
         call system('rm '//trim(zFile))
         
      end if
      
      !!
      !! Read other grid info from uFile
      !!
      ierr = NF90_INQ_VARID (id_ncu, 'lon', id_lon)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'lat', id_lat)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'lev', id_lev)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'ap_bnds', id_a)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'b_bnds', id_b)
      CALL err(ierr)
      
      ierr = NF90_GET_VAR (id_ncu, id_lon, vlon, start=[1], count=[IMT])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_lat, vlat2, start=[1], count=[JMT+1])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_lev, vlev, start=[1], count=[KM])
      CALL err(ierr)
      
      ierr = NF90_GET_VAR (id_ncu, id_a, aa(1:KM), start=[2,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_a, aa(0:KM-1), start=[1,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_b, bb(1:KM), start=[2,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_b, bb(0:KM-1), start=[1,1], count=[1,KM])
      CALL err(ierr)
      
      !!
      !! Define dy and calculate dx
      !!
      call calc_dxdy()
      
      
   END IF
   
   ierr = NF90_GET_VAR (id_ncu, id_u, uxyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep],     &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncv, id_v, vxyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep], &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_nct, id_t, txyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep], &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncp, id_p, pxy(:,0:JMT),  &
   &                                        start=[  1,     1, istep], &
   &                                        count=[IMT, JMT+1,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncq, id_q, qxyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep], &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   
         
   IF ( (imon == 12 .AND. iday == 31 .AND. ihour == 18) .OR. llast) THEN
      
      ierr = NF90_CLOSE(id_ncu)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncv)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncp)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_nct)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncq)
      CALL err(ierr)
      
      if (tweak_zmean >= 1) then
         call system('rm '//trim(uFile))
         call system('rm '//trim(vFile))
         call system('rm '//trim(qFile))
         call system('rm '//trim(tFile))
         call system('rm '//trim(pFile))
      end if
      
   END IF
   
   
   !!
   !! A-grid -> C-grid
   !!
   call a2cgrid_t(pxy(:,:),pt(:,:))
   call a2cgrid_u(pxy(:,:),pu(:,:))
   call a2cgrid_v(pxy(:,:),pv(:,:))
   call a2cgrid_t(zxy(:,:),geo_i(:,:,km))
   
   rho_i(:,:,0) = aa(km) + bb(km) * pt(:,:)
   
   
   DO jk=1,KM
      
      il = km - jk + 1
      ku = km - jk + 1
      kb = km - jk
      
      !!
      !! A-grid -> C-grid 
      !! u,v -> uflux,vflux
      !!
      call a2cgrid_t(txyz(:,:,il),tem(:,:,jk))
      call a2cgrid_t(qxyz(:,:,il),sal(:,:,jk))
      
      rho_i(:,:,jk) = aa(kb) + bb(kb) * pt(:,:)
      
      da = aa(kb) - aa(ku)
      db = bb(kb) - bb(ku)
      dzt(:,:,jk,2) = da + db * pt(:,:)
      vol(:,:,jk) = dzt(:,:,jk,2) * dxdy(:,:) / dg 
      
      call a2cgrid_u(uxyz(:,:,il),uxyz(:,:,il))
      uflux(:,:,jk) = uxyz(:,:,il) * dy(:,:) * (da + db * pu(:,:)) / dg
      call a2cgrid_v(vxyz(:,:,il),vxyz(:,:,il))
      vflux(:,:,jk) = vxyz(:,:,il) * dx(:,:) * (da + db * pv(:,:)) / dg
      
   END DO
   
   
   DEALLOCATE( txyz, qxyz, uxyz, vxyz, pxy, pv, pu, pt )
   
   !!
   !! Geopotential at surface
   !!
   geo_i(:,:,KM) = geo_i(:,:,km) * dg
   
   !!
   !! Calculate geopotential
   !!
   call int_geo()
   
   !!
   !! Calculate DSE, MSE
   !!
   geo(:,:,1:km) = 0.5 * (geo_i(:,:,1:km) + geo_i(:,:,0:km-1))
   rho(:,:,1:km) = 0.5 * (rho_i(:,:,1:km) + rho_i(:,:,0:km-1))
   
   
   !!
   !! Calculate vertical mass flux
   !!
   call calc_wflux()
   
!$OMP END MASTER


   RETURN
   
   END SUBROUTINE
   
   
   
   SUBROUTINE get_data_noresm(run)
   !!---------------------------------------------------------------------------
   !!
   !! Subroutine to get data from the CAM-Oslo model, part of NorESM1.
   !! Data is stored every 6h with one file every 6 months, with ps every 10 years. 
   !!
   !!---------------------------------------------------------------------------
   
   CHARACTER, INTENT(IN) :: run*10
   
!$OMP MASTER
   
   dtstep = FLOAT(hourstep) * 3600.  
   
   !!
   !! Allocate dzt, dxdy and dx the first time step
   !! and reset them
   !!
   IF (iyear == yearstart .AND. &
   &    imon == monstart  .AND. &
   &    iday == daystart  .AND. &
   &   ihour == hourstart ) THEN
      
      lfirst = .true.
      
      !!
      !! Private arrays that will be kept throughout the integrations
      !!
      ALLOCATE( aa(0:KM), bb(0:KM), zxy(IMT,0:JMT), steps(12,31,4), &
      &         steps2(12,31,4) ) 
      
      dzt(:,:,:,:) = 0.
      dxdy(:,:) = 0.
      dx(:,:) = 0.
      dy(:,:) = 0.
      vlat2(:) = 0.
      aa(:) = 0.
      bb(:) = 0.
      
      istep = 0
      istep2 = 0
      DO ji=1,12
         DO jj=1,idmax(ji,1850)
            DO jk=1,4
               IF (jk == 1 .AND. jj == 1 .AND. &
                  (ji == 1 .OR. ji == 7) ) THEN
                  istep = 0
               END IF
               istep = istep + 1
               steps(ji,jj,jk) = istep
               istep2 = istep2 + 1
               steps2(ji,jj,jk) = istep2
            END DO
         END DO
      END DO
      
   ELSE
      
      lfirst = .false.
      
   END IF
   
   istep  = steps (imon,iday,INT(ihour/6)+1)
   SELECT CASE (TRIM(run))
      CASE ('ESM1-HISTR')
         istep2 = steps2(imon,iday,INT(ihour/6)+1) + 1460 * MOD(iyear,10)
      CASE ('ESM1-RCP85')
         IF (MOD(iyear,10) >= 6) THEN
            istep2 = steps2(imon,iday,INT(ihour/6)+1) + 1460 * (MOD(iyear,10)-6)
         ELSE IF (MOD(iyear,10) <= 5) THEN
            istep2 = steps2(imon,iday,INT(ihour/6)+1) + 1460 * (MOD(iyear,10)+4)
         END IF
   END SELECT
   
   PRINT*,iyear,imon,iday,ihour,istep,istep2
   
   !!
   !! Temporary arrays
   !!
   
   allocate ( txyz(IMT,0:JMT,KM), &
   &          qxyz(IMT,0:JMT,KM), &
   &          uxyz(IMT,0:JMT,KM), &
   &          vxyz(IMT,0:JMT,KM), &
   &           pxy(IMT,0:JMT),    &
   &            pt(imt,jmt),      &
   &            pu(imt,jmt),      &
   &            pv(imt,0:jmt) )
   
   geo_i(:,:,:) = 0.
   txyz(:,:,:) = 0.
   qxyz(:,:,:) = 0.
   uxyz(:,:,:) = 0.
   vxyz(:,:,:) = 0.
   pxy(:,:) = 0.
   
   
   !!
   !! Store previous time step and reset current
   !!
   dzt(:,:,:,1) = dzt(:,:,:,2)
   dzt(:,:,:,2) = 0.
   
   IF (lverbose) THEN
      PRINT*,' Day, month, step ',iday, imon, istep
   END IF
      
   IF (istep == 1 .OR. lfirst) THEN
      
      !!
      !! File names
      !!
      SELECT CASE (TRIM(run))
         CASE ('ESM1-HISTR')
            IF (imon <= 6) THEN
               string='0000010100-0000063018.nc'
            ELSE IF (imon >= 7) THEN
               string='0000070100-0000123118.nc'
            END IF
            WRITE (string(1:4),   '(i4.4)') iyear
            WRITE (string(12:15), '(i4.4)') iyear
            dataprefix='6hrLev_NorESM1-M_historical_r1i1p1_'//TRIM(string)
            string='/0000/'
            WRITE (string(2:5),'(i4.4)') iyear
            uFile=TRIM(inDataDir)//TRIM(string)//'ua_'//TRIM(dataprefix)
            vFile=TRIM(inDataDir)//TRIM(string)//'va_'//TRIM(dataprefix)
            tFile=TRIM(inDataDir)//TRIM(string)//'ta_'//TRIM(dataprefix)
            qFile=TRIM(inDataDir)//TRIM(string)//'hus_'//TRIM(dataprefix)
            zFile=TRIM(topoDir)//'orog_fx_NorESM1-M_historical_r0i0p0.nc'
            topoFile=TRIM(topoDir)//'areacella_fx_NorESM1-M_historical_r0i0p0.nc'
            
         
         CASE ('ESM1-RCP85')
            IF (imon <= 6) THEN
               string='0000010100-0000063018.nc'
            ELSE IF (imon >= 7) THEN
               string='0000070100-0000123118.nc'
            END IF
            WRITE (string(1:4),   '(i4.4)') iyear
            WRITE (string(12:15), '(i4.4)') iyear
            dataprefix='6hrLev_NorESM1-M_rcp85_r1i1p1_'//TRIM(string)
            string='0000'
            WRITE (string(1:4),'(i4.4)') iyear
            uFile=TRIM(inDataDir)//TRIM(string)//'/ua_'//TRIM(dataprefix)
            vFile=TRIM(inDataDir)//TRIM(string)//'/va_'//TRIM(dataprefix)
            tFile=TRIM(inDataDir)//TRIM(string)//'/ta_'//TRIM(dataprefix)
            qFile=TRIM(inDataDir)//TRIM(string)//'/hus_'//TRIM(dataprefix)
            zFile=TRIM(topoDir)//'orog_fx_NorESM1-M_rcp85_r0i0p0.nc'
            topoFile=TRIM(topoDir)//'areacella_fx_NorESM1-M_rcp85_r0i0p0.nc'
            
         
      END SELECT
        
      
      !!
      !! Re-grid data
      !! 
      if (tweak_zmean >= 1) then
         
         tmpFile = trim(tmpDataDir)//trim(prefix)//'tmp.nc'
         
         uFile2 = trim(tmpDataDir)//trim(prefix)//'ua.nc'
         vFile2 = trim(tmpDataDir)//trim(prefix)//'va.nc'
         tFile2 = trim(tmpDataDir)//trim(prefix)//'ta.nc'
         qFile2 = trim(tmpDataDir)//trim(prefix)//'hus.nc'
         
         call cdo_interp(uFile,uFile2,topoFile,tmpFile,rstring,'ua',1,1,1)
         call cdo_interp(vFile,vFile2,topoFile,tmpFile,rstring,'va,ps',1,1,1)
         call cdo_interp(tFile,tFile2,topoFile,tmpFile,rstring,'ta',1,0,0)
         call cdo_interp(qFile,qFile2,topoFile,tmpFile,rstring,'hus',1,0,0)
         
         uFile = uFile2
         vFile = vFile2
         tFile = tFile2
         qFile = qFile2
         
      end if
      
      !!
      !! Open the file and read the identifiers for the variables
      !!
      IF (lverbose) THEN
         PRINT*,uFile
      END IF
      ierr = NF90_OPEN( uFile, NF90_SHARE, id_ncu)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'ua', id_u) 
      CALL err(ierr)
      
      IF (lverbose) THEN
         PRINT*,vFile
      END IF
      ierr = NF90_OPEN( vFile, NF90_SHARE, id_ncv)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncv, 'va', id_v) 
      CALL err(ierr)
      
      ierr = NF90_INQ_VARID (id_ncv, 'ps', id_p) 
      CALL err(ierr)
      
      IF (lverbose) THEN
         PRINT*,tFile
      END IF
      ierr = NF90_OPEN( tFile, NF90_SHARE, id_nct)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_nct, 'ta', id_t) 
      CALL err(ierr)
      
      IF (lverbose) THEN
         PRINT*,qFile
      END IF
      ierr = NF90_OPEN( qFile, NF90_SHARE, id_ncq)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncq, 'hus', id_q) 
      CALL err(ierr)
      
      
      !!
      !! Check dimensions
      !!
      ierr = NF90_INQ_DIMID (id_ncu, 'lon', id_x)
      CALL err(ierr)
      ierr = NF90_INQ_DIMID (id_ncu, 'lat', id_y)
      CALL err(ierr)
      ierr = NF90_INQ_DIMID (id_ncu, 'lev', id_ml)
      CALL err(ierr)
      
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_x, LEN=IMT2)
      CALL err(ierr)
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_y, LEN=JMT2)
      CALL err(ierr)
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_ml, LEN=KM2)
      CALL err(ierr)
      
      !!
      !! Test so that IMT, JMT, KM that has been set agrees with the data
      !!
      IF (IMT2 /= IMT) THEN
         PRINT*,' Error: IMT in the data is not as you have set it! '
         PRINT*,IMT,IMT2
         STOP
      END IF
      IF (JMT2 /= JMT+1) THEN
         PRINT*,' Error: JMT in the data is not as you have set it! '
         PRINT*,JMT,JMT2
         STOP
      END IF
      IF (KM2 /= KM) THEN
         PRINT*,' Error: KM in the data is not as you have set it! '
         PRINT*,KM,KM2
         STOP
      END IF
      
   END IF
   
   IF (lfirst) THEN
      
      IF (lverbose) THEN
         PRINT*,' Reading lon, lat, lev, A(k), B(k) '
      END IF
      
      ierr = NF90_INQ_VARID (id_ncu, 'lon', id_lon)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'lat', id_lat)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'lev', id_lev)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'a_bnds', id_a)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'b_bnds', id_b)
      CALL err(ierr)
      
      ierr = NF90_GET_VAR (id_ncu, id_lon, vlon, start=[1], count=[IMT])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_lat, vlat2(0:JMT), start=[1], count=[JMT+1])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_lev, vlev, start=[1], count=[KM])
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'p0', id_p0)
      CALL err(ierr)
      
      ierr = NF90_GET_VAR (id_ncu, id_a, aa(1:KM), start=[2,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_a, aa(0:KM-1), start=[1,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_b, bb(1:KM), start=[2,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_b, bb(0:KM-1), start=[1,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_p0, p0)
      CALL err(ierr)
      
      !! 
      !! dx, dy, dxdy 
      !! dx at v-points, dy at u-points and dxdy at T-points
      !!
      call calc_dxdy()
      
      
      !!
      !! Read orography (Z file)
      !!
      if (tweak_zmean >= 1) then
         zFile2 = trim(tmpDataDir)//trim(prefix)//'z0.nc'
         call cdo_interp(zFile,zFile2,topoFile,tmpFile,rstring,'orog',0,0,0)
         zFile = zFile2
      end if
      
      ierr = NF90_OPEN( zFile, NF90_SHARE, id_ncz)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncz, 'orog', id_z) 
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncz, id_z, zxy(:,0:JMT), start=[  1,     1], &
      &                                                count=[IMT, JMT+1])
      CALL err(ierr)
      
      ierr = NF90_CLOSE(id_ncz)
      CALL err(ierr)
      
      if (tweak_zmean >= 1) then
         call system('rm '//trim(zFile))
      end if
      
      
   END IF
   
   ! u [m/s]
   ierr = NF90_GET_VAR (id_ncu, id_u, uxyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep],     &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   
   ! v [m/s]
   ierr = NF90_GET_VAR (id_ncv, id_v, vxyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep], &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   
   ! t [K]
   ierr = NF90_GET_VAR (id_nct, id_t, txyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep], &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   
   ! ps [Pa]
   ierr = NF90_GET_VAR (id_ncv, id_p, pxy(:,0:JMT),  &
   &                                        start=[  1,     1,  istep], &
   &                                        count=[IMT, JMT+1,      1])
   CALL err(ierr)
   
   ! q [kg/kg]
   ierr = NF90_GET_VAR (id_ncq, id_q, qxyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep], &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   
         
   IF ( (ihour == 18 .AND. iday == idmax(imon,iyear) .AND. &
   &    (imon == 6 .OR. imon == 12))  .OR. llast) THEN
      
      ierr = NF90_CLOSE(id_ncu)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncv)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_nct)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncq)
      CALL err(ierr)
      
      if (tweak_zmean >= 1) then
         call system('rm '//trim(uFile))
         call system('rm '//trim(vFile))
         call system('rm '//trim(qFile))
         call system('rm '//trim(tFile))
      end if
      
   END IF
   
   
   
   !!
   !! A-grid -> C-grid
   !!
   call a2cgrid_t(pxy(:,:),pt(:,:))
   call a2cgrid_t(zxy(:,:),geo_i(:,:,km))
   call a2cgrid_u(pxy(:,:),pu(:,:))
   call a2cgrid_v(pxy(:,:),pv(:,:))
   
   rho_i(:,:,0) = aa(km) * p0 + bb(km) * pt(:,:)
   
   DO jk=1,KM
      
      il = km - jk + 1
      ku = km - jk + 1
      kb = km - jk
      
      rho_i(:,:,jk) = aa(kb) * p0 + bb(kb) * pt(:,:)
      
      da = aa(kb) - aa(ku)
      db = bb(kb) - bb(ku)
      
      dzt(:,:,jk,2) = da * p0 + db * pt(:,:)
      vol(:,:,jk) = dzt(:,:,jk,2) * dxdy(:,:) / dg
      
      call a2cgrid_t(txyz(:,:,il),tem(:,:,jk))
      call a2cgrid_t(qxyz(:,:,il),sal(:,:,jk))
      
      !!
      !! u,v -> uflux,vflux
      !!
      call a2cgrid_u(uxyz(:,:,il),uxyz(:,:,il))
      uflux(:,:,jk) = uxyz(:,:,il) * dy(:,:) * (da * p0 + db * pu(:,:)) / dg
      call a2cgrid_v(vxyz(:,:,il),vxyz(:,:,il))
      vflux(:,:,jk) = vxyz(:,:,il) * dx(:,:) * (da * p0 + db * pv(:,:)) / dg
      
   END DO
   
   DEALLOCATE( txyz, qxyz, uxyz, vxyz, pxy, pt, pu, pv )
   
   
   !!
   !! Calculate geopotential
   !!
   call int_geo()
   
   geo(:,:,1:km) = 0.5 * (geo_i(:,:,0:km-1) + geo_i(:,:,1:km))
   rho(:,:,1:km) = 0.5 * (rho_i(:,:,0:km-1) + rho_i(:,:,1:km))
   
   !!
   !! Calculate vertical mass flux
   !!
   call calc_wflux()
   
!$OMP END MASTER
   

   RETURN
   
   END SUBROUTINE
   
   
   subroutine get_data_csiro(run)
   !!---------------------------------------------------------------------------
   !!
   !! Subroutine to get data from CSIRO-Mk3.6.0.
   !! Data is stored every 6h with one file every year (1 Jan 06 -> 1 Jan 00). 
   !!
   !!---------------------------------------------------------------------------
   
   character, intent(in) :: run*5
   
!$OMP MASTER
   
   dtstep = FLOAT(hourstep) * 3600.  
   
   !!
   !! Allocate dzt, dxdy and dx the first time step
   !! and reset them
   !!
   IF (iyear == yearstart .AND. &
   &    imon == monstart  .AND. &
   &    iday == daystart  .AND. &
   &   ihour == hourstart ) THEN
      
      lfirst = .true.
      
      !!
      !! Private arrays that will be kept throughout the integrations
      !!
      ALLOCATE( aa(0:KM), bb(0:KM), steps(12,31,4), zxy(imt,0:jmt) ) 
      
      dzt(:,:,:,:) = 0.
      dxdy(:,:) = 0.
      dx(:,:) = 0.
      dy(:,:) = 0.
      vlat2(:) = 0.
      aa(:) = 0.
      bb(:) = 0.
      
      istep = 1459
      DO ji=1,12
         DO jj=1,idmax(ji,1850)
            DO jk=1,4
               IF (ji == 1 .AND. jj == 1 .AND. jk == 2) THEN
                  istep = 0
               END IF
               istep = istep + 1
               steps(ji,jj,jk) = istep
            END DO
         END DO
      END DO
      
   ELSE
      
      lfirst = .false.
      
   END IF
   
   istep = steps(imon,iday,INT(ihour/6)+1)
   
   !!
   !! Temporary arrays
   !!
   
   allocate( txyz(imt,0:jmt,km), &
   &         qxyz(imt,0:jmt,km), &
   &         uxyz(imt,0:jmt,km), &
   &         vxyz(imt,0:jmt,km), &
   &         pxy(imt,0:jmt),      &
   &         pt(imt,1:jmt),      &
   &         pu(imt,1:jmt),      &
   &         pv(imt,0:jmt) )
   
   geo_i(:,:,:) = 0.
   txyz(:,:,:) = 0.
   qxyz(:,:,:) = 0.
   uxyz(:,:,:) = 0.
   vxyz(:,:,:) = 0.
   pxy(:,:) = 0.
   
   !!
   !! Store previous time step and reset current
   !!
   dzt(:,:,:,1) = dzt(:,:,:,2)
   dzt(:,:,:,2) = 0.
   
   !IF (lverbose) THEN
      PRINT*,' Day, month, step ',iday, imon, istep
   !END IF
      
   IF (istep == 1 .OR. lfirst) THEN
      
      !!
      !! File names
      !! One file per year
      !!
      SELECT CASE (TRIM(run))
          
         CASE ('HISTR')
            
            string='000001010600-000001010000.nc'
            if (imon == 1 .AND. iday == 1 .AND. ihour == 0) THEN
               write (string(1:4),   '(i4)') iyear-1
               write (string(14:17), '(i4)') iyear
            else
               write (string(1:4),   '(i4)') iyear
               write (string(14:17), '(i4)') iyear+1
            end if
            dataprefix='6hrLev_CSIRO-Mk3-6-0_historical_r1i1p1_'//TRIM(string)
            string='/0000/'
            if (imon == 1 .AND. iday == 1 .AND. ihour == 0) THEN
               write (string(2:5),       '(i4)') iyear-1
            else
               write (string(2:5),       '(i4)') iyear
            end if
            uFile = trim(inDataDir)//trim(string)//'ua_'//trim(dataprefix)
            vFile = trim(inDataDir)//trim(string)//'va_'//trim(dataprefix)
            tFile = trim(inDataDir)//trim(string)//'ta_'//trim(dataprefix)
            qFile = trim(inDataDir)//trim(string)//'hus_'//trim(dataprefix)
            zFile = trim(topoDir)//'orog_fx_CSIRO-Mk3-6-0_historical_r0i0p0.nc'
            topoFile = trim(topoDir)//'areacella_fx_CSIRO-Mk3-6-0_historical_r0i0p0.nc'
         
      END SELECT
            
      !!
      !! Interpolate data
      !!
      if (tweak_zmean >= 1) then
         
         tmpFile = trim(tmpDataDir)//trim(prefix)//'tmp.nc'
         
         uFile2 = trim(tmpDataDir)//trim(prefix)//'ua.nc'
         vFile2 = trim(tmpDataDir)//trim(prefix)//'va.nc'
         tFile2 = trim(tmpDataDir)//trim(prefix)//'ta.nc'
         qFile2 = trim(tmpDataDir)//trim(prefix)//'hus.nc'
         zFile2 = trim(tmpDataDir)//trim(prefix)//'z0.nc'
         
         call cdo_interp(uFile,uFile2,topoFile,tmpFile,rstring,'ua,ps',1,1,1)
         call cdo_interp(vFile,vFile2,topoFile,tmpFile,rstring,'va',1,1,1)
         call cdo_interp(tFile,tFile2,topoFile,tmpFile,rstring,'ta',1,0,0)
         call cdo_interp(qFile,qFile2,topoFile,tmpFile,rstring,'hus',1,0,0)
         
         uFile = uFile2
         vFile = vFile2
         tFile = tFile2
         qFile = qFile2
         
      end if
      
      
      !!
      !! Open the files and read the identifiers for the variables
      !!
      
      !! U file
      IF (lverbose) THEN
         PRINT*,uFile
      END IF
      ierr = NF90_OPEN( uFile, NF90_SHARE, id_ncu)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'ua', id_u) 
      CALL err(ierr)
      
      !! ps from U file
      ierr = NF90_INQ_VARID (id_ncu, 'ps', id_p) 
      CALL err(ierr)
      
      !! V file
      IF (lverbose) THEN
         PRINT*,vFile
      END IF
      ierr = NF90_OPEN( vFile, NF90_SHARE, id_ncv)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncv, 'va', id_v) 
      CALL err(ierr)
      
      !! T file
      IF (lverbose) THEN
         PRINT*,tFile
      END IF
      ierr = NF90_OPEN( tFile, NF90_SHARE, id_nct)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_nct, 'ta', id_t) 
      CALL err(ierr)
      
      !! Q file
      IF (lverbose) THEN
         PRINT*,qFile
      END IF
      ierr = NF90_OPEN( qFile, NF90_SHARE, id_ncq)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncq, 'hus', id_q) 
      CALL err(ierr)
      
      
      !!
      !! Check dimensions of uFile
      !!
      ierr = NF90_INQ_DIMID (id_ncu, 'lon', id_x)
      CALL err(ierr)
      ierr = NF90_INQ_DIMID (id_ncu, 'lat', id_y)
      CALL err(ierr)
      ierr = NF90_INQ_DIMID (id_ncu, 'lev', id_ml)
      CALL err(ierr)
      
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_x, LEN=IMT2)
      CALL err(ierr)
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_y, LEN=JMT2)
      CALL err(ierr)
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_ml, LEN=KM2)
      CALL err(ierr)
      
      !! Test so that IMT, JMT, KM that has been set agrees with the data
      IF (IMT2 /= IMT) THEN
         PRINT*,' Error: IMT in the data is not as you have set it! '
         PRINT*,IMT,IMT2
         STOP
      END IF
      IF (JMT2 /= JMT+1) THEN
         PRINT*,' Error: JMT in the data is not as you have set it! '
         PRINT*,JMT,JMT2
         STOP
      END IF
      IF (KM2 /= KM) THEN
         PRINT*,' Error: KM in the data is not as you have set it! '
         PRINT*,KM,KM2
         STOP
      END IF
      
   END IF
   
   IF (lfirst) THEN
      
      IF (lverbose) THEN
         PRINT*,' Reading lon, lat, lev, A(k), B(k) '
      END IF
      
      !!
      !! Read topography and grid info from separate files
      !!
      ierr = NF90_OPEN( topoFile, NF90_SHARE, id_nca)
      CALL err(ierr)
      
      ierr = NF90_INQ_VARID (id_ncu, 'lon', id_lon)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'lat', id_lat)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'lev', id_lev)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'a_bnds', id_a)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'b_bnds', id_b)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'p0', id_p0)
      CALL err(ierr)
      
      ierr = NF90_GET_VAR (id_ncu, id_lon, vlon, start=[1], count=[IMT])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_lat, vlat2(0:jmt), start=[1], count=[JMT+1])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_lev, vlev, start=[1], count=[KM])
      CALL err(ierr)
      
      ierr = NF90_GET_VAR (id_ncu, id_a, aa(1:KM), start=[2,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_a, aa(0:KM-1), start=[1,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_b, bb(1:KM), start=[2,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_b, bb(0:KM-1), start=[1,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_p0, p0)
      CALL err(ierr)
      
      aa(0:km) = aa(0:km) * p0
      
      ierr = NF90_CLOSE( id_nca)
      CALL err(ierr)
      
      !!
      !! Define dy and calculate dx
      !!
      call calc_dxdy()
      
      if (lverbose) then
         
         print*,'dx'
         print*,dx(1,:)
         print*,'dy'
         print*,dy(1,:)
         print*,'dxdy'
         print*,dxdy(1,:)
         
      end if
      
      !!
      !! Orography (Z file)
      !!
      if (tweak_zmean >= 1) then
         
         call cdo_interp(zFile,zFile2,topoFile,tmpFile,rstring,'orog',0,0,0)
         zFile = zFile2
      
      end if   
      
      ierr = NF90_OPEN( zFile, NF90_SHARE, id_ncz)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncz, 'orog', id_z) 
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncz, id_z, zxy(:,:),    start=[  1,   1], &
      &                                               count=[IMT, JMT+1])
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncz)
      CALL err(ierr)
      
      if (tweak_zmean >= 1) then
         
         call system('rm -f '//trim(zFile))
         
      end if
      
      
   END IF
   
   
   ierr = NF90_GET_VAR (id_ncu, id_u, uxyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep],     &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncv, id_v, vxyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep], &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_nct, id_t, txyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep], &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncu, id_p, pxy(:,0:JMT),  &
   &                                        start=[  1,     1, istep], &
   &                                        count=[IMT, JMT+1,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncq, id_q, qxyz(:,0:JMT,:), &
   &                                        start=[  1,     1,  1, istep], &
   &                                        count=[IMT, JMT+1, KM,     1])
   CALL err(ierr)
   
         
   IF ( (imon == 12 .AND. iday == 31 .AND. ihour == 18) .OR. llast) THEN
      
      ierr = NF90_CLOSE(id_ncu)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncv)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_nct)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncq)
      CALL err(ierr)
      
      if (tweak_zmean >= 1) then
         
         call system('rm -f '//trim(uFile))
         call system('rm -f '//trim(vFile))
         call system('rm -f '//trim(tFile))
         call system('rm -f '//trim(qFile))
      
      end if
      
   END IF
   
   !!
   !! A-grid -> C-grid
   !!
   call a2cgrid_t(pxy(:,:),pt(:,:))
   call a2cgrid_t(zxy(:,:),geo_i(:,:,km))
   call a2cgrid_u(pxy(:,:),pu(:,:))
   call a2cgrid_v(pxy(:,:),pv(:,:))
   
   rho_i(:,:,0) = aa(km) + bb(km) * pt(:,:)
   geo_i(:,:,km) = geo_i(:,:,km) * dg
   
   DO jk=1,KM
      
      il = km - jk + 1
      ku = km - jk + 1
      kb = km - jk
      
      rho_i(:,:,jk) = aa(kb) + bb(kb) * pt(:,:)
      
      da = aa(kb) - aa(ku)
      db = bb(kb) - bb(ku)
      
      dzt(:,:,jk,2) = da + db * pt(:,:)
      vol(:,:,jk) = dzt(:,:,jk,2) * dxdy(:,:) / dg
      
      call a2cgrid_t(txyz(:,:,il),tem(:,:,jk))
      call a2cgrid_t(qxyz(:,:,il),sal(:,:,jk))
      
      !!
      !! u,v -> uflux,vflux
      !!
      call a2cgrid_u(uxyz(:,:,il),uxyz(:,:,il))
      uflux(:,:,jk) = uxyz(:,:,il) * dy(:,:) * (da + db * pu(:,:)) / dg
      call a2cgrid_v(vxyz(:,:,il),vxyz(:,:,il))
      vflux(:,:,jk) = vxyz(:,:,il) * dx(:,:) * (da + db * pv(:,:)) / dg
      
   END DO
   
   deallocate( txyz, qxyz, uxyz, vxyz, pxy, pt, pu, pv )
   
   !!
   !! Calculate geopotential
   !!
   
   rho_i(:,:,0) = 100.
   
   call int_geo()
   
   geo(:,:,1:km) = 0.5 * (geo_i(:,:,1:km) + geo_i(:,:,0:km-1))
   rho(:,:,1:km) = 0.5 * (rho_i(:,:,1:km) + rho_i(:,:,0:km-1))
   
   !!
   !! Calculate vertical mass flux
   !!
   call calc_wflux()
   
!$OMP END MASTER


   return
   
   end subroutine
   
   
   SUBROUTINE get_data_giss(run)
   !!---------------------------------------------------------------------------
   !!
   !! Subroutine to get data from GISS-E2-R.
   !! Data is stored every 6h with one file covering 1 Jan 06Z -> 1 Jul 00Z and
   !! one for 1 Jul 06Z -> 1 Jan 00Z
   !!
   !!---------------------------------------------------------------------------
   
   CHARACTER, INTENT(IN) :: run*5
   
!$OMP MASTER
   
   dtstep = FLOAT(hourstep) * 3600.  
   
   !!
   !! Allocate dzt, dxdy and dx the first time step
   !! and reset them
   !!
   IF (iyear == yearstart .AND. &
   &    imon == monstart  .AND. &
   &    iday == daystart  .AND. &
   &   ihour == hourstart ) THEN
      
      lfirst = .true.
      
      !!
      !! Private arrays that will be kept throughout the integrations
      !!
      ALLOCATE( dzt(IMT,JMT,KM,2), dxdy(IMT,JMT), dx(IMT,0:JMT), dy(IMT,JMT), &
      &         vlat2(0:JMT), aa(0:KM), bb(0:KM), zxy(IMT,JMT), steps(12,31,4) ) 
      
      dzt(:,:,:,:) = 0.
      dxdy(:,:) = 0.
      dx(:,:) = 0.
      dy(:,:) = 0.
      vlat2(:) = 0.
      aa(:) = 0.
      bb(:) = 0.
      
      istep = 735
      DO ji=1,12
         DO jj=1,idmax(ji,1850)
            DO jk=1,4
               IF ( (ji == 1 .OR. ji == 7) .AND. jj == 1 .AND. jk == 2) THEN
                  istep = 0
               END IF
               istep = istep + 1
               steps(ji,jj,jk) = istep
            END DO
         END DO
      END DO
      
   ELSE
      
      lfirst = .false.
      
   END IF
   
   istep = steps(imon,iday,INT(ihour/6)+1)
   
   !!
   !! Temporary arrays
   !!
   ALLOCATE( geo_i(IMT,JMT,0:KM) )
   
   ALLOCATE( txyz(IMT,JMT,KM), &
   &         qxyz(IMT,JMT,KM), &
   &         uxyz(IMT,JMT,KM), &
   &         vxyz(IMT,JMT,KM), &
   &         pxy(IMT,JMT) )
   
   geo_i(:,:,:) = 0.
   txyz(:,:,:) = 0.
   qxyz(:,:,:) = 0.
   uxyz(:,:,:) = 0.
   vxyz(:,:,:) = 0.
   pxy(:,:) = 0.
   
   !!
   !! Store previous time step and reset current
   !!
   dzt(:,:,:,1) = dzt(:,:,:,2)
   dzt(:,:,:,2) = 0.
   
   !IF (lverbose) THEN
      PRINT*,' Day, month, step ',iday, imon, istep
   !END IF
      
   IF (istep == 1 .OR. lfirst) THEN
      
      !!
      !! File names
      !! One file per year
      !!
      SELECT CASE (TRIM(run))
         CASE ('HISTR')
            string='000001010600-000001010000.nc'
            IF (imon <= 6) THEN
               IF (iday == 1 .AND. ihour == 0) THEN
                  string='000007010600-000001010000.nc'
                  WRITE (string(1:4),   '(i4)') iyear-1
                  WRITE (string(14:17), '(i4)') iyear
               ELSE 
                  string='000001010600-000007010000.nc'
                  WRITE (string(1:4),   '(i4)') iyear
                  WRITE (string(14:17), '(i4)') iyear
               END IF
            ELSE IF (imon >= 7) THEN
               IF (iday == 1 .AND. ihour == 0) THEN
                  string='000001010600-000007010000.nc'
                  WRITE (string(1:4),   '(i4)') iyear
                  WRITE (string(14:17), '(i4)') iyear
               ELSE
                  string='000007010600-000001010000.nc'
                  WRITE (string(1:4),   '(i4)') iyear
                  WRITE (string(14:17), '(i4)') iyear+1
               END IF
            END IF
            
            dataprefix='6hrLev_GISS-E2-R_historical_r6i1p1_'//TRIM(string)
            
            uFile=TRIM(inDataDir)//'/'//TRIM(string(1:4))//'/ua_'//TRIM(dataprefix)
            vFile=TRIM(inDataDir)//'/'//TRIM(string(1:4))//'/va_'//TRIM(dataprefix)
            pFile=TRIM(inDataDir)//'/'//TRIM(string(1:4))//'/ps_'//TRIM(dataprefix)
            tFile=TRIM(inDataDir)//'/'//TRIM(string(1:4))//'/ta_'//TRIM(dataprefix)
            qFile=TRIM(inDataDir)//'/'//TRIM(string(1:4))//'/hus_'//TRIM(dataprefix)
            zFile=TRIM(topoDir)//'/orog_fx_GISS-E2-R_historical_r0i0p0.nc'
            topoFile=TRIM(topoDir)//'/areacella_fx_GISS-E2-R_historical_r0i0p0.nc'
         
      END SELECT
            
      
      !!
      !! Open the files and read the identifiers for the variables
      !!
      
      !! U file
      IF (lverbose) THEN
         PRINT*,uFile
      END IF
      ierr = NF90_OPEN( uFile, NF90_SHARE, id_ncu)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'ua', id_u) 
      CALL err(ierr)
      
      !! ps from U file
      IF (lverbose) THEN
         PRINT*,pFile
      END IF
      ierr = NF90_OPEN( pFile, NF90_SHARE, id_ncp)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncp, 'ps', id_p) 
      CALL err(ierr)
      
      !! V file
      IF (lverbose) THEN
         PRINT*,vFile
      END IF
      ierr = NF90_OPEN( vFile, NF90_SHARE, id_ncv)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncv, 'va', id_v) 
      CALL err(ierr)
      
      !! T file
      IF (lverbose) THEN
         PRINT*,tFile
      END IF
      ierr = NF90_OPEN( tFile, NF90_SHARE, id_nct)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_nct, 'ta', id_t) 
      CALL err(ierr)
      
      !! Q file
      IF (lverbose) THEN
         PRINT*,qFile
      END IF
      ierr = NF90_OPEN( qFile, NF90_SHARE, id_ncq)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncq, 'hus', id_q) 
      CALL err(ierr)
      
      
      !!
      !! Check dimensions of uFile
      !!
      ierr = NF90_INQ_DIMID (id_ncu, 'lon', id_x)
      CALL err(ierr)
      ierr = NF90_INQ_DIMID (id_ncu, 'lat', id_y)
      CALL err(ierr)
      ierr = NF90_INQ_DIMID (id_ncu, 'lev', id_ml)
      CALL err(ierr)
      
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_x, LEN=IMT2)
      CALL err(ierr)
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_y, LEN=JMT2)
      CALL err(ierr)
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_ml, LEN=KM2)
      CALL err(ierr)
      
      !! Test so that IMT, JMT, KM that has been set agrees with the data
      IF (IMT2 /= IMT) THEN
         PRINT*,' Error: IMT in the data is not as you have set it! '
         PRINT*,IMT,IMT2
         STOP
      END IF
      IF (JMT2 /= JMT) THEN
         PRINT*,' Error: JMT in the data is not as you have set it! '
         PRINT*,JMT,JMT2
         STOP
      END IF
      IF (KM2 /= KM) THEN
         PRINT*,' Error: KM in the data is not as you have set it! '
         PRINT*,KM,KM2
         STOP
      END IF
      
   END IF
   
   IF (lfirst) THEN
      
      IF (lverbose) THEN
         PRINT*,' Reading lon, lat, lev, A(k), B(k) '
      END IF
      
      !!
      !! Read topography and grid info from separate files
      !!
      ierr = NF90_OPEN( topoFile, NF90_SHARE, id_nca)
      CALL err(ierr)
      
      ierr = NF90_INQ_VARID (id_ncu, 'lon', id_lon)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'lat', id_lat)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'lev', id_lev)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'a_bnds', id_a)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'b_bnds', id_b)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'p0', id_p0)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_nca, 'areacella', id_area)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'lat_bnds', id_lat2)
      CALL err(ierr)
      
      ierr = NF90_GET_VAR (id_ncu, id_lon, vlon, start=[1], count=[IMT])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_lat, vlat, start=[1], count=[JMT])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_lev, vlev, start=[1], count=[KM])
      CALL err(ierr)
      
      ierr = NF90_GET_VAR (id_ncu, id_a, aa(1:KM), start=[2,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_a, aa(0:KM-1), start=[1,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_b, bb(1:KM), start=[2,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_b, bb(0:KM-1), start=[1,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_p0, p0)
      CALL err(ierr)
      
      ierr = NF90_GET_VAR (id_ncu, id_lat2, vlat2(0:JMT-1), start=[1,1], count=[1,JMT])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_lat2, vlat2(1:JMT), start=[2,1], count=[1,JMT])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_nca, id_area, dxdy(:,:), start=[1,1], count=[IMT,JMT])
      CALL err(ierr)
      
      ierr = NF90_CLOSE( id_nca)
      CALL err(ierr)
      
      !!
      !! Define dy and calculate dx
      !!
      DO jj=1,JMT
         dy(:,jj) = (vlat2(jj) - vlat2(jj-1)) * deg
      END DO
      DO jj=1,JMT-1
         dx(:,jj) =  0.5 * ( dxdy(:,jj+1)/dy(:,jj+1) + dxdy(:,jj)/dy(:,jj) ) 
      END DO
      dx(:,0) = 0.
      dx(:,JMT) = 0.
      
      
      
      IF (lverbose) THEN
         PRINT*,'dx(1,:)',dx(1,:)
         PRINT*,'dy(1,:)',dy(1,:)
         PRINT*,'dxdy(1,:)',dxdy(1,:)
      END IF
      
      
      !!
      !! Read orography (Z file)
      !!
      ierr = NF90_OPEN( zFile, NF90_SHARE, id_ncz)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncz, 'orog', id_z) 
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncz, id_z, zxy(:,:),    start=[  1,   1], &
      &                                               count=[IMT, JMT])
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncz)
      CALL err(ierr)
      
      
      
      
   END IF
   
   ierr = NF90_GET_VAR (id_ncu, id_u, uxyz(:,:,:), &
   &                                        start=[  1,   1,  1, istep],     &
   &                                        count=[IMT, JMT, KM,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncv, id_v, vxyz(:,:,:), &
   &                                        start=[  1,   1,  1, istep], &
   &                                        count=[IMT, JMT, KM,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_nct, id_t, txyz(:,:,:), &
   &                                        start=[  1,   1,  1, istep], &
   &                                        count=[IMT, JMT, KM,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncp, id_p, pxy(:,:),  &
   &                                        start=[  1,   1, istep], &
   &                                        count=[IMT, JMT,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncq, id_q, qxyz(:,:,:), &
   &                                        start=[  1,   1,  1, istep], &
   &                                        count=[IMT, JMT, KM,     1])
   CALL err(ierr)
   
     
   IF ( ((imon == 1 .OR. imon == 7) .AND. iday == 1 .AND. ihour == 0) .OR. llast) THEN
      
      ierr = NF90_CLOSE(id_ncu)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncv)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncp)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_nct)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncq)
      CALL err(ierr)
      
   END IF
   
   DO jk=1,KM
      
      il=KM-jk+1
      
      !!
      !! A-grid -> C-grid 
      !! u,v -> uflux,vflux
      !!
      DO jj=1,JMT
         
         jp = jj+1
         
         DO ji=1,IMT
            
            ip = ji+1
            IF (ip == IMT+1) THEN
               ip = 1
            END IF
            
            tem(ji,jj,jk)    =  txyz(ji,jj,il)
            sal(ji,jj,jk)    =  qxyz(ji,jj,il)
            
            da               =  0.5 * (aa(il-1) + aa(il))
            db               =  0.5 * (bb(il-1) + bb(il))
            rho  (ji,jj,jk)  =  da * p0 + db * pxy(ji,jj)
            
            da  = aa(il-1) - aa(il) 
            db  = bb(il-1) - bb(il) 
            dzt  (ji,jj,jk,2)=  da * p0 + db * pxy(ji,jj)
            vol  (ji,jj,jk)  =  dzt(ji,jj,jk,2) * dxdy(ji,jj) / dg 
            
            pp  = 0.5 * (pxy(ip,jj) + pxy(ji,jj))
            dp  = da * p0 + db * pp
            uflux(ji,jj,jk) =  0.5 *( uxyz(ip,jj,il) + uxyz(ji,jj,il) ) / &
            &                  dg * dy(ji,jj) * dp
            
            IF ( jj /= JMT) THEN
               pp  = 0.5 * (pxy(ji,jp) + pxy(ji,jj))
               dp  = da * p0 + db * pp
               vflux(ji,jj,jk) =  0.5 *( vxyz(ji,jp,il) + vxyz(ji,jj,il) ) / &
               &                  dg * dx(ji,jj) * dp 
            END IF
            
         END DO
         
      END DO
      
      !!
      !! No fluxes at poles
      !!
      vflux(:,0,:) = 0.
      vflux(:,JMT,:) = 0.
      
   END DO
   
   DEALLOCATE( txyz, qxyz, uxyz, vxyz )
   
   geo_i(:,:,KM) = zxy(:,:) * dg
   
   !!
   !! Calculate geopotential
   !!
   DO jj = 1,JMT
      DO jk = KM,1,-1   !Integrate from surface to top
         il = KM-jk+1
         DO ji = 1,IMT
                     
            ! Virtual temperature at layer mid-point
            tv = ( 1.0 + 0.61 * sal(ji,jj,jk) ) * tem(ji,jj,jk)
            ! Pressure at bottom interface 
            pc = aa(il-1) * p0 + bb(il-1) * pxy(ji,jj)
            ! Pressure at top interface 
            pm = aa(il) * p0 + bb(il) * pxy(ji,jj)
            ! Geopotential at top interface
            geo_i(ji,jj,jk-1) = geo_i(ji,jj,jk) + Rd * tv * LOG (pc/pm)
            
         END DO
      END DO
   END DO
   
   
   !!
   !! Calculate DSE, MSE
   !!
   DO jk=1,KM
      geo(:,:,jk) = 0.5 * (geo_i(:,:,jk) + geo_i(:,:,jk-1))
   END DO
   
   
   DEALLOCATE ( geo_i, pxy )
   
   !!
   !! Calculate vertical mass flux
   !!
   wflux(:,:,0) = 0.
   
   DO jj=1,JMT
      
      jm = jj-1
      
      DO jk=1,KM
         
         DO ji=1,IMT
            
            im = ji-1
            IF (im == 0) THEN
               im = IMT
            END IF
            
            wflux(ji,jj,jk) = wflux(ji,jj,jk-1) - &
            &                         ( uflux(ji,jj,jk) - uflux(im,jj,jk) + &
            &                           vflux(ji,jj,jk) - vflux(ji,jm,jk) + &
            &                          (dzt(ji,jj,jk,2) - dzt(ji,jj,jk,1))* &
            &                           dxdy(ji,jj) / dtstep / dg )
            
         END DO
      END DO
   END DO
!$OMP END MASTER


   RETURN
   
   END SUBROUTINE
   
   
   SUBROUTINE get_data_cnrm(run)
   !!---------------------------------------------------------------------------
   !!
   !! Subroutine to get data from ARPEGE-Climatv5.2.1, atmospheric model of CNRM-CM5.
   !! Data is stored every 6h with one file every month. 
   !!
   !!---------------------------------------------------------------------------
   
   CHARACTER, INTENT(IN) :: run*9
   
!$OMP MASTER
   
   dtstep = FLOAT(hourstep) * 3600.  
   
   !!
   !! Allocate dzt, dxdy and dx the first time step
   !! and reset them
   !!
   IF (iyear == yearstart .AND. &
   &    imon == monstart  .AND. &
   &    iday == daystart  .AND. &
   &   ihour == hourstart ) THEN
      
      lfirst = .true.
      
      !!
      !! Private arrays that will be kept throughout the integrations
      !!
      ALLOCATE( dzt(IMT,JMT,KM,2), dxdy(IMT,JMT), dx(IMT,0:JMT), dy(IMT,JMT), &
      &         vlat2(0:JMT), aa(0:KM), bb(0:KM), zxy(IMT,0:JMT) ) 
      
      dzt(:,:,:,:) = 0.
      dxdy(:,:) = 0.
      dx(:,:) = 0.
      dy(:,:) = 0.
      vlat2(:) = 0.
      aa(:) = 0.
      bb(:) = 0.
      
      istep = 123
      DO ji=1,12
         DO jj=1,idmax(ji,iyear)
            DO jk=1,4
               IF ( jj == 1 .AND. jk == 2 ) THEN
                  istep = 0
               END IF
               istep = istep + 1
               steps(ji,jj,jk) = istep 
            END DO
         END DO
      END DO
      
   ELSE
      
      lfirst = .false.
      
   END IF
   
   istep = steps(imon,iday,INT(FLOAT(ihour)/6.)+1)
   
   !!
   !! Temporary arrays
   !!
   ALLOCATE( geo_i(IMT,JMT,0:KM) )
   
   ALLOCATE( txyz(IMT,JMT,KM), &
   &         qxyz(IMT,JMT,KM), &
   &         uxyz(IMT,JMT,KM), &
   &         vxyz(IMT,JMT,KM), &
   &         pxy (IMT,JMT) )
   
   geo_i(:,:,:) = 0.
   txyz(:,:,:) = 0.
   qxyz(:,:,:) = 0.
   uxyz(:,:,:) = 0.
   vxyz(:,:,:) = 0.
   pxy(:,:) = 0.
   
   !!
   !! Store previous time step and reset current
   !!
   dzt(:,:,:,1) = dzt(:,:,:,2)
   dzt(:,:,:,2) = 0.
   
   IF (lverbose) THEN
      PRINT*,' Day, month, step ',iday, imon, istep
   END IF
      
   IF (istep == 1 .OR. lfirst) THEN
      
      !!
      !! File names
      !! One file per year
      !!
      SELECT CASE (TRIM(run))
         
         CASE ('CM5-HISTR')
            
            string='000000000600-00000000000000.nc'
            WRITE (string(1:4),   '(i4.4)') iyear
            WRITE (string(5:6),   '(i2.2)') imon
            WRITE (string(7:8),   '(i2.2)') 1
            IF (imon == 12) THEN
               WRITE (string(14:17), '(i4.4)') iyear+1
               WRITE (string(5:6),   '(i2.2)') 1
            ELSE
               WRITE (string(14:17), '(i4.4)') iyear
               WRITE (string(18:19), '(i2.2)') imon + 1
            END IF
            WRITE (string(20:21), '(i2.2)') 1
            
            IF ( iday == 1 .AND. ihour == 6 ) THEN
               WRITE (string(5:6),   '(i2.2)') imon-1
               IF ( imon == 1 ) THEN
                  WRITE (string(1:4),   '(i4.4)') iyear-1
                  WRITE (string(5:6),   '(i2.2)') 12   
               END IF
            END IF
            
            dataprefix='6hrLev_CNRM-CM5_historical_r1i1p1_'//TRIM(string)
            string='/0000/'
            WRITE (string(2:5),       '(i4)') iyear
            uFile=TRIM(inDataDir)//TRIM(string)//'ua_'//TRIM(dataprefix)
            vFile=TRIM(inDataDir)//TRIM(string)//'va_'//TRIM(dataprefix)
            pFile=TRIM(inDataDir)//TRIM(string)//'ps_'//TRIM(dataprefix)
            tFile=TRIM(inDataDir)//TRIM(string)//'ta_'//TRIM(dataprefix)
            qFile=TRIM(inDataDir)//TRIM(string)//'hus_'//TRIM(dataprefix)
            zFile=TRIM(topoDir)//'orog_fx_CNRM-CM5_historical_r0i0p0.nc'
            topoFile=TRIM(topoDir)//'areacella_fx_CNRM-CM5_historical_r0i0p0.nc'
            
         CASE ('CM5-RCP85')
            string='0000010100-0000123123.nc'
            WRITE (string(1:4),   '(i4)') iyear
            WRITE (string(12:15), '(i4)') iyear
            dataprefix='6hrLev_GFDL-CM3_rcp85_r1i1p1_'//TRIM(string)
            string='/0000/'
            WRITE (string(2:5),       '(i4)') iyear
            uFile=TRIM(inDataDir)//TRIM(string)//'ua_'//TRIM(dataprefix)
            vFile=TRIM(inDataDir)//TRIM(string)//'va_'//TRIM(dataprefix)
            pFile=TRIM(inDataDir)//TRIM(string)//'ps_'//TRIM(dataprefix)
            tFile=TRIM(inDataDir)//TRIM(string)//'ta_'//TRIM(dataprefix)
            qFile=TRIM(inDataDir)//TRIM(string)//'hus_'//TRIM(dataprefix)
            zFile=TRIM(topoDir)//'orog_fx_GFDL-CM3_historical_r0i0p0.nc'
            topoFile=TRIM(topoDir)//'areacella_fx_GFDL-CM3_historical_r0i0p0.nc'
         
      END SELECT
      
      
      !!
      !! Re-grid if necessary
      !!
      IF (tweak_zmean >= 1) THEN
         SELECT CASE (tweak_zmean)
            CASE (1)
               string = 'r144x72'
            CASE (2)
               string = 'r72x36'
         END SELECT
         
         tmpFile = TRIM(tmpDataDir)//TRIM(prefix)//'tmp.nc'
         tmpFile2 = TRIM(tmpDataDir)//TRIM(prefix)//'tmp2.nc'
         pFile2 = TRIM(tmpDataDir)//TRIM(prefix)//'ps.nc'
         uFile2 = TRIM(tmpDataDir)//TRIM(prefix)//'ua.nc'
         vFile2 = TRIM(tmpDataDir)//TRIM(prefix)//'va.nc'
         tFile2 = TRIM(tmpDataDir)//TRIM(prefix)//'ta.nc'
         qFile2 = TRIM(tmpDataDir)//TRIM(prefix)//'hus.nc'
         
         PRINT*,' Isolate ps '
         CALL SYSTEM('ncks -O -v ps,lon_bnds,lat_bnds '//TRIM(pFile)//' '//TRIM(tmpFile))
         PRINT*,' Include grid cell area '
         CALL SYSTEM('ncks -A -v areacella '//TRIM(topoFile)//' '//TRIM(tmpFile))
         PRINT*,' Interpolate to '//TRIM(string)//' grid'
         CALL SYSTEM('cdo -O remapbil,'//TRIM(string)//' '//TRIM(tmpFile)//' '//TRIM(pFile2)) 
         CALL SYSTEM('rm '//TRIM(tmpFile))
         
         PRINT*,' Isolate ua '
         CALL SYSTEM('ncks -O -v ua,lon_bnds,lat_bnds,lev_bnds '&
         &//TRIM(uFile)//' '//TRIM(tmpFile))
         PRINT*,' Include grid cell area '
         CALL SYSTEM('ncks -A -v areacella '//TRIM(topoFile)//' '//TRIM(tmpFile))
         PRINT*,' Interpolate to '//TRIM(string)//' grid'
         CALL SYSTEM('cdo -O remapbil,'//TRIM(string)//' '//TRIM(tmpFile)//' '//TRIM(uFile2)) 
         CALL SYSTEM('rm '//TRIM(tmpFile))
         
         PRINT*,' Isolate va '
         CALL SYSTEM('ncks -O -v va,lon_bnds,lat_bnds,lev_bnds '&
         &//TRIM(vFile)//' '//TRIM(tmpFile))
         PRINT*,' Include grid cell area '
         CALL SYSTEM('ncks -A -v areacella '//TRIM(topoFile)//' '//TRIM(tmpFile))
         PRINT*,' Interpolate to '//TRIM(string)//' grid'
         CALL SYSTEM('cdo -O remapbil,'//TRIM(string)//' '//TRIM(tmpFile)//' '//TRIM(vFile2)) 
         CALL SYSTEM('rm '//TRIM(tmpFile))
         
         PRINT*,' Isolate ta '
         CALL SYSTEM('ncks -O -v ta,lon_bnds,lat_bnds,lev_bnds '&
         &//TRIM(tFile)//' '//TRIM(tmpFile))
         PRINT*,' Include grid cell area '
         CALL SYSTEM('ncks -A -v areacella '//TRIM(topoFile)//' '//TRIM(tmpFile))
         PRINT*,' Interpolate to '//TRIM(string)//' grid'
         CALL SYSTEM('cdo -O remapbil,'//TRIM(string)//' '//TRIM(tmpFile)//' '//TRIM(tFile2)) 
         CALL SYSTEM('rm '//TRIM(tmpFile))
         
         PRINT*,' Isolate hus '
         CALL SYSTEM('ncks -O -v hus,lon_bnds,lat_bnds,lev_bnds '&
         &//TRIM(qFile)//' '//TRIM(tmpFile))
         PRINT*,' Include grid cell area '
         CALL SYSTEM('ncks -A -v areacella '//TRIM(topoFile)//' '//TRIM(tmpFile))
         PRINT*,' Interpolate to '//TRIM(string)//' grid'
         CALL SYSTEM('cdo -O remapbil,'//TRIM(string)//' '//TRIM(tmpFile)//' '//TRIM(qFile2)) 
         CALL SYSTEM('rm '//TRIM(tmpFile))
         
         pFile = pFile2
         uFile = uFile2
         vFile = vFile2
         tFile = tFile2
         qFile = qFile2
      
      END IF
      
      !!
      !! Open the files and read the identifiers for the variables
      !!
      
      !! U file
      IF (lverbose) THEN
         PRINT*,uFile
      END IF
      ierr = NF90_OPEN( uFile, NF90_SHARE, id_ncu)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'ua', id_u) 
      CALL err(ierr)
      
      !! P file
      IF (lverbose) THEN
         PRINT*,pFile
      END IF
      ierr = NF90_OPEN( pFile, NF90_SHARE, id_ncp)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncp, 'ps', id_p) 
      CALL err(ierr)
      
      !! V file
      IF (lverbose) THEN
         PRINT*,vFile
      END IF
      ierr = NF90_OPEN( vFile, NF90_SHARE, id_ncv)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncv, 'va', id_v) 
      CALL err(ierr)
      
      !! T file
      IF (lverbose) THEN
         PRINT*,tFile
      END IF
      ierr = NF90_OPEN( tFile, NF90_SHARE, id_nct)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_nct, 'ta', id_t) 
      CALL err(ierr)
      
      !! Q file
      IF (lverbose) THEN
         PRINT*,qFile
      END IF
      ierr = NF90_OPEN( qFile, NF90_SHARE, id_ncq)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncq, 'hus', id_q) 
      CALL err(ierr)
      
      
      !!
      !! Check dimensions of uFile
      !!
      ierr = NF90_INQ_DIMID (id_ncu, 'lon', id_x)
      CALL err(ierr)
      ierr = NF90_INQ_DIMID (id_ncu, 'lat', id_y)
      CALL err(ierr)
      ierr = NF90_INQ_DIMID (id_ncu, 'lev', id_ml)
      CALL err(ierr)
      
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_x, LEN=IMT2)
      CALL err(ierr)
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_y, LEN=JMT2)
      CALL err(ierr)
      ierr = NF90_INQUIRE_DIMENSION (id_ncu, id_ml, LEN=KM2)
      CALL err(ierr)
      
      !! Test so that IMT, JMT, KM that has been set agrees with the data
      IF (IMT2 /= IMT) THEN
         PRINT*,' Error: IMT in the data is not as you have set it! '
         PRINT*,IMT,IMT2
         STOP
      END IF
      IF (JMT2 /= JMT) THEN
         PRINT*,' Error: JMT in the data is not as you have set it! '
         PRINT*,JMT,JMT2
         STOP
      END IF
      IF (KM2 /= KM) THEN
         PRINT*,' Error: KM in the data is not as you have set it! '
         PRINT*,KM,KM2
         STOP
      END IF
      
   END IF
   
   IF (lfirst) THEN
      
      IF (lverbose) THEN
         PRINT*,' Reading lon, lat, lev, A(k), B(k) '
      END IF
      
      !!
      !! Read topography and grid info from separate files
      !!
      ierr = NF90_INQ_VARID (id_ncu, 'lon', id_lon)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'lat', id_lat)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'lev', id_lev)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'a_bnds', id_a)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'b_bnds', id_b)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'p0', id_p0)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'lon_bnds', id_lon2)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncu, 'lat_bnds', id_lat2)
      CALL err(ierr)
      
      ierr = NF90_GET_VAR (id_ncu, id_lon, vlon, start=[1], count=[IMT])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_lat, vlat, start=[1], count=[JMT])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_lev, vlev, start=[1], count=[KM])
      CALL err(ierr)
      
      ierr = NF90_GET_VAR (id_ncu, id_a, aa(1:KM), start=[2,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_a, aa(0:KM-1), start=[1,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_b, bb(1:KM), start=[2,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_b, bb(0:KM-1), start=[1,1], count=[1,KM])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_p0, p0)
      CALL err(ierr)
      
      ierr = NF90_GET_VAR (id_ncu, id_lat2, vlat2(0:JMT-1), start=[1,1], count=[1,JMT])
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncu, id_lat2, vlat2(1:JMT), start=[2,1], count=[1,JMT])
      CALL err(ierr)
      
      !!
      !! Define dy and calculate dx
      !!
      DO jj=1,JMT
         dy(:,jj) = (vlat2(jj) - vlat2(jj-1)) * deg
      END DO
      DO jj=0,JMT
         DO ji=1,IMT
            im = ji-1
            IF (im == 0) THEN
               im = IMT
            END IF
            dlon = vlon(ji) - vlon(im)
            IF (dlon < 0.) THEN
               dlon = dlon + 360.
            END IF
            dx(:,jj) = ABS( dlon * deg * COS (vlat2(jj) * radian) )
         END DO
      END DO
      DO jj=1,JMT
         jm = jj-1
         dxdy(:,jj) =  0.5 * (dx(:,jj) + dx(:,jm)) * dy(:,jj)
      END DO
      
      
      IF (lverbose) THEN
         PRINT*,'dx(1,:)',dx(1,:)
         PRINT*,'dy(1,:)',dy(1,:)
         PRINT*,'dxdy(1,:)',dxdy(1,:)
      END IF
      
      !!
      !! Re-grid if necessary
      !!
      IF (tweak_zmean >= 1) THEN
         PRINT*,' Isolate orog '
         tFile2 = TRIM(tmpDataDir)//'orog.nc'
         CALL SYSTEM('ncks -O -v orog,lon_bnds,lat_bnds '&
         &//TRIM(tFile)//' '//TRIM(tmpFile))
         PRINT*,' Include grid cell area '
         CALL SYSTEM('ncks -A -v areacella '//TRIM(topoFile)//' '//TRIM(tmpFile))
         PRINT*,' Interpolate to '//TRIM(string)//' grid'
         CALL SYSTEM('cdo -O remapbil,'//TRIM(string)//' '//TRIM(tmpFile)//' '//TRIM(tFile2)) 
         CALL SYSTEM('rm '//TRIM(tmpFile))
         zFile = tFile2
      END IF
      
      !!
      !! Read orography (Z file)
      !!
      ierr = NF90_OPEN( zFile, NF90_SHARE, id_ncz)
      CALL err(ierr)
      ierr = NF90_INQ_VARID (id_ncz, 'orog', id_z) 
      CALL err(ierr)
      ierr = NF90_GET_VAR (id_ncz, id_z, zxy(:,:),    start=[  1,   1], &
      &                                               count=[IMT, JMT])
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncz)
      CALL err(ierr)
      IF (tweak_tend >=1) THEN
         CALL SYSTEM('rm '//TRIM(zFile))
      END IF
      
   END IF
   
   
   ierr = NF90_GET_VAR (id_ncu, id_u, uxyz(:,1:JMT,:), &
   &                                        start=[  1,   1,  1, istep],     &
   &                                        count=[IMT, JMT, KM,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncv, id_v, vxyz(:,1:JMT,:), &
   &                                        start=[  1,   1,  1, istep], &
   &                                        count=[IMT, JMT, KM,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_nct, id_t, txyz(:,1:JMT,:), &
   &                                        start=[  1,   1,  1, istep], &
   &                                        count=[IMT, JMT, KM,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncp, id_p, pxy(:,1:JMT),  &
   &                                        start=[  1,   1, istep], &
   &                                        count=[IMT, JMT,     1])
   CALL err(ierr)
   ierr = NF90_GET_VAR (id_ncq, id_q, qxyz(:,1:JMT,:), &
   &                                        start=[  1,   1,  1, istep], &
   &                                        count=[IMT, JMT, KM,     1])
   CALL err(ierr)
   
         
   IF ( (iday == 1 .AND. ihour == 0) .OR. llast) THEN
      
      ierr = NF90_CLOSE(id_ncu)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncv)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncp)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_nct)
      CALL err(ierr)
      ierr = NF90_CLOSE(id_ncq)
      CALL err(ierr)
      
      IF (tweak_zmean >= 1) THEN
         CALL SYSTEM('rm '//TRIM(pFile))
         CALL SYSTEM('rm '//TRIM(uFile))
         CALL SYSTEM('rm '//TRIM(vFile))
         CALL SYSTEM('rm '//TRIM(qFile))
         CALL SYSTEM('rm '//TRIM(tFile))
      END IF
      
   END IF
   
   
   DO jk=1,KM
      
      il=KM-jk+1
      
      !!
      !! A-grid -> C-grid 
      !! u,v -> uflux,vflux
      !!
      DO jj=1,JMT
         
         jp = jj+1
         
         DO ji=1,IMT
            
            ip = ji+1
            IF (ip == IMT+1) THEN
               ip = 1
            END IF
            
            tem(ji,jj,jk)    =  txyz(ji,jj,il)
            sal(ji,jj,jk)    =  qxyz(ji,jj,il)
            
            pp               =  pxy(ji,jj)
            da               =  0.5 * (aa(il-1) + aa(il))
            db               =  0.5 * (bb(il-1) + bb(il))
            rho  (ji,jj,jk)  =  da * p0 + db * pp
            
            da  = aa(il-1) - aa(il) 
            db  = bb(il-1) - bb(il) 
            dzt  (ji,jj,jk,2)=  da * p0 + db * pp
            vol  (ji,jj,jk)  =  dzt(ji,jj,jk,2) * dxdy(ji,jj) / dg 
            
            pp  = 0.5 * (pxy(ip,jj) + pxy(ji,jj))
            dp  = da * p0 + db * pp
            uflux(ji,jj,jk) =  0.5 *( uxyz(ip,jj,il) + uxyz(ji,jj,il) ) / &
            &                  dg * dy(ji,jj) * dp
            
            IF ( jj /= JMT) THEN
               pp  = 0.5 * (pxy(ji,jp) + pxy(ji,jj))
               dp  = da * p0 + db * pp
               vflux(ji,jj,jk) =  0.5 *( vxyz(ji,jp,il) + vxyz(ji,jj,il) ) / &
               &                  dg * dx(ji,jj) * dp
            END IF
            
         END DO
         
      END DO
      
      !!
      !! No fluxes at poles
      !!
      vflux(:,0,:) = 0.
      vflux(:,JMT,:) = 0.
      
   END DO
   
   DEALLOCATE( txyz, qxyz, uxyz, vxyz )
   
   geo_i(:,:,KM) = zxy(:,:) * dg
   
   !!
   !! Calculate geopotential
   !!
   DO jj = 1,JMT
      DO jk = KM,1,-1   !Integrate from surface to top
         il = KM-jk+1
         DO ji = 1,IMT
                     
            ! Virtual temperature at layer mid-point
            tv = ( 1.0 + 0.61 * sal(ji,jj,jk) ) * tem(ji,jj,jk)
            ! Pressure at bottom interface 
            pc = aa(il-1) * p0 + bb(il-1) * pxy(ji,jj)
            ! Pressure at top interface 
            pm = aa(il) * p0 + bb(il) * pxy(ji,jj)
            ! Geopotential at top interface k-1
            geo_i(ji,jj,jk-1) = geo_i(ji,jj,jk) + Rd * tv * LOG (pc/pm)
            
            
         END DO
      END DO
   END DO
   
   
   !!
   !! Calculate DSE, MSE
   !!
   DO jk=1,KM
      geo(:,:,jk) = 0.5 * (geo_i(:,:,jk) + geo_i(:,:,jk-1))
   END DO
   
   DEALLOCATE ( geo_i, pxy )
   
   PRINT*,geo(1,60,:)
   PRINT*,rho(1,60,:)
   
   !!
   !! Calculate vertical mass flux
   !!
   wflux(:,:,0) = 0.
   
   DO jj=1,JMT
      
      jm = jj-1
      
      DO jk=1,KM
         
         DO ji=1,IMT
            
            im = ji-1
            IF (im == 0) THEN
               im = IMT
            END IF
            
            wflux(ji,jj,jk) = wflux(ji,jj,jk-1) - &
            &                         ( uflux(ji,jj,jk) - uflux(im,jj,jk) + &
            &                           vflux(ji,jj,jk) - vflux(ji,jm,jk) + &
            &                          (dzt(ji,jj,jk,2) - dzt(ji,jj,jk,1))* &
            &                           dxdy(ji,jj) / dtstep / dg )
            
         END DO
      END DO
   END DO
!$OMP END MASTER


   RETURN
   
   END SUBROUTINE
   
   
END MODULE

