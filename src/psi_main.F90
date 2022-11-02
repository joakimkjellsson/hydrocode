PROGRAM PSIXX_IFS
!!------------------------------------------------------------------------------
!!
!!
!!  Subroutine to read ERA-Interim data, compute stream functions, 
!!  and store computations to netCDF files.
!!
!!
!!  Original code for the oceans (NEMO ocean model output) by Kristofer Doos
!!  Re-written for the atmosphere (ERA-Interim) by Joakim Kjellsson
!!
!!
!!  v1.0 January 2012
!!
!!  EDIT February 2013: v2.0 Added OpenMP parallelization and 
!!                           GRIB->netCDF conversion with CDO
!!  EDIT April 2013:    v2.1 Added support for MERRA reanalysis on p-levels
!!
!!  EDIT June 2013:     v2.2 Time mean and zonal mean
!!                           and EC-Earth and GFDL CM3/AMIP data
!!
!!  EDIT July 2013:     v1.0 Project renamed to Glowing Octo Shame 1.0
!!                           Re-structured code, unified different versions
!!  
!!  EDIT August 2013:   v1.1 Added pre-processing flags to control what kind 
!!                           of stream functions are outputed. 
!!  
!!  EDIT August 2013:   v1.1.1 Fixed wierd OpenMP bug that sometimes caused 
!!                             nmean in the time-mean code to be zero. 
!!                             
!!  EDIT Sept. 2013:    v1.2 Added support for CCSM4, IPSL-CM5A-LR,
!!                           CanESM2, NorESM1
!!  EDIT Dec. 2013:     v1.2 Added support for CSIRO-Mk3-6-0
!!
!!  EDIT Apr 2022:      v1.3 Add ERA-5 (read from shared dir on Mistral, DKRZ)
!!
!!------------------------------------------------------------------------------
   
   USE netcdf
   USE omp_lib
   
   USE mod_data
   USE mod_vars
   USE mod_supp
   
   USE mod_grid
   
   IMPLICIT none

!!------------------------------------------------------------------------------

   INTEGER                                       ::  MR    = 501,              &
   &                                                 MR2   = 40,               &
   &                                                 LBAS  = 3,                &
   &                                                 LBAS2 = 1,                &
   &                                                 NRST  = 5
   
   INTEGER                                       ::  ji, jj, jk, jl, jn, jr,   &
   &                                                 ip, im, jp, jm, ku, kb,   &
   &                                                 il, m1, imm,              &
   &                                                 jzr, jzr2, ir2, id_rts2,  &
   &                                                 ir, ierr, id_lon, id_lat, &
   &                                                 id_lev, id_rts, id_nc,    &
   &                                                 id_mra2, id_nc2,          &
   &                                                 id_mra, id_mrb, id_psiyzv,&
   &                                                 id_psiyzw,                &
   &                                                 id_psiyr, id_psirz,       &
   &                                                 id_psixyu, id_psixyv,     &
   &                                                 id_psixzu, id_psixzw,     & 
   &                                                 id_psixru, id_psixrw,     &
   &                                                 id_psirr, id_rtsyz,       &
   &                                                 id_volyr, id_psiyr2,      &
   &                                                 id_volrr, id_lba, id_nst, &
   &                                                 id_nsta, id_nstb,         & 
   &                                                 id_tim, id_rtsyr2,id_lba2,&
   &                                                 id_vlat, id_rtsyr,id_vlon,&
   &                                                 id_psirr2, id_year,id_mon,&
   &                                                 ints, intsend,            &
   &                                                 ints2, nsteps, i_lh,      &
   &                                                 nthreads, ithread,        &
   &                                                 nmean, nmean2,            &
   &                                                 ical, nav
   
   INTEGER*8                                     ::  ntime
   INTEGER, ALLOCATABLE, DIMENSION(:)            ::  mm,mm2
   INTEGER, DIMENSION(6)                         ::  id_dims, start, count
   INTEGER, DIMENSION(12)                        ::  mois
   INTEGER, ALLOCATABLE, DIMENSION(:,:)          ::  izov
   INTEGER, ALLOCATABLE, DIMENSION(:,:,:)        ::  nyz, nyz_av
   INTEGER, ALLOCATABLE, DIMENSION(:,:,:)        ::  nyr, nyr_av, nyr2, nyr2_av
   

   REAL*4                                        ::  rmin = 0.,                &
   &                                                 rmax = 1100.,             &
   &                                                 tmin = 173.,              &
   &                                                 tmax = 323.,              &
   &                                                 smin = 0.,                &
   &                                                 smax = 25.,               &
   &                                                 dmin = 150.,              &
   &                                                 dmax = 1350.,             &
   &                                                 mmin = 150.,              &
   &                                                 mmax = 1400.,             &
   &                                                 amin = 0.6,               &
   &                                                 amax = 1.6
   
   REAL*8                                        ::  mk
   
   REAL*4                                        ::  dr, dt, ds, dd, dm,   &
   &                                                 dr2, dt2, ds2, dd2, dm2,  &
   &                                                 lh, da, da2

   
   REAL*4, ALLOCATABLE, DIMENSION(:,:,:)         ::  psiz, psiz_av,            &
   &                                                 psiz2, psiz2_av,          &
   &                                                 psxzu, psxzu_av,          &
   &                                                 psxzw, psxzw_av,          &
   &                                                 psxyu, psxyu_av,          &
   &                                                 psxyv, psxyv_av,          &
   &                                                 tem2, sal2, rho2, geo2,   &
   &                                                 dse2, mse2, alpha2
   REAL*4, ALLOCATABLE, DIMENSION(:,:,:,:)       ::  psir, psir_av, psir2,     &
   &                                                 psir2_av, pszr, pszr_av,  &
   &                                                 psxru, psxru_av,          &
   &                                                 psxrw, psxrw_av,          &
   &                                                 tsryz, tsryz_av, tsryr,   &
   &                                                 tsryr_av, tsryr2,         &
   &                                                 tsryr2_av, volyr, volyr_av
   REAL*4, ALLOCATABLE, DIMENSION(:,:,:,:,:)     ::  volrr, volrr_av, psrr,    &
   &                                                 psrr_av, psrr2, psrr2_av
   
   
   !!
   !! Fields
   !!
   REAL*4, ALLOCATABLE, DIMENSION(:)             ::  vars
   
   REAL*4, ALLOCATABLE, DIMENSION(:,:,:)         ::  dse, mse, alpha
   
   REAL*4, ALLOCATABLE, DIMENSION(:,:,:)         ::  umean, vmean, wmean,      &
   &                                                 tmean, smean, rmean,      &
   &                                                 gmean, dmean
   REAL*4, ALLOCATABLE, DIMENSION(:,:)           ::  umean2, vmean2, wmean2,   &
   &                                                 tmean2, smean2, rmean2,   &
   &                                                 gmean2, dmean2
   
   !!
   !! Grid
   !!
   REAL*4, ALLOCATABLE, DIMENSION(:)             ::  vmra,vmra2
   REAL*4, ALLOCATABLE, DIMENSION(:,:)           ::  vrts, vrts2


   CHARACTER                                     ::  outDataDir*100,           &
   &                                                 ncFile*100,               &
   &                                                 project*100,              &
   &                                                 ntimestring*10,           &
   &                                                 yearstring*4

   LOGICAL                                       ::  logp, lfirst
   
!!------------------------------------------------------------------------------

   ! Namelists must be declared first
   ! Intel allows you to declare it later, but GNU requries it to be done
   ! before any execution
   NAMELIST /TIME/    yearstart,monstart,daystart,hourstart,hourstep,intsend
   NAMELIST /DIR/     inDataDir,outDataDir,tmpDataDir,topoDir
   NAMELIST /FILE/    project,prefix   
   NAMELIST /GRID/    IMT,JMT,KM,MR,MR2,NRST,LBAS
   NAMELIST /COORD/   rmin,rmax, &
   &                  tmin,tmax, &
   &                  smin,smax, &
   &                  dmin,dmax, &
   &                  mmin,mmax, &
   &                  amin,amax
   NAMELIST /MISC/    lverbose, logp, tweak_tmean, tweak_zmean, tweak_tend, tweak_freq
   
   OPEN (8, FILE='hydro_parameters.in', STATUS='OLD', DELIM='APOSTROPHE')
   READ (8, nml=TIME)
   READ (8, nml=DIR)
   READ (8, nml=FILE)
   READ (8, nml=GRID)
   READ (8, nml=COORD)
   READ (8, nml=MISC)
   CLOSE(8)
   
   DATA mois/31,28,31,30,31,30,31,31,30,31,30,31/
   
   lfirst = .true.
   
   101 FORMAT ('     ',i3,'    ',i3,'   ',i3,'    ',i3,'    ',i3,'    ',i3,'    ',i3)
   102 FORMAT (a75)
   103 FORMAT (a75)
   104 FORMAT ('     t = [',f6.1,',',f6.1,']   s = [',f6.1,',',f6.1,           &
   &            ']   r = [',f6.1,',',f6.1,']' )
   105 FORMAT ('     ',i4,'-',i2,'-',i2,' at ',i2,':00')
   106 FORMAT ('     ',i6)
   201 FORMAT ('  YYYYMMDDHH : ',i10)
   301 FORMAT (' integrating ntime = ',i12)
   302 FORMAT (' Number of threads = ',i3)
   
!!------------------------------------------------------------------------------

   PRINT*,'--------------------------------------------------------------------'
   PRINT*,'--------------------------------------------------------------------'
   PRINT*,' Starting up...                                                     '
   PRINT*,'--------------------------------------------------------------------'
   PRINT*,'--------------------------------------------------------------------'
   PRINT*,'                                                                    '
   PRINT*,' Today is                                                           '
   CALL system('date')
   PRINT*,'                                                                    '
   PRINT*,' Parameters                                                         '
   PRINT*,'    IMT    JMT    KM     MR     MR2    LBAS  NRST                   '
   PRINT 101, IMT,JMT,KM,MR,MR2,LBAS,NRST
   PRINT*,'                                                                    '
   PRINT*,' Stream functions are computed with                                 '
   PRINT 104, tmin, tmax, smin, smax, rmin, rmax
   PRINT*,'                                                                    '
   IF (logp) THEN
      PRINT*,' Using log(pressure) instead of pressure'
      print*,' '
   END IF
   IF (lverbose) THEN
      PRINT*,' I will be pestering you with extra info since you asked for it ;-)'
      print*,' '
   END IF
   
!!------------------------------------------------------------------------------
   

!!
!! Allocation of shared variables
!!
   
   !! Wind, temperature, humidity, pressure, geopotential, DSE, MSE
   
   call init_vars()
   
   allocate( dse(IMT,JMT,KM), &
      &          mse(IMT,JMT,KM), alpha(IMT,JMT,KM) )
   
   IF ( tweak_tmean /= 0 ) THEN
      !! Time averaged data
      ALLOCATE( umean(IMT,JMT,KM), vmean(IMT,0:JMT,KM), wmean(IMT,JMT,0:KM),   &
&               tmean(IMT,JMT,KM), smean(IMT,JMT,KM), rmean(IMT,JMT,KM),       &
&               gmean(IMT,JMT,KM), dmean(IMT,JMT,KM) )
   END IF
   
   IF ( tweak_zmean == -1 .OR. tweak_zmean == -2 ) THEN
      !! Zonally averaged data
      ALLOCATE( umean2(JMT,KM), vmean2(0:JMT,KM), wmean2(JMT,0:KM),            &
&               tmean2(JMT,KM), smean2(JMT,KM),   rmean2(JMT,KM),              &
&               gmean2(JMT,KM), dmean2(JMT,KM) )
   END IF
   
   IF ( tweak_tend /= 0 ) THEN
      !! Local tendencies can only be calculated if we store two time steps
      ALLOCATE( tem2(IMT,JMT,KM), sal2(IMT,JMT,KM), rho2(IMT,JMT,KM),          &
&               geo2(IMT,JMT,KM), dse2(IMT,JMT,KM), mse2(IMT,JMT,KM),          & 
&               alpha2(IMT,JMT,KM) )
   END IF
   
   !! Grid variables: lon, lat, level, dx, dx*dy, dz, etc
   ALLOCATE( vmra(MR), vmra2(MR2),             &
&            izov(IMT,JMT), vrts(MR,NRST), vrts2(MR2,NRST) )
      
   !! Stream function in x-z coordinates
   ALLOCATE( psxzu(IMT,KM,LBAS), psxzu_av(IMT,KM,LBAS),                        &
   &         psxzw(IMT,KM,LBAS), psxzw_av(IMT,KM,LBAS) )
   
   !! Stream function in x-r coordinates
   ALLOCATE( psxru(IMT,MR,LBAS,NRST), psxru_av(IMT,MR,LBAS,NRST),              &
   &         psxrw(IMT,MR,LBAS,NRST), psxrw_av(IMT,MR,LBAS,NRST) )
   
   !! Stream function in x-y coordinates
   ALLOCATE( psxyu(IMT,JMT,LBAS), psxyu_av(IMT,JMT,LBAS),                      &
   &         psxyv(IMT,JMT,LBAS), psxyv_av(IMT,JMT,LBAS) )
   
   !! Stream function in y-z coordinates
   ALLOCATE( psiz (JMT,KM,LBAS), psiz_av (JMT,KM,LBAS),                        &
&            psiz2(JMT,KM,LBAS), psiz2_av(JMT,KM,LBAS) )
   
   !! Stream function in y-r coordinates
   ALLOCATE( psir (JMT,MR,LBAS,NRST), psir_av (JMT,MR,LBAS,NRST),              &
&            psir2(JMT,MR,LBAS2,NRST), psir2_av(JMT,MR,LBAS2,NRST) ) 
   
   !! Stream function in z-r coordinates
   ALLOCATE( pszr(MR,KM,LBAS,NRST), pszr_av(MR,KM,LBAS,NRST) )
   
   !! Stream function in r-r coordinates
   ALLOCATE( psrr (MR,MR,LBAS,NRST,NRST), psrr_av (MR,MR,LBAS,NRST,NRST),      &
&            psrr2(MR,MR,LBAS,NRST,NRST), psrr2_av(MR,MR,LBAS,NRST,NRST) )
             
   !! Atmospheric data in y-z coordinates 
   ALLOCATE( tsryz(JMT,KM,LBAS,NRST), tsryz_av(JMT,KM,LBAS,NRST),              &
&            nyz  (JMT,KM,LBAS),      nyz_av  (JMT,KM,LBAS) )
   
   !! Atmospheric data in x-y-r coordinates
   ALLOCATE( tsryr(IMT,JMT,MR2,NRST), tsryr_av(IMT,JMT,MR2,NRST),              &
&            nyr  (IMT,JMT,MR2),      nyr_av  (IMT,JMT,MR2) )
   
   !! Atmospheric data in y-r coordinates
   ALLOCATE( tsryr2(JMT,MR,NRST,NRST), tsryr2_av(JMT,MR,NRST,NRST),            &
&            nyr2  (JMT,MR,NRST),      nyr2_av  (JMT,MR,NRST) )
   
   !! Total mass in y-r coordinates
   ALLOCATE( volyr(JMT,MR,LBAS,NRST), volyr_av(JMT,MR,LBAS,NRST) )
   
   !! Total mass in r-r coordinates
   ALLOCATE( volrr(MR,MR,LBAS,NRST,NRST), volrr_av(MR,MR,LBAS,NRST,NRST) )
   
   
   !!
   !! Initialize shared variables
   !!
   psiz(:,:,:)      =  0.   !Overturning in (y,eta) with vflux
   psiz2(:,:,:)     =  0.   !Overturning in (y,eta) with wflux
   psir(:,:,:,:)    =  0.   !Overturning in (y,r) with vflux
   pszr(:,:,:,:)    =  0.   !Overturning in (y,r) with wflux
   psxyu(:,:,:)     =  0.   !Barotropic (x,y) with uflux
   psxyv(:,:,:)     =  0.   !Barotropic (x,y) with vflux
   psxzu(:,:,:)     =  0.   !Zonal overturning (x,eta) with uflux
   psxzw(:,:,:)     =  0.   !Zonal overturning (x,eta) with wflux
   psxru(:,:,:,:)   =  0.   !Zonal overturning (x,r) with uflux
   psxrw(:,:,:,:)   =  0.   !Zonal overturning (x,r) with wflux
   psrr(:,:,:,:,:)  =  0.   !Vertically integrated flux in (r,r)
   tsryz(:,:,:,:)   =  0.   !Zonally integrated r in model coordinates
   tsryr(:,:,:,:)   =  0.
   tsryr2(:,:,:,:)  =  0.   !Zonally integrated r in r coordinates
   nyz(:,:,:)       =  0    !Number of zonal points
   nyr(:,:,:)       =  0        
   nyr2(:,:,:)      =  0
   volyr(:,:,:,:)   =  0.
   volrr(:,:,:,:,:) =  0.   !Grid box mass in (r,r)
   
   !!
   !! Earth constants
   !!
   radian = pi / 180.d0
   deg    = radius * radian ! ~ 111000 metre

   !! 
   !! Coordinates
   !!
   dt = (tmax-tmin) / FLOAT (MR-1)
   ds = (smax-smin) / FLOAT (MR-1)
   dr = (rmax-rmin) / FLOAT (MR-1)
   dd = (dmax-dmin) / FLOAT (MR-1)
   dm = (mmax-mmin) / FLOAT (MR-1)
   da = (amax-amin) / FLOAT (MR-1)
   
   dt2 = (tmax-tmin) / FLOAT (MR2-1)
   ds2 = (smax-smin) / FLOAT (MR2-1)
   dr2 = (rmax-rmin) / FLOAT (MR2-1)
   dd2 = (dmax-dmin) / FLOAT (MR2-1)
   dm2 = (mmax-mmin) / FLOAT (MR2-1)
   da2 = (amax-amin) / FLOAT (MR2-1)

   
   !!
   !!  Set up the grid 
   !!
   
   DO jm=1,MR
      vmra(jm) = FLOAT(jm)
   END DO
   
   DO jm=1,MR2
      vmra2(jm) = FLOAT(jm)
   END DO
   
   DO jzr=1,NRST
      
      SELECT CASE(jzr)
      
      CASE(1)
         DO jm=1,MR
            vrts(jm,jzr) = smin + ds * FLOAT(jm-1)         !Humidity
         END DO
         
         vrts2(:,jzr) = smin + ds2 * (vmra2(:) - 1.)
         
      CASE(2)
         DO jm=1,MR
            vrts(jm,jzr) = dmin + dd * FLOAT(jm-1)         !Dry static energy
         END DO
         
         vrts2(:,jzr) = dmin + dd2 * (vmra2(:) - 1.)
         
      CASE(3)
         DO jm=1,MR
            vrts(jm,jzr) = mmin + dm * FLOAT(jm-1)         !Moist static energy
         END DO
         
         vrts2(:,jzr) = mmin + dm2 * (vmra2(:) - 1.)
      
      CASE(4)                                              !Specific volume
         
         vrts (:,jzr) = amin + da  * (vmra (:) - 1.)
         vrts2(:,jzr) = amin + da2 * (vmra2(:) - 1.)
      
      CASE(5)
         DO jm=1,MR
            vrts(jm,jzr) = rmin + dr * FLOAT(jm-1)         !Pressure
         END DO
         
         vrts2(:,jzr) = rmin + dr2 * (vmra2(:) - 1.)
      
      CASE(6)
         DO jm=1,MR
            vrts(jm,jzr) = tmin + dt * FLOAT(jm-1)         !Temperature
         END DO
         
         vrts2(:,jzr) = tmin + dt2 * (vmra2(:) - 1.)
      
      END SELECT
   
   END DO
   
   
   
   
   !!
   !! Start counting
   !!
   iyear = yearstart
   imon  = monstart
   iday  = daystart
   ihour = hourstart
   
   PRINT*,'                                                                    '
   PRINT*,'--------------------------------------------------------------------'
   PRINT*,'                                                                    '
   SELECT CASE (TRIM(project))
      
      CASE ('CanESM2-HISTR')
         PRINT*,' Analysing CanESM historical simulation, T63L35 '
         ical = 1   !365-day year 
      CASE ('CanESM2-RCP85')
         PRINT*,' Analysing CanESM RCP 8.5 simulation, T63L35 '
         ical = 1   !365-day year
      
      CASE ('CCSM-HISTR')
         PRINT*,' Analysing CCSM4 simulation, f9L26 '
         ical = 1   !365-day year 
      CASE ('CCSM-RCP85')
         PRINT*,' Analysing CCSM4 RCP8.5 simulation, f9L26 '
         ical = 1
      
      CASE ('CNRM-HISTR')
         PRINT*,' Analysing CNRM-CM5 historical simulation, T127,L31'
         ical = 0
      CASE ('CNRM-RCP85')
         PRINT*,' Analysing CNRM-CM5 RCP8.5 simulation, T127L31 '
         ical = 0
      
      case ('CSIRO-HISTR')
         print*,' Analysing CSIRO-Mk3-6-0 historical simulation, T63,L18 '
         ical = 1 
      case ('CSIRO-RCP85')
         print*,' Analysing CSIRO-Mk3-6-0 RCP8.5 simulation, T63,L18 '
         ical = 1
      
      CASE ('ERA')
         PRINT*,' Analysing ERA-Interim reanalysis data, T255L60 '
         ical = 0
      CASE ('ERA5') 
         PRINT*,' Analysing ERA-5 reanalysis data, Tco639L137 '
         ical = 0
      
      CASE ('GFDL-AMIP')
         PRINT*,' Analysing GFDL AMIP control simulation '
         ical = 0
      CASE ('GFDL-HISTR')
         PRINT*,' Analysing GFDL historical simulation, C48L48 '
         ical = 1   !365-day year
      CASE ('GFDL-RCP85')
         PRINT*,' Analysing GFDL RCP 8.5 simulation, C48L48 '
         ical = 1   !365-day year
      
      CASE ('GISS-HISTR')
         PRINT*,' Analysing GISS-E2-R historical simulation '
         ical = 1   !365-day year
      
      CASE ('HAD-HIST')
         PRINT*,' Analysing HadGEM historical simulation, N96L38 '
         ical = 2   !360-day year (all months 30 days)
      
      CASE ('IFS-SHC')
         PRINT*,' Analysing EC-Earth historical simulation, T159L62 '
         ical = 0  !Leap year
      CASE ('IFS-SS41')
         PRINT*,' Analysing EC-Earth RCP4.5 simulation, T159L62 '
         ical = 0
      CASE ('IFS-SS81')
         PRINT*,' Analysing EC-Earth RCP8.5 simulation, T159L62 ' 
         ical = 0
      
      CASE ('IPSL-HISTR')
         PRINT*,' Analysing IPSL-CM5A historical simulation, 144x143xL39 '
         ical = 1   !365-day year 
      CASE ('IPSL-RCP85')
         PRINT*,' Analysing IPSL-CM5A RCP 8.5 simulation, 144x143xL39 '
         ical = 1   !365-day year
         
      
      CASE ('MERRA')
         PRINT*,' Analysing MERRA reanalysis data '
         ical = 0
      
      
      CASE ('NorESM1-HISTR')
         PRINT*,' Analysing NorESM1-M historical simulation, f19L26 '
         ical = 1   !365-day year       
      CASE ('NorESM1-RCP85')
         PRINT*,' Analysing NorESM1-M RCP 8.5 simulation, f19L26 '
         ical = 1   !365-day year              
      
      CASE DEFAULT
         PRINT*,' Run: ',TRIM(project)
         PRINT*,' This project is not implemented yet. Stopping. '
         STOP
   END SELECT
   
   SELECT CASE (ical)
      CASE(0)
         PRINT*,' Calendar has leap years '
      CASE(1)
         PRINT*,' Calendar has no leap years (all years 365 days)'
      CASE(2)
         PRINT*,' Calendar has 360-day years (all months 30 days)'
   END SELECT
   
   !!
   !! Length of each month
   !!   
   DO ji=1000,3000
      DO jk=1,12
         IF (ical == 2) THEN
            idmax(jk,ji) = 30   ! 360-day year
         ELSE
            idmax(jk,ji) = mois(jk)
         END IF
      END DO
      IF (ical == 0) THEN
         IF( (MOD(ji,4) == 0 .AND. MOD(ji,100) /= 0) .OR.  &
                                               &  MOD(ji,400) == 0 ) THEN
            idmax(2,ji) = 29
         END IF
      END IF
   END DO
      
   PRINT*,'                                                                    '
   PRINT*,'--------------------------------------------------------------------'
   PRINT*,'                                                                    '
   PRINT*,' Integrations start on                                              '
   PRINT 105, yearstart,monstart,daystart,hourstart
   PRINT*,' and will go on for                                                 '
   PRINT 106, intsend
   PRINT*,' time steps                                                         '
   PRINT*,' each of                                                            '
   PRINT 106, hourstep
   PRINT*,' hours                                                              '
   PRINT*,' '
#ifdef psixy
   PRINT*,' Calculating x-y diagnostics '
#endif
#ifdef psixz
   PRINT*,' Calculating x-z diagnostics '
#endif
#ifdef psixr
   PRINT*,' Calculating x-r diagnostics '
#endif
#ifdef psiyz
   PRINT*,' Calculating y-z diagnostics '
#endif
#ifdef psiyr
   PRINT*,' Calculating y-r diagnostics '
#endif
#ifdef psirr
   PRINT*,' Calculating r-r diagnostics '
#endif
   PRINT*,'--------------------------------------------------------------------'
   
   SELECT CASE (tweak_freq)
      
      CASE (0)
         PRINT*,' Storing output every month '
         PRINT*,'--------------------------------------------------------------------'
      CASE (1:)
         PRINT*,' Storing every ',tweak_freq,' time steps'
         PRINT*,'--------------------------------------------------------------------'
   
   END SELECT
   
   SELECT CASE (tweak_zmean)
      
      CASE (0)
      
         PRINT*,' Using full zonal resolution '
         PRINT*,'--------------------------------------------------------------------'      
         
      CASE (-1)
      
         PRINT*,' Using a zonal mean '
         PRINT*,'--------------------------------------------------------------------'
         
      CASE (-2)
      
         PRINT*,' Using deviations from zonal mean '
         PRINT*,'--------------------------------------------------------------------'
         
      CASE (1:)
         
         write (lonstring,*) imt
         write (latstring,*) jmt+1
         rstring = 'r'//trim(adjustl(lonstring))//'x'//trim(adjustl(latstring))
         print*,' Data will be re-gridded to a ',trim(rstring),' resolution '
         print*,'--------------------------------------------------------------------'
         
      CASE DEFAULT
      
         PRINT*,' tweak_zmean should be -2,-1,0,1,etc. '
         PRINT*,'--------------------------------------------------------------------'
         STOP
         
   END SELECT
      
   SELECT CASE (tweak_tmean)
      
      CASE (0)
         PRINT*,' Using full temporal resolution '
         PRINT*,'--------------------------------------------------------------------'      
      CASE (-1) 
         PRINT*,' Using daily averages '
         PRINT*,'--------------------------------------------------------------------'
      CASE (-2)
         PRINT*,' Using monthly averages '
         PRINT*,'--------------------------------------------------------------------'
      CASE (-3)
         PRINT*,' Using annual averages '
         PRINT*,'--------------------------------------------------------------------'
      CASE (:-4)
         PRINT*,' tweak_tmean should be positive or 0,-1,-2,-3 '
         PRINT*,'--------------------------------------------------------------------'
         STOP
      CASE DEFAULT
         PRINT*,' Averaging every ',tweak_tmean,'th time step'
         PRINT*,'--------------------------------------------------------------------'
         
   END SELECT
   
   SELECT CASE (tweak_tend)
      
      CASE(0)
         PRINT*,' Calculating only advective part of stream functions '
         imm = 6
         PRINT*,'--------------------------------------------------------------------'
      CASE DEFAULT
         PRINT*,' Storing tendencies in psrr2 '
         imm = 7
         PRINT*,'--------------------------------------------------------------------'
   
   END SELECT
   
     
   

   !!
   !! OpenMP parallel region begins here
   !! All threads can access the field variables and stream functions
   !! but loop indices and temporary variables are kept private
   !!

!$OMP PARALLEL DEFAULT(shared) &
!$OMP & PRIVATE(jl, jk, jj, ji, im, ip, jm, jp, kb, ku, il,&
!$OMP & mm, mm2, ints, nsteps, id_nc2, &
!$OMP & ierr, id_nc, id_psiyzv, id_psiyzw, ir, ir2, jzr, jzr2, jn, m1,&
!$OMP & jr, vars, ithread, nmean2, mk, nav)
    
!$OMP MASTER  
   nthreads = OMP_GET_NUM_THREADS() 
   PRINT*,' The code will the parallelized using OpenMP threading              '
   PRINT 302, nthreads
      
!$OMP END MASTER

   !!------------------------------
   !!
   !!  START OF MAIN TIME LOOP
   !!
   !!------------------------------
   
   nav = 0
   
   DO ints=1,intsend

!$OMP BARRIER            
!$OMP MASTER   
   !!
   !! Update clock
   !!
   IF (ints /= 1) THEN
      
      ihour = ihour + hourstep
      IF( ihour >= 24) THEN
      
         ihour = 0
         iday  = iday + 1
     
         IF( iday > idmax(imon,iyear)) THEN
        
            iday = iday - idmax(imon,iyear)
            imon = imon + 1
         
            IF( imon == 13) THEN
            
               imon  = 1
               iyear = iyear + 1
         
            END IF
         END IF
      END IF
      
   END IF
   
   IF ( (ihour == 0 .AND. iday == 1 .AND. imon == 1) .OR. ints == 1 ) THEN !Happy new year
      
      ints2 = 0
      
      WRITE(yearstring,'(i4)') iyear
      
      !!
      !! Create a netCDF file.
      !! Define dimensions.
      !!
      WRITE ( ntimestring(1:4), '(i4)' ) iyear
      WRITE ( ntimestring(5:6), '(i2.2)' ) imon
      WRITE ( ntimestring(7:8), '(i2.2)' ) iday
      WRITE ( ntimestring(9:10), '(i2.2)' ) ihour
      
      ncFile = TRIM(outDataDir)//TRIM(prefix)//'_'//TRIM(yearstring)//'.nc'
      
      PRINT*,' Opening new output stream: ',ncFile
      
      ierr = NF90_CREATE( ncFile, NF90_SHARE, id_nc )
      CALL err(ierr)
      
      ierr = NF90_DEF_DIM( id_nc, 'longitude', IMT, id_lon )
      CALL err(ierr)
      ierr = NF90_DEF_DIM( id_nc, 'latitude', JMT, id_lat )
      CALL err(ierr)
      ierr = NF90_DEF_DIM( id_nc, 'mod_lev', KM, id_lev )
      CALL err(ierr)
      ierr = NF90_DEF_DIM( id_nc, 'lbas', LBAS, id_lba )
      CALL err(ierr)
      ierr = NF90_DEF_DIM( id_nc, 'lbas2', LBAS2, id_lba2 )
      CALL err(ierr)
      ierr = NF90_DEF_DIM( id_nc, 'coords', NRST, id_nst )
      CALL err(ierr)
      ierr = NF90_DEF_DIM( id_nc, 'coord1', NRST, id_nsta )
      CALL err(ierr)
      ierr = NF90_DEF_DIM( id_nc, 'coord2', NRST, id_nstb )
      CALL err(ierr)
      ierr = NF90_DEF_DIM( id_nc, 'mra_lev', MR, id_mra )
      CALL err(ierr)
      ierr = NF90_DEF_DIM( id_nc, 'mrb_lev', MR, id_mrb )
      CALL err(ierr)
      ierr = NF90_DEF_DIM( id_nc, 'mra_lev2', MR2, id_mra2 )
      CALL err(ierr)
      ierr = NF90_DEF_DIM( id_nc, 'time', NF90_UNLIMITED, id_tim )
      CALL err(ierr)
      
      
      !!
      !! Define variables
      !!
      ierr = NF90_DEF_VAR(id_nc, 'longitude', NF90_FLOAT, id_lon, id_vlon)
      CALL err(ierr)
      
      ierr = NF90_DEF_VAR(id_nc, 'latitude', NF90_FLOAT, id_lat, id_vlat)
      CALL err(ierr)
      
      ierr = NF90_DEF_VAR(id_nc, 'year', NF90_INT, id_tim, id_year)
      CALL err(ierr)
      
      ierr = NF90_DEF_VAR(id_nc, 'month', NF90_INT, id_tim, id_mon)
      CALL err(ierr)
      
      id_dims(1:2) = [id_mra, id_nst]
      ierr = NF90_DEF_VAR(id_nc, 'rts_lev', NF90_FLOAT, id_dims(1:2), id_rts)
      CALL err(ierr)
      
      id_dims(1:2) = [id_mra2, id_nst]
      ierr = NF90_DEF_VAR(id_nc, 'rts_lev2', NF90_FLOAT, id_dims(1:2), id_rts2)
      CALL err(ierr)
#ifdef psiyz      
      id_dims(1:4) = [id_lat, id_lev, id_lba, id_tim]
      ierr = NF90_DEF_VAR(id_nc, 'psiyzv', NF90_FLOAT, id_dims(1:4), id_psiyzv)
      CALL err(ierr)
      
      id_dims(1:4) = [id_lat, id_lev, id_lba, id_tim]
      ierr = NF90_DEF_VAR(id_nc, 'psiyzw', NF90_FLOAT, id_dims(1:4), id_psiyzw)
      CALL err(ierr)
      
      id_dims(1:5) = [id_lat, id_lev, id_lba, id_nst, id_tim]
      ierr = NF90_DEF_VAR(id_nc, 'rtsyz', NF90_FLOAT, id_dims(1:5), id_rtsyz)
      CALL err(ierr)
#endif
#ifdef psizr      
      id_dims(1:5) = [id_mra, id_lev, id_lba, id_nst, id_tim]
      ierr = NF90_DEF_VAR(id_nc, 'psirz', NF90_FLOAT, id_dims(1:5), id_psirz)
      CALL err(ierr)
#endif
#ifdef psiyr      
      id_dims(1:5) = [id_lat, id_mra, id_lba, id_nst, id_tim]
      ierr = NF90_DEF_VAR(id_nc, 'psiyr', NF90_FLOAT, id_dims(1:5), id_psiyr)
      CALL err(ierr)
      
      id_dims(1:5) = [id_lat, id_mra, id_lba2, id_nst, id_tim]
      ierr = NF90_DEF_VAR(id_nc, 'psiyr2', NF90_FLOAT, id_dims(1:5), id_psiyr2)
      CALL err(ierr)
      
      id_dims(1:5) = [id_lon, id_lat, id_mra2, id_nst, id_tim]
      ierr = NF90_DEF_VAR(id_nc, 'rtsyr', NF90_FLOAT, id_dims(1:5), id_rtsyr)
      CALL err(ierr)
      
      id_dims(1:5) = [id_lat, id_mra, id_nst, id_nst, id_tim]
      ierr = NF90_DEF_VAR(id_nc, 'rtsyr_zm', NF90_FLOAT, id_dims(1:5), id_rtsyr2)
      CALL err(ierr)
      
      id_dims(1:5) = [id_lat, id_mra, id_lba, id_nst, id_tim]
      ierr = NF90_DEF_VAR(id_nc, 'volyr', NF90_FLOAT, id_dims(1:5), id_volyr)
      CALL err(ierr)
#endif
#ifdef psirr      
      id_dims(1:6) = [id_mra, id_mrb, id_lba, id_nsta, id_nstb, id_tim]
      ierr = NF90_DEF_VAR(id_nc, 'psirr', NF90_FLOAT, id_dims(1:6), id_psirr)
      CALL err(ierr)
      
      id_dims(1:6) = [id_mra, id_mrb, id_lba, id_nsta, id_nstb, id_tim]
      ierr = NF90_DEF_VAR(id_nc, 'psirr2', NF90_FLOAT, id_dims(1:6), id_psirr2)
      CALL err(ierr)
      
      id_dims(1:6) = [id_mra, id_mrb, id_lba, id_nsta, id_nstb, id_tim]
      ierr = NF90_DEF_VAR(id_nc, 'volrr', NF90_FLOAT, id_dims(1:6), id_volrr)
      CALL err(ierr)
#endif
#ifdef psixz      
      id_dims(1:4) = [id_lon, id_lev, id_lba, id_tim]
      ierr = NF90_DEF_VAR(id_nc, 'psixzu', NF90_FLOAT, id_dims(1:4), id_psixzu)
      CALL err(ierr)
      
      id_dims(1:4) = [id_lon, id_lev, id_lba, id_tim]
      ierr = NF90_DEF_VAR(id_nc, 'psixzw', NF90_FLOAT, id_dims(1:4), id_psixzw)
      CALL err(ierr)
#endif
#ifdef psixy      
      id_dims(1:4) = [id_lon, id_lat, id_lba, id_tim]
      ierr = NF90_DEF_VAR(id_nc, 'psixyu', NF90_FLOAT, id_dims(1:4), id_psixyu)
      CALL err(ierr)
      
      id_dims(1:4) = [id_lon, id_lat, id_lba, id_tim]
      ierr = NF90_DEF_VAR(id_nc, 'psixyv', NF90_FLOAT, id_dims(1:4), id_psixyv)
      CALL err(ierr)
#endif
#ifdef psixr      
      id_dims(1:5) = [id_lon, id_mra, id_lba, id_nst, id_tim]
      ierr = NF90_DEF_VAR(id_nc, 'psixru', NF90_FLOAT, id_dims(1:5), id_psixru)
      CALL err(ierr)
      
      id_dims(1:5) = [id_lon, id_mra, id_lba, id_nst, id_tim]
      ierr = NF90_DEF_VAR(id_nc, 'psixrw', NF90_FLOAT, id_dims(1:5), id_psixrw)
      CALL err(ierr)
#endif      
      
   
      !!
      !! Set attributes
      !!
#ifdef psiyz
      ierr = NF90_PUT_ATT(id_nc, id_psiyzv, 'name',  &
      & 'overturning stream function in model coordinates' )
      CALL err(ierr)
      ierr = NF90_PUT_ATT(id_nc, id_psiyzv, 'units',  &
      & 'kilograms per second' )
      CALL err(ierr)
      ierr = NF90_PUT_ATT(id_nc, id_psiyzv, '_FillValue', &
      & undef)
      CALL err(ierr)
      
      ierr = NF90_PUT_ATT(id_nc, id_psiyzw, 'name',  &
      & 'overturning stream function in model coordinates' )
      CALL err(ierr)
      ierr = NF90_PUT_ATT(id_nc, id_psiyzw, 'units',  &
      & 'kilograms per second' )
      CALL err(ierr)
      ierr = NF90_PUT_ATT(id_nc, id_psiyzw, '_FillValue', &
      & undef)
      CALL err(ierr)
      
      ierr = NF90_PUT_ATT(id_nc, id_rtsyz, 'name',  &
      & 'zonally averaged T,q,p,s,h in model coordinates' )
      CALL err(ierr)
      ierr = NF90_PUT_ATT(id_nc, id_rtsyz, '_FillValue', &
      & undef)
      CALL err(ierr)
#endif
#ifdef psiyr      
      ierr = NF90_PUT_ATT(id_nc, id_psiyr, 'name',  &
      & 'overturning stream function in various coordinates' )
      CALL err(ierr)
      ierr = NF90_PUT_ATT(id_nc, id_psiyr, 'units',  &
      & 'kilograms per second' )
      CALL err(ierr)
      
      ierr = NF90_PUT_ATT(id_nc, id_rtsyr, 'name',  &
      & 'full 3D T,q,p,s,h in various coordinates' )
      CALL err(ierr)
      ierr = NF90_PUT_ATT(id_nc, id_rtsyr, '_FillValue', &
      & undef)
      CALL err(ierr)
      
      ierr = NF90_PUT_ATT(id_nc, id_rtsyr2, 'name',  &
      & 'zonal mean T,q,p,s,h in various coordinates' )
      CALL err(ierr)
      ierr = NF90_PUT_ATT(id_nc, id_rtsyr2, '_FillValue', &
      & undef)
      CALL err(ierr)
      
      ierr = NF90_PUT_ATT(id_nc, id_volyr, 'name',  &
      & 'total mass in y,r coordinates' )
      CALL err(ierr)
      ierr = NF90_PUT_ATT(id_nc, id_volyr, 'units',  &
      & 'kilograms' )
      CALL err(ierr)
#endif
#ifdef psizr      
      ierr = NF90_PUT_ATT(id_nc, id_psirz, 'name',  &
      & 'stream function in r and model coordinates' )
      CALL err(ierr)
      ierr = NF90_PUT_ATT(id_nc, id_psirz, 'units',  &
      & 'kilograms per second' )
      CALL err(ierr)
#endif
#ifdef psirr      
      ierr = NF90_PUT_ATT(id_nc, id_psirr, 'name',  &
      & 'advective part of r,r stream function' )
      CALL err(ierr)
      ierr = NF90_PUT_ATT(id_nc, id_psirr, 'units',  &
      & 'kilograms per second' )
      CALL err(ierr)
      
      ierr = NF90_PUT_ATT(id_nc, id_psirr2, 'name',  &
      & 'local time derivative part of r,r stream function ' )
      CALL err(ierr)
      ierr = NF90_PUT_ATT(id_nc, id_psirr2, 'units',  &
      & 'kilograms per second' )
      CALL err(ierr)
      
      ierr = NF90_PUT_ATT(id_nc, id_volrr, 'name',  &
      & 'total mass in r,r coordinates' )
      CALL err(ierr)
      ierr = NF90_PUT_ATT(id_nc, id_volrr, 'units',  &
      & 'kilograms' )
      CALL err(ierr)
#endif
#ifdef psixy      
      ierr = NF90_PUT_ATT(id_nc, id_psixyu, 'name', &
      & 'barotropic stream function in x,y coordinates' )
      CALL err(ierr)
      ierr = NF90_PUT_ATT(id_nc, id_psixyu, 'units', &
      & 'kilograms per second' )
      
      ierr = NF90_PUT_ATT(id_nc, id_psixyv, 'name', &
      & 'zonal stream function in x,y coordinates' )
      CALL err(ierr)
      ierr = NF90_PUT_ATT(id_nc, id_psixyv, 'units', &
      & 'kilograms per second' )
#endif
#ifdef psixr      
      ierr = NF90_PUT_ATT(id_nc, id_psixru, 'name', &
      & 'zonal stream function in x,r coordinates' )
      CALL err(ierr)
      ierr = NF90_PUT_ATT(id_nc, id_psixru, 'units', &
      & 'kilograms per second' )
      
      ierr = NF90_PUT_ATT(id_nc, id_psixrw, 'name', &
      & 'zonal stream function in x,r coordinates' )
      CALL err(ierr)
      ierr = NF90_PUT_ATT(id_nc, id_psixrw, 'units', &
      & 'kilograms per second' )
#endif
#ifdef psixz      
      ierr = NF90_PUT_ATT(id_nc, id_psixzu, 'name', &
      & 'zonal stream function in x,z coordinates' )
      CALL err(ierr)
      ierr = NF90_PUT_ATT(id_nc, id_psixzu, 'units', &
      & 'kilograms per second' )
      
      ierr = NF90_PUT_ATT(id_nc, id_psixzw, 'name', &
      & 'zonal stream function in x,r coordinates' )
      CALL err(ierr)
      ierr = NF90_PUT_ATT(id_nc, id_psixzw, 'units', &
      & 'kilograms per second' )
#endif      
      
      !!
      !! End define
      !!
      ierr = NF90_ENDDEF(id_nc)
      CAll err(ierr)
   
   END IF
   
!$OMP END MASTER
!$OMP BARRIER
   
   IF ( (iday == 1 .AND. ihour == 0. .AND. tweak_freq == 0) .OR. &
   &     ints == 1 .OR. &
   &    (nsteps == tweak_freq .AND. tweak_freq /= 0) ) THEN 
      nsteps = 0
   END IF   

!$OMP MASTER
   IF ( nsteps == 0 ) THEN 
   
   !!
   !! Reset fields
   !!
   psiz_av   = 0. 
   psiz2_av  = 0. 
   psir_av   = 0. 
   psir2_av  = 0.
   pszr_av   = 0.
   psrr_av   = 0.
   psrr2_av  = 0.
   psxyu_av  = 0.
   psxyv_av  = 0.
   psxzu_av  = 0.
   psxzw_av  = 0.
   psxru_av  = 0.
   psxrw_av  = 0.
   tsryz_av  = 0.
   tsryr_av  = 0.
   tsryr2_av = 0. 
   volrr_av  = 0.
   volyr_av  = 0. 
   nyz_av    = 0
   nyr_av    = 0
   nyr2_av   = 0
   
   END IF
   
   IF ( ints == intsend ) THEN
      llast = .true.
   ELSE
      llast = .false.
   END IF
   
   psiz  = 0.
   psiz2 = 0. 
   psir  = 0.
   psir2 = 0. 
   pszr  = 0.
   psrr  = 0.
   psrr2 = 0. 
   psxyu = 0.
   psxyv = 0.
   psxzu = 0.
   psxzw = 0.
   psxru = 0.
   psxrw = 0.
   tsryz = 0.
   tsryr = 0.
   tsryr2 = 0. 
   volrr   = 0.
   volyr   = 0. 
   nyz   = 0
   nyr   = 0
   nyr2  = 0
   
   uflux = 0.
   vflux = 0.
   wflux = 0.
   
   ntime = 1000000 * iyear + 10000 * imon + 100 * iday + ihour
   PRINT 301, ntime
   
   !!
   !! Store previous time step
   !!
   IF ( tweak_tend /= 0 .and. ints /= 1) THEN
      IF (  tweak_tmean == 0 .OR. &
         & (tweak_tmean == -1 .AND. ihour == 0 )                            .OR.&
         & (tweak_tmean == -2 .AND. ihour == 0 .AND. iday == 1 )            .OR.&
         & (tweak_tmean == -3 .AND. ihour == 0 .AND. iday == 1 .AND. imon == 1 ).OR.&
         &  nmean == 0 .and. lfirst .neqv. .true. )  THEN
         PRINT*,' Storing variables from last stored step '
         tem2 = tem
         sal2 = sal
         rho2 = rho
         geo2 = geo
         dse2 = dse
         mse2 = mse
         alpha2 = alpha
      END IF
   END IF
   
   !!
   !! Call routine for getting the data
   !!
   
   if (lverbose) then
         print*,' Reading new time step '
   end if
   
   SELECT CASE (TRIM(project))
      
      CASE ('CanESM2-HISTR')
         CALL get_data_canesm('ESM2-HISTR')
      CASE ('CanESM2-RCP85')
         CALL get_data_canesm('ESM2-RCP85')
      
      CASE ('CCSM-HISTR')
         CALL get_data_ccsm('HISTR')
      CASE ('CCSM-RCP85')
         CALL get_data_ccsm('RCP85')
      
      CASE ('CNRM-HISTR')
         CALL get_data_cnrm('HISTR')
      CASE ('CNRM-RCP85')
         CALL get_data_cnrm('RCP85')
      
      CASE ('CSIRO-HISTR')
         CALL get_data_csiro('HISTR')
      case ('CSIRO-RCP85')
         call get_data_csiro('RCP85')
      
      CASE ('ERA')
         CALL get_data_era()
      CASE ('ERA5')
         CALL get_data_era5()
      
      CASE ('GFDL-AMIP')
         CALL get_data_gfdl('AMIP-CONT')
      CASE ('GFDL-HISTR')
         CALL get_data_gfdl('CM3-HISTR')
      CASE ('GFDL-RCP85')
         CALL get_data_gfdl('CM3-RCP85')
      
      CASE ('GISS-HISTR')
         CALL get_data_giss('HISTR')
      
      CASE ('HAD-HIST')
         CALL get_data_hadgem('HAD-HISTR')
      
      CASE ('IFS-SHC')
         CALL get_data_ifs('IFS-HISTR')
      CASE ('IFS-SS41')
         CALL get_data_ifs('IFS-RCP45')
      CASE ('IFS-SS81')
         CALL get_data_ifs('IFS-RCP85')
      
      CASE ('IPSL-HISTR')
         CALL get_data_ipsl('IPSL-HISTR')
      CASE ('IPSL-RCP85')
         CALL get_data_ipsl('IPSL-RCP85')
      
      CASE ('MERRA')
         CALL get_data_merra()
      
      CASE ('NorESM1-HISTR')
         CALL get_data_noresm('ESM1-HISTR')
      CASE ('NorESM1-RCP85')
         CALL get_data_noresm('ESM1-RCP85')   
      
      CASE DEFAULT
         PRINT*,TRIM(project)
         PRINT*,' This project has not been implemented yet '
         STOP
         
   END SELECT
      
!$OMP END MASTER
!$OMP BARRIER
IF (lverbose) THEN
   PRINT*,'new field',ihour,iday,imon,iyear,OMP_GET_THREAD_NUM()
END IF

!$OMP BARRIER

   IF ( ints == 1 ) THEN

!$OMP MASTER 
      !!
      !! Define regions
      !!
      DO jj=1,JMT
      
         IF(     (vlat(jj) >= -90. .AND. vlat(jj) < -60.) .OR. &   !Polar regions
            &       (vlat(jj) >  60.  .AND. vlat(jj) <= 90. ) ) THEN
            izov(:,jj) = 3
         ELSEIF( (vlat(jj) >= -60. .AND. vlat(jj) < -30.) .OR. &   !Midlatitudes
            &       (vlat(jj) >  30.  .AND. vlat(jj) <= 60. ) ) THEN
            izov(:,jj) = 2
         ELSEIF( (vlat(jj) >= -30. .AND. vlat(jj) <= 30. ) ) THEN   !Tropics
            izov(:,jj) = 1
         END IF
    
      END DO
           
      !!
      !! Fill coordinate variables
      !!
      IF (lverbose) THEN
         PRINT*,'vlon',vlon
         PRINT*,'vlat',vlat
         PRINT*,'vlev',vlev
      END IF
      ierr = NF90_PUT_VAR(id_nc, id_vlon, vlon)
      CALL err(ierr)
      ierr = NF90_PUT_VAR(id_nc, id_vlat, vlat)
      CALL err(ierr)
      ierr = NF90_PUT_VAR(id_nc, id_rts, vrts)
      CALL err(ierr)
      ierr = NF90_PUT_VAR(id_nc, id_rts2, vrts2)
      CALL err(ierr)
!$OMP END MASTER
!$OMP BARRIER
      
      IF ( tweak_tmean /= 0 ) THEN   
!$OMP SECTIONS         
         umean = 0.
!$OMP SECTION
         vmean = 0.
!$OMP SECTION
         wmean = 0.
!$OMP SECTION
         tmean = 0.
!$OMP SECTION
         smean = 0.
!$OMP SECTION
         gmean = 0.
!$OMP SECTION
         rmean = 0.
!$OMP SECTION
         dmean = 0.
!$OMP SECTION
         nmean = 0
!$OMP END SECTIONS
!$OMP BARRIER
      END IF
         
   END IF
   
   IF ( tweak_tmean == 0 ) THEN
      
      IF (lverbose) THEN
         PRINT*,' Using full time resolution of data ',OMP_GET_THREAD_NUM()
      END IF
   
   ELSE IF (tweak_tmean /= 0 .AND. ints /= 1) THEN

!$OMP SECTIONS
!$OMP SECTION      
      umean = umean + uflux
      IF (lverbose) THEN
         PRINT*,'Adding uflux to time average. nmean, thread = ',nmean,OMP_GET_THREAD_NUM()
      END IF
!$OMP SECTION
      vmean = vmean + vflux
      IF (lverbose) THEN
         PRINT*,'Adding vflux to time average. nmean, thread = ',nmean,OMP_GET_THREAD_NUM()
      END IF
!$OMP SECTION
      wmean = wmean + wflux
      IF (lverbose) THEN
         PRINT*,'Adding wflux to time average. nmean, thread = ',nmean,OMP_GET_THREAD_NUM()
      END IF
!$OMP SECTION
      tmean = tmean + tem
      IF (lverbose) THEN
         PRINT*,'Adding tem to time average. nmean, thread = ',nmean,OMP_GET_THREAD_NUM()
      END IF
!$OMP SECTION
      smean = smean + sal
      IF (lverbose) THEN
         PRINT*,'Adding sal to time average. nmean, thread = ',nmean,OMP_GET_THREAD_NUM()
      END IF
!$OMP SECTION
      gmean = gmean + geo
      IF (lverbose) THEN
         PRINT*,'Adding geo to time average. nmean, thread = ',nmean,OMP_GET_THREAD_NUM()
      END IF
!$OMP SECTION
      rmean = rmean + rho
      IF (lverbose) THEN
         PRINT*,'Adding rho to time average. nmean, thread = ',nmean,OMP_GET_THREAD_NUM()
      END IF
!$OMP SECTION
      dmean = dmean + vol
      IF (lverbose) THEN
         PRINT*,'Adding vol to time average. nmean, thread = ',nmean,OMP_GET_THREAD_NUM()
      END IF
!$OMP SECTION
      nmean = nmean + 1
      IF (lverbose) THEN
         PRINT*,'Adding nmean to time average. nmean, thread = ',nmean,OMP_GET_THREAD_NUM()
      END IF
!$OMP END SECTIONS
!$OMP BARRIER
      
      !!
      !! Check if its time to average yet
      !!
      !$OMP BARRIER
      IF ( (tweak_tmean == -1 .AND. ihour == (24-hourstep) ) .OR.              &
         & (tweak_tmean == -2 .AND. ihour == (24-hourstep) .AND.               &
         &  iday == idmax(imon,iyear) )                     .OR.               &
         & (tweak_tmean == -3 .AND. ihour == (24-hourstep) .AND.               &
         &  iday == idmax(imon,iyear) .AND. imon == 12 )    .OR.               &
         &  ints == intsend                                 .OR.               &
         &  nmean == tweak_tmean )  THEN
         
         nav = nav + 1
         
         PRINT*,' Time to average fields. nmean, nav thread = ',& 
         &      nmean,nav,OMP_GET_THREAD_NUM()
!$OMP SECTIONS         
!$OMP SECTION
         umean = umean / FLOAT(nmean)
         IF (lverbose) THEN
            PRINT*,' Averaging uflux after ',nmean,' time steps',OMP_GET_THREAD_NUM()
         END IF
!$OMP SECTION
         vmean = vmean / FLOAT(nmean)
         IF (lverbose) THEN
            PRINT*,' Averaging vflux after ',nmean,' time steps',OMP_GET_THREAD_NUM()
         END IF
!$OMP SECTION
         wmean = wmean / FLOAT(nmean)
         IF (lverbose) THEN
            PRINT*,' Averaging wflux after ',nmean,' time steps',OMP_GET_THREAD_NUM()
         END IF
!$OMP SECTION
         tmean = tmean / FLOAT(nmean)
         IF (lverbose) THEN
            PRINT*,' Averaging tem after ',nmean,' time steps',OMP_GET_THREAD_NUM()
         END IF
!$OMP SECTION
         smean = smean / FLOAT(nmean)
         IF (lverbose) THEN
            PRINT*,' Averaging sal after ',nmean,' time steps',OMP_GET_THREAD_NUM()
         END IF
!$OMP SECTION
         gmean = gmean / FLOAT(nmean)
         IF (lverbose) THEN
            PRINT*,' Averaging geo after ',nmean,' time steps',OMP_GET_THREAD_NUM()
         END IF
!$OMP SECTION
         rmean = rmean / FLOAT(nmean)
         IF (lverbose) THEN
            PRINT*,' Averaging rho after ',nmean,' time steps',OMP_GET_THREAD_NUM()
         END IF
!$OMP SECTION
         dmean = dmean / FLOAT(nmean)
         IF (lverbose) THEN
            PRINT*,' Averaging vol after ',nmean,' time steps',OMP_GET_THREAD_NUM()
         END IF
!$OMP END SECTIONS
!$OMP BARRIER

!$OMP SINGLE
         uflux = umean
         vflux = vmean
         wflux = wmean
         tem   = tmean
         sal   = smean
         geo   = gmean
         rho   = rmean
         vol   = dmean
!$OMP END SINGLE

!$OMP SINGLE
         umean = 0.
         vmean = 0.
         wmean = 0.
         tmean = 0.
         smean = 0.
         gmean = 0.
         rmean = 0.
         dmean = 0.
         nmean = 0
         IF (lverbose) THEN
            PRINT*,' Resetting means. nmean, thread = ',nmean,OMP_GET_THREAD_NUM()
         END IF
!$OMP END SINGLE

      ELSE
      
         IF (lverbose) THEN
            PRINT*,' Time step added to time average. nmean, thread ',nmean,OMP_GET_THREAD_NUM()
         END IF
         CYCLE
      
      END IF
         
   END IF   
   
   IF ( tweak_zmean == 0 ) THEN
      
      IF (lverbose) THEN
         PRINT*,' Using full zonal resolution of data, thread ',OMP_GET_THREAD_NUM()
      END IF
   
   ELSE IF ( (tweak_zmean == -1 .OR. tweak_zmean == -2) .AND. ints /= 1 ) THEN
      
!$OMP SECTIONS
!$OMP SECTION
      umean2 = 0.
!$OMP SECTION
      vmean2 = 0.
!$OMP SECTION
      wmean2 = 0.
!$OMP SECTION
      tmean2 = 0.
!$OMP SECTION
      smean2 = 0.
!$OMP SECTION
      rmean2 = 0.
!$OMP SECTION
      gmean2 = 0.
!$OMP SECTION
      dmean2 = 0.
!$OMP END SECTIONS
!$OMP BARRIER
            
!$OMP DO      
      DO jk=1,KM
         DO jj=1,JMT
            nmean2 = 0
            DO ji=1,IMT
               umean2(jj,jk) = umean2(jj,jk) + uflux(ji,jj,jk)
               vmean2(jj,jk) = vmean2(jj,jk) + vflux(ji,jj,jk)
               wmean2(jj,jk) = wmean2(jj,jk) + wflux(ji,jj,jk)
               tmean2(jj,jk) = tmean2(jj,jk) + tem(ji,jj,jk) 
               smean2(jj,jk) = smean2(jj,jk) + sal(ji,jj,jk)
               rmean2(jj,jk) = rmean2(jj,jk) + rho(ji,jj,jk)
               gmean2(jj,jk) = gmean2(jj,jk) + geo(ji,jj,jk)
               dmean2(jj,jk) = dmean2(jj,jk) + vol(ji,jj,jk)
               nmean2 = nmean2 + 1
            END DO
            umean2(jj,jk) = umean2(jj,jk) / FLOAT(nmean2)
            vmean2(jj,jk) = vmean2(jj,jk) / FLOAT(nmean2)
            wmean2(jj,jk) = wmean2(jj,jk) / FLOAT(nmean2)
            tmean2(jj,jk) = tmean2(jj,jk) / FLOAT(nmean2)
            smean2(jj,jk) = smean2(jj,jk) / FLOAT(nmean2)
            rmean2(jj,jk) = rmean2(jj,jk) / FLOAT(nmean2)
            gmean2(jj,jk) = gmean2(jj,jk) / FLOAT(nmean2)
            dmean2(jj,jk) = dmean2(jj,jk) / FLOAT(nmean2)
         END DO
      END DO
!$OMP END DO
      
   END IF
      
   IF (tweak_zmean == -1 .AND. ints /= 1) THEN
   !! 
   !! Zonal mean
   !!
!$OMP DO
      DO jk=1,KM
         DO jj=1,JMT
            DO ji=1,IMT
               uflux(ji,jj,jk) = umean2(jj,jk)
               vflux(ji,jj,jk) = vmean2(jj,jk)
               wflux(ji,jj,jk) = wmean2(jj,jk)
               tem(ji,jj,jk)   = tmean2(jj,jk)
               sal(ji,jj,jk)   = smean2(jj,jk)
               rho(ji,jj,jk)   = rmean2(jj,jk)
               geo(ji,jj,jk)   = gmean2(jj,jk)
               vol(ji,jj,jk)   = dmean2(jj,jk)
            END DO
         END DO
      END DO
!$OMP END DO
   END IF
   
   IF (tweak_zmean == -2 .AND. ints /= 1) THEN
   !!
   !! Zonal asymmetries
   !!
!$OMP DO
      DO jk=1,KM
         DO jj=1,JMT
            DO ji=1,IMT
               uflux(ji,jj,jk) = uflux(ji,jj,jk) - umean2(jj,jk)
               vflux(ji,jj,jk) = vflux(ji,jj,jk) - vmean2(jj,jk)
               wflux(ji,jj,jk) = wflux(ji,jj,jk) - wmean2(jj,jk)
               tem(ji,jj,jk)   = tem(ji,jj,jk)   - tmean2(jj,jk)
               sal(ji,jj,jk)   = sal(ji,jj,jk)   - smean2(jj,jk)
               rho(ji,jj,jk)   = rho(ji,jj,jk)   - rmean2(jj,jk)
               geo(ji,jj,jk)   = geo(ji,jj,jk)   - gmean2(jj,jk)
               vol(ji,jj,jk)   = vol(ji,jj,jk)   - dmean2(jj,jk)
            END DO
         END DO
      END DO
!$OMP END DO
   END IF
   
!$OMP DO
   DO jk=1,KM
      DO jj=1,JMT
         DO ji=1,IMT
            
            IF (tem(ji,jj,jk) /= undef) THEN
            
               ! Dry static energy at layer mid-points
               dse(ji,jj,jk) = cp * tem(ji,jj,jk) + geo(ji,jj,jk) ![J/kg]
               ! Moist static energy
               mse(ji,jj,jk) = dse(ji,jj,jk) + Lv * sal(ji,jj,jk) ![J/kg]
               ! Specific volume
               alpha(ji,jj,jk) = Rd * tem(ji,jj,jk) / rho(ji,jj,jk) ![m3/kg]
               
               !if (ji == 112 .and. jj == 49) then
               !   print*,'jk,tem,sal,geo,',jk,tem(ji,jj,jk),sal(ji,jj,jk),geo(ji,jj,jk)
               !end if
               !
               ! Rescale to proper units  
               !
               geo(ji,jj,jk) = geo(ji,jj,jk) / 1000.       ![m2/s2] -> [km m/s2]
               sal(ji,jj,jk) = sal(ji,jj,jk) * 1000.       ![kg/kg] -> [g/kg]
               IF (logp) THEN
                  rho(ji,jj,jk) = log(rho(ji,jj,jk))       ![Pa]    -> [ln(Pa)]
               ELSE
                  rho(ji,jj,jk) = rho(ji,jj,jk) * 0.01     ![Pa]    -> [hPa]
               END IF
      
               dse(ji,jj,jk) = dse(ji,jj,jk) / 1000.       ![J/kg]  -> [kJ/kg]
               mse(ji,jj,jk) = mse(ji,jj,jk) / 1000.       ![J/kg]  -> [kJ/kg] 
               
               if (sal(ji,jj,jk) > 30) then 
                   print*,'sal, ji, jj, jk',sal(ji,jj,jk),ji,jj,jk
                   stop
               end if
               if (dse(ji,jj,jk) < 210) then 
                   print*,'dse, ji, jj, jk',dse(ji,jj,jk),ji,jj,jk,tem(ji,jj,jk),geo(ji,jj,jk)
                   stop
               end if
               
            ELSE
               
               dse(ji,jj,jk) = undef
               mse(ji,jj,jk) = undef
               alpha(ji,jj,jk) = undef
               
            END IF
            
         END DO
      END DO
   END DO
!$OMP END DO
   
   if (lverbose) then
   ! Control print values
   !ji = INT(IMT/2)
   !jj = INT(JMT/2)
   ji = 112
   jj = 49
   do jk = 1,km
       print*,'t,q,p,z,dse',tem(ji,jj,jk),sal(ji,jj,jk),rho(ji,jj,jk),geo(ji,jj,jk),dse(ji,jj,jk)
   end do
   !print*,'tem(ji,jj,km)',tem(ji,jj,:)
   !print*,'sal(ji,jj,km)',sal(ji,jj,:)
   !print*,'rho(ji,jj,km)',rho(ji,jj,:)
   !print*,'geo(ji,jj,km)',geo(ji,jj,:)
   !print*,'dse(ji,jj,km)',dse(ji,jj,:)
   endif   
   
   ! We should skip first step if we are doing time means
   ! or we are using time tendencies
   !IF (ints == 1 .or. (tweak_tend /= 0 .and. nav <= 1)) then
   IF (ints == 1 .or. (tweak_tend /= 0 .and. ints <= 1)) then
      IF (lverbose) THEN
         PRINT*,' Cycling',OMP_GET_THREAD_NUM()
      END IF
      CYCLE
   END IF
   
   ALLOCATE (vars(NRST), mm(0:imm), mm2(0:imm))
   vars(:) = 0.
   mm(:) = 0
   mm2(:) = 0
   

PRINT*,' Start main loop '   
!$OMP DO REDUCTION(+:psrr,psrr2,volyr,volrr,psir,psir2,psiz,psiz2,tsryz,tsryr,tsryr2,nyz,nyr,nyr2)
   DO jk=1,KM
      
      ithread = OMP_GET_THREAD_NUM()
      IF (lverbose) THEN
         PRINT*,' Main loop. thread, level = ',ithread,jk
      END IF
      ku = jk-1
      kb = jk+1
         
      IF (ku == 0) THEN
         ku = 1
      END IF
      IF (kb == KM+1) THEN
         kb = KM
      END IF
      
      DO jj=1,JMT
         
         jp = jj+1
         jm = jj-1
         
         IF (jm == 0) THEN
            jm=1
         END IF
         IF (jp == JMT+1) THEN
            jp=JMT
         END IF

         DO ji=1,IMT
               
            im = ji-1
            ip = ji+1
            
            IF (im == 0) THEN
               im = IMT
            END IF
            IF (ip == IMT+1) THEN
               ip = 1
            END IF
             
            IF (LBAS == 1) THEN
               il = 1
            ELSE
               il = izov(ji,jj) !Find the sub-domain of this grid point
            END IF
               
            IF (il /= 0 .AND. tem(ji,jj,jk) /= undef) THEN
            
               if (sal(ji,jj,jk) > 30) then 
                   print*,'sal, ji, jj, jk',sal(ji,jj,jk),ji,jj,jk
                   stop
               end if
               if (dse(ji,jj,jk) < 210) then 
                   print*,'dse, ji, jj, jk',dse(ji,jj,jk),ji,jj,jk
                   stop
               end if
               
               
#ifdef psiyz               
               !!
               !! yz overturning in model coordinates
               !!
               !print*,'psiz(jj,jk,il)',psiz(jj,jk,il)
               psiz(jj,jk,il) = psiz(jj,jk,il) + 0.5 * &
                              & ( vflux(ji,jj-1,jk) + vflux(ji,jj,jk) )
               
               psiz2(jj,jk,il) = psiz2(jj,jk,il) + wflux(ji,jj,jk)
#endif
#ifdef psixy               
               !!
               !! Barotropic stream function
               !!
               psxyu(ji,jj,il) = psxyu(ji,jj,il) - uflux(ji,jj,jk)
               psxyv(ji,jj,il) = psxyv(ji,jj,il) + vflux(ji,jj,jk)
#endif
#ifdef psixz               
               !!
               !! xz overturning
               !!
               psxzu(ji,jk,il) = psxzu(ji,jk,il) + 0.5 * &
               &                 ( uflux(ji,jj,jk) + uflux(im,jj,jk) )
               psxzw(ji,jk,il) = psxzw(ji,jk,il) + 0.5 * &
               &                 ( wflux(ji,jj,jk) + wflux(ji,jj,jk-1) )
#endif
               
               !!
               !! Find grid point ir
               !!
               DO jzr=1,NRST
                  !print*,'jzr',jzr
                  SELECT CASE(jzr)
                  
                     CASE(1)
                        !print*,'sal',ji,jj,jk,sal(ji,jj,jk),smin
                        mm(0) = NINT((sal(ji,jj,jk)-smin)/ds) + 1 
                        !print*,'sal mm ',mm(0)
                        mm(1) = NINT((sal(ip,jj,jk)-smin)/ds) + 1 
                        mm(2) = NINT((sal(im,jj,jk)-smin)/ds) + 1 
                        mm(3) = NINT((sal(ji,jp,jk)-smin)/ds) + 1 
                        mm(4) = NINT((sal(ji,jm,jk)-smin)/ds) + 1 
                        mm(5) = NINT((sal(ji,jj,ku)-smin)/ds) + 1 
                        mm(6) = NINT((sal(ji,jj,kb)-smin)/ds) + 1
                        IF (tweak_tend /= 0) THEN
                           mm(7) = NINT((sal2(ji,jj,jk)-smin)/ds) + 1
                        END IF
                        vars(1) = sal(ji,jj,jk)
                     CASE(2)
                        mm(0) = NINT((dse(ji,jj,jk)-dmin)/dd) + 1 
                        mm(1) = NINT((dse(ip,jj,jk)-dmin)/dd) + 1 
                        mm(2) = NINT((dse(im,jj,jk)-dmin)/dd) + 1 
                        mm(3) = NINT((dse(ji,jp,jk)-dmin)/dd) + 1 
                        mm(4) = NINT((dse(ji,jm,jk)-dmin)/dd) + 1 
                        mm(5) = NINT((dse(ji,jj,ku)-dmin)/dd) + 1 
                        mm(6) = NINT((dse(ji,jj,kb)-dmin)/dd) + 1
                        IF (tweak_tend /= 0) THEN
                           mm(7) = NINT((dse2(ji,jj,jk)-dmin)/dd) + 1
                        END IF
                        vars(2) = dse(ji,jj,jk)
                        !do jn = 0,imm
                        !   if (mm(jn) < 0) then
                        !   print*,'mm < 0',jn,ji,jj,jk,mm(jn),dse(ji,jj,jk),geo(ji,jj,jk),rho(ji,jj,jk),tem(ji,jj,jk)
                        !   print*,dse(ji,jj,:)
                        !   !stop
                        !   endif
                        !enddo
                        
                     CASE(3)
                        mm(0) = NINT((mse(ji,jj,jk)-mmin)/dm) + 1 
                        mm(1) = NINT((mse(ip,jj,jk)-mmin)/dm) + 1 
                        mm(2) = NINT((mse(im,jj,jk)-mmin)/dm) + 1 
                        mm(3) = NINT((mse(ji,jp,jk)-mmin)/dm) + 1 
                        mm(4) = NINT((mse(ji,jm,jk)-mmin)/dm) + 1 
                        mm(5) = NINT((mse(ji,jj,ku)-mmin)/dm) + 1 
                        mm(6) = NINT((mse(ji,jj,kb)-mmin)/dm) + 1
                        IF (tweak_tend /= 0) THEN
                           mm(7) = NINT((mse2(ji,jj,jk)-mmin)/dm) + 1
                        END IF
                        vars(3) = mse(ji,jj,jk)
                     CASE(4)
                        mm(0) = NINT((alpha(ji,jj,jk)-amin)/da) + 1 
                        mm(1) = NINT((alpha(ip,jj,jk)-amin)/da) + 1 
                        mm(2) = NINT((alpha(im,jj,jk)-amin)/da) + 1 
                        mm(3) = NINT((alpha(ji,jp,jk)-amin)/da) + 1 
                        mm(4) = NINT((alpha(ji,jm,jk)-amin)/da) + 1 
                        mm(5) = NINT((alpha(ji,jj,ku)-amin)/da) + 1 
                        mm(6) = NINT((alpha(ji,jj,kb)-amin)/da) + 1
                        IF (tweak_tend /= 0) THEN
                           mm(7) = NINT((alpha2(ji,jj,jk)-amin)/da) + 1
                        END IF
                        vars(4) = alpha(ji,jj,jk)
                     CASE(5)
                        mm(0) = NINT((rho(ji,jj,jk)-rmin)/dr) + 1 
                        mm(1) = NINT((rho(ip,jj,jk)-rmin)/dr) + 1 
                        mm(2) = NINT((rho(im,jj,jk)-rmin)/dr) + 1 
                        mm(3) = NINT((rho(ji,jp,jk)-rmin)/dr) + 1 
                        mm(4) = NINT((rho(ji,jm,jk)-rmin)/dr) + 1 
                        mm(5) = NINT((rho(ji,jj,ku)-rmin)/dr) + 1 
                        mm(6) = NINT((rho(ji,jj,kb)-rmin)/dr) + 1
                        IF (tweak_tend /= 0) THEN
                           mm(7) = NINT((rho2(ji,jj,jk)-rmin)/dr) + 1
                        END IF
                        vars(5) = rho(ji,jj,jk)
                     CASE(6)
                        mm(0) = NINT((tem(ji,jj,jk)-tmin)/dt) + 1 
                        mm(1) = NINT((tem(ip,jj,jk)-tmin)/dt) + 1 
                        mm(2) = NINT((tem(im,jj,jk)-tmin)/dt) + 1 
                        mm(3) = NINT((tem(ji,jp,jk)-tmin)/dt) + 1 
                        mm(4) = NINT((tem(ji,jm,jk)-tmin)/dt) + 1 
                        mm(5) = NINT((tem(ji,jj,ku)-tmin)/dt) + 1 
                        mm(6) = NINT((tem(ji,jj,kb)-tmin)/dt) + 1
                        IF (tweak_tend /= 0) THEN
                           mm(7) = NINT((tem2(ji,jj,jk)-tmin)/dt) + 1
                        END IF
                        vars(6) = tem(ji,jj,jk)
                        
                  END SELECT
                  
                  !print*,'mm',mm
                  DO jn=0,imm
                     mm(jn) = MAX(1 ,mm(jn))
                     mm(jn) = MIN(MR,mm(jn))  
                  END DO 
                                                      
                  !!
                  !! Atmospheric variables in pressure coordinates
                  !!
                  ir2 = NINT( (rho(ji,jj,jk)-rmin)/dr2 ) + 1
                  ir2 = MAX (1, ir2)
                  ir2 = MIN (MR2,ir2)
                  
                  tsryr(ji,jj,ir2,1:NRST) = tsryr(ji,jj,ir2,1:NRST) + &
                  &                         vars(1:NRST)
                  nyr(ji,jj,ir2) = nyr(ji,jj,ir2) + 1
                  
                  
                  !!
                  !! Zonal-mean distributions of variables is various coordinates
                  !!
                  tsryr2(jj,mm(0),jzr,1:NRST) = tsryr2(jj,mm(0),jzr,1:NRST) + vars(1:NRST)
                  nyr2(jj,mm(0),jzr) = nyr2(jj,mm(0),jzr) + 1   
 
                  !!
                  !! Zonal overturning
                  !!
#ifdef psixr                  
                  psxru(ji,mm(0),il,jzr) = psxru(ji,mm(0),il,jzr) + 0.5 * &
                  &                 ( uflux(ji,jj,jk) + uflux(im,jj,jk) )
                  psxrw(ji,mm(0),il,jzr) = psxrw(ji,mm(0),il,jzr) + 0.5 * &
                  &                 ( wflux(ji,jj,jk) + wflux(ji,jj,jk-1) )
#endif
#ifdef psiyr                  
                  !!
                  !! Overturning with generalized vertical coordinate (psir)
                  !! and with generalized horizontal coordinate with vertical
                  !! model coordinates (pszr).
                  !!
                  psir(jj,mm(0),il,jzr) = psir(jj,mm(0),il,jzr) + 0.5 * &
                  &                  ( vflux(ji,jj,jk) + vflux(ji,jj-1,jk) )
#endif
#ifdef psizr
                  pszr(mm(0),jk,il,jzr) = pszr(mm(0),jk,il,jzr) + 0.5 * &
                  &                  ( wflux(ji,jj,jk) + wflux(ji,jj,jk-1) )
#endif
#ifdef psiyr                  
                  !!
                  !! Grid box mass in y-r coordinates
                  !!                  
                  volyr(jj,mm(0),il,jzr) = volyr(jj,mm(0),il,jzr) + vol(ji,jj,jk)
#endif
#ifdef psirr                  
                  DO jzr2=1,NRST
                     
                     SELECT CASE(jzr2)
                     
                     CASE(1)
                        mm2(0) = NINT((sal(ji,jj,jk)-smin)/ds) + 1 
                        mm2(1) = NINT((sal(ip,jj,jk)-smin)/ds) + 1 
                        mm2(2) = NINT((sal(im,jj,jk)-smin)/ds) + 1 
                        mm2(3) = NINT((sal(ji,jp,jk)-smin)/ds) + 1 
                        mm2(4) = NINT((sal(ji,jm,jk)-smin)/ds) + 1 
                        mm2(5) = NINT((sal(ji,jj,ku)-smin)/ds) + 1 
                        mm2(6) = NINT((sal(ji,jj,kb)-smin)/ds) + 1
                        IF (tweak_tend /= 0) THEN
                           mm2(7) = NINT((sal2(ji,jj,jk)-smin)/ds) + 1
                        END IF
                        
                     CASE(2)
                        mm2(0) = NINT((dse(ji,jj,jk)-dmin)/dd) + 1 
                        mm2(1) = NINT((dse(ip,jj,jk)-dmin)/dd) + 1 
                        mm2(2) = NINT((dse(im,jj,jk)-dmin)/dd) + 1 
                        mm2(3) = NINT((dse(ji,jp,jk)-dmin)/dd) + 1 
                        mm2(4) = NINT((dse(ji,jm,jk)-dmin)/dd) + 1 
                        mm2(5) = NINT((dse(ji,jj,ku)-dmin)/dd) + 1 
                        mm2(6) = NINT((dse(ji,jj,kb)-dmin)/dd) + 1
                        IF (tweak_tend /= 0) THEN
                           mm2(7) = NINT((dse2(ji,jj,jk)-dmin)/dd) + 1
                        END IF
                        
                     CASE(3)
                        mm2(0) = NINT((mse(ji,jj,jk)-mmin)/dm) + 1 
                        mm2(1) = NINT((mse(ip,jj,jk)-mmin)/dm) + 1 
                        mm2(2) = NINT((mse(im,jj,jk)-mmin)/dm) + 1 
                        mm2(3) = NINT((mse(ji,jp,jk)-mmin)/dm) + 1 
                        mm2(4) = NINT((mse(ji,jm,jk)-mmin)/dm) + 1 
                        mm2(5) = NINT((mse(ji,jj,ku)-mmin)/dm) + 1 
                        mm2(6) = NINT((mse(ji,jj,kb)-mmin)/dm) + 1
                        IF (tweak_tend /= 0) THEN
                           mm2(7) = NINT((mse2(ji,jj,jk)-mmin)/dm) + 1
                        END IF
                        
                     CASE(4)
                        mm2(0) = NINT((alpha(ji,jj,jk)-amin)/da) + 1 
                        mm2(1) = NINT((alpha(ip,jj,jk)-amin)/da) + 1 
                        mm2(2) = NINT((alpha(im,jj,jk)-amin)/da) + 1 
                        mm2(3) = NINT((alpha(ji,jp,jk)-amin)/da) + 1 
                        mm2(4) = NINT((alpha(ji,jm,jk)-amin)/da) + 1 
                        mm2(5) = NINT((alpha(ji,jj,ku)-amin)/da) + 1 
                        mm2(6) = NINT((alpha(ji,jj,kb)-amin)/da) + 1 
                        IF (tweak_tend /= 0) THEN
                           mm2(7) = NINT((alpha2(ji,jj,jk)-amin)/da) + 1   
                        END IF  
                     
                     CASE(5)
                        mm2(0) = NINT((rho(ji,jj,jk)-rmin)/dr) + 1 
                        mm2(1) = NINT((rho(ip,jj,jk)-rmin)/dr) + 1 
                        mm2(2) = NINT((rho(im,jj,jk)-rmin)/dr) + 1 
                        mm2(3) = NINT((rho(ji,jp,jk)-rmin)/dr) + 1 
                        mm2(4) = NINT((rho(ji,jm,jk)-rmin)/dr) + 1 
                        mm2(5) = NINT((rho(ji,jj,ku)-rmin)/dr) + 1 
                        mm2(6) = NINT((rho(ji,jj,kb)-rmin)/dr) + 1
                        IF (tweak_tend /= 0) THEN
                           mm2(7) = NINT((rho2(ji,jj,jk)-rmin)/dr) + 1
                        END IF
                     
                     CASE(6)
                        mm2(0) = NINT((tem(ji,jj,jk)-tmin)/dt) + 1 
                        mm2(1) = NINT((tem(ip,jj,jk)-tmin)/dt) + 1 
                        mm2(2) = NINT((tem(im,jj,jk)-tmin)/dt) + 1 
                        mm2(3) = NINT((tem(ji,jp,jk)-tmin)/dt) + 1 
                        mm2(4) = NINT((tem(ji,jm,jk)-tmin)/dt) + 1 
                        mm2(5) = NINT((tem(ji,jj,ku)-tmin)/dt) + 1 
                        mm2(6) = NINT((tem(ji,jj,kb)-tmin)/dt) + 1
                        IF (tweak_tend /= 0) THEN
                           mm2(7) = NINT((tem2(ji,jj,jk)-tmin)/dt) + 1
                        END IF
                      
                     END SELECT
                     
                     
                     DO jn=0,imm
                        mm2(jn) = MAX(1 ,mm2(jn))
                        mm2(jn) = MIN(MR,mm2(jn))  
                     END DO   
                     
                     !!
                     !! Stream function using two generalized coordinates (psrr)
                     !!
                     !! The mass flux between two adjacent grid boxes in space or
                     !! time is treated as a flux between two grid boxes in 
                     !! r-r space (not necessarily adjacent). 
                     !! The flux is put along a linear line between the two boxes
                     !!      
                     DO jn=1,6  ! loop over directions
                        IF (mm2(jn) > mm2(0)) THEN  
                           
                           !! Calculate slope
                           !! Warning: This may be infinit if mm(jn) = mm(0)
                           !mk = (DBLE(mm2(jn))-DBLE(mm2(0))) / (DBLE(mm(jn))-DBLE(mm(0)))
                           
                           DO m1=mm2(0),mm2(jn)-1 
                              
                              !! Calculate jzr value along slope
                              !ir = INT( (DBLE(m1)-DBLE(mm2(0)))/mk + DBLE(mm(0)) )
                              ir = INT( mm(0) ) ! step along constant at mm(0)
                              ir = max(ir,1)
                              ir = min(ir,mr)
                              
                              SELECT CASE(jn)
                              
                              CASE(1)
                                 psrr(ir,m1,il,jzr,jzr2) = &
                              &  psrr(ir,m1,il,jzr,jzr2) + uflux(ji,jj,jk)
                              CASE(2)
                                 psrr(ir,m1,il,jzr,jzr2) = &
                              &  psrr(ir,m1,il,jzr,jzr2) - uflux(im,jj,jk)
                              CASE(3)
                                 psrr(ir,m1,il,jzr,jzr2) = &
                              &  psrr(ir,m1,il,jzr,jzr2) + vflux(ji,jj,jk)
                              CASE(4)
                                 psrr(ir,m1,il,jzr,jzr2) = &
                              &  psrr(ir,m1,il,jzr,jzr2) - vflux(ji,jm,jk)
                              CASE(5)
                                 psrr(ir,m1,il,jzr,jzr2) = &
                              &  psrr(ir,m1,il,jzr,jzr2) - wflux(ji,jj,ku)
                              CASE(6)
                                 psrr(ir,m1,il,jzr,jzr2) = &
                              &  psrr(ir,m1,il,jzr,jzr2) + wflux(ji,jj,jk)
                              END SELECT
                              
                           END DO
                        END IF
                     END DO
                     
                     IF (tweak_tend /= 0) THEN
                        !! Calculate slope
                        !mk = (DBLE(mm2(0))-DBLE(mm2(7))) / (DBLE(mm(0))-DBLE(mm(7)))
                        
                        IF (mm2(0) > mm2(7)) THEN
                           DO m1=mm2(7),mm2(0)-1
                              !! Calculate jzr value
                              !ir = INT( (DBLE(m1)-DBLE(mm2(7)))/mk + DBLE(mm(7)) )
                              ir = INT( mm(0) )
                              ir = max(ir,1)
                              ir = min(ir,mr)
                              !print*,'incr',ir,m1
                              psrr2(ir,m1,il,jzr,jzr2) = &
                           &  psrr2(ir,m1,il,jzr,jzr2) + &
                           &  vol(ji,jj,jk) / (FLOAT(hourstep) * 3600.)
                           END DO
                        ELSE IF (mm2(0) < mm2(7)) THEN
                           DO m1=mm2(0),mm2(7)-1
                              !! Calculate jzr value
                              !ir = INT( (DBLE(m1)-DBLE(mm2(7)))/mk + DBLE(mm(7)) )
                              ir = INT( mm(0) )
                              ir = max(ir,1)
                              ir = min(ir,mr)
                              !print*,'decr',ir,m1
                              psrr2(ir,m1,il,jzr,jzr2) = &
                           &  psrr2(ir,m1,il,jzr,jzr2) - &
                           &  vol(ji,jj,jk) / (FLOAT(hourstep) * 3600.)
                           END DO
                        END IF
                     END IF
                     
                     
                     
                     !!
                     !! Mass distribution in generalized coordinates
                     !!
                     volrr(mm(0),mm2(0),il,jzr,jzr2) = volrr(mm(0),mm2(0),il,jzr,jzr2)+&
                     &                             vol(ji,jj,jk)
                  
                  END DO
#endif                  
               END DO
               
               tsryz(jj,jk,il,1:NRST) = tsryz(jj,jk,il,1:NRST) + vars(1:NRST)
               nyz(jj,jk,il) = nyz(jj,jk,il) + 1 
               
            END IF
         END DO    !end i loop
      END DO       !end j loop
   END DO          !end k loop
!$OMP END DO 
!$OMP BARRIER   
   DEALLOCATE (vars, mm, mm2)     
   
   
   !!
   !! Compute the average and mask where there is no data
   !!
   IF (lverbose) THEN
      PRINT*,' Perform averages and mask empty points. Thread = ',OMP_GET_THREAD_NUM()
   END IF
!$OMP DO 
   DO jzr=1,NRST
      
#ifdef psiyr
      DO jk=1,MR2
         DO jj=1,JMT
            DO ji=1,IMT
               IF (nyr(ji,jj,jk) /= 0) THEN
                  tsryr(ji,jj,jk,jzr) = DBLE(tsryr(ji,jj,jk,jzr)) / DBLE(nyr(ji,jj,jk))
               ELSE
                  tsryr(ji,jj,jk,jzr) = 0.
               END IF
            END DO
         END DO
      END DO
#endif
#ifdef psiyz
      DO jl=1,LBAS
         DO jk=1,KM
            DO jj=1,JMT
               IF (nyz(jj,jk,jl) /= 0) THEN
                  tsryz(jj,jk,jl,jzr) = DBLE(tsryz(jj,jk,jl,jzr)) / DBLE(nyz(jj,jk,jl))
               ELSE
                  tsryz(jj,jk,jl,jzr) = 0.
               END IF
            END DO
         END DO
      END DO
#endif
#ifdef psiyr
      DO jzr2=1,NRST
         DO jk=1,MR
            DO jj=1,JMT
               IF (nyr2(jj,jk,jzr) /= 0) THEN
                  tsryr2(jj,jk,jzr,jzr2) = DBLE(tsryr2(jj,jk,jzr,jzr2)) &
                  & / DBLE(nyr2(jj,jk,jzr))
               ELSE
                  tsryr2(jj,jk,jzr,jzr2) = 0.
               END IF
            END DO
         END DO
      END DO
#endif
   END DO
!$OMP END DO  
   
!   ! Integrate psi(y,r)
!$OMP SINGLE
#ifdef psiyr   
   DO jzr=1,NRST
      DO jk=MR-1,1,-1
         psir(:,jk,:,jzr) = psir(:,jk,:,jzr) + psir(:,jk+1,:,jzr)
      END DO
   END DO
#endif
!   
!   jzr = 3 !MSE
!   psir2(:,:,1,jzr) = psir(:,:,1,jzr) * 10**(-9)
!   psir2(:,:,2,jzr) = psir(:,:,1,jzr) * 10**(-9)
!   DO jk=1,MR
!      DO jj=1,JMT
!         IF ( psir2(jj,jk,1,jzr) < 30. ) THEN
!            psir2(jj,jk,1,jzr) = 0.
!         END IF
!         print*,jj,jk,psir2(jj,jk,1,jzr)
!         IF ( psir2(jj,jk,2,jzr) > -30. ) THEN
!            psir2(jj,jk,2,jzr) = 0.
!         END IF
!      END DO
!   END DO
!$OMP END SINGLE   

!$OMP SINGLE
#ifdef psiyz
   DO jk=2,KM
      psiz(:,jk,:) = psiz(:,jk,:) + psiz(:,jk-1,:)
   END DO
#endif
!$OMP END SINGLE

!$OMP SINGLE
#ifdef psiyz
   DO jj=2,JMT
      psiz2(jj,:,:) = psiz2(jj,:,:) + psiz2(jj-1,:,:)
   END DO
#endif
!$OMP END SINGLE

   
!$OMP SINGLE
!    Find zonal mean specific humidity               
!   DO jzr=1,NRST
!      DO jr=1,MR
!         DO jj=2,JMT
!            lh = tsryr2(jj,jr,jzr,1)
!            IF (lh /= undef) THEN
!                              
!               i_lh = INT ((lh-smin)/ds)
!               print*,lh,smin,i_lh
!               
!               psrr2(i_lh,jr,1,1,jzr) = psrr2(i_lh,jr,1,1,jzr) + &
!             & psir2(jj,jr,1,jzr) - psir2(jj-1,jr,1,jzr)
!               print*,lh,i_lh,psrr2(i_lh,jr,1,jzr)
!               print*,jj,jr,jzr,psir2(jj,jr,1,jzr)
!            END IF
!         END DO
!      END DO
!   END DO   
!$OMP END SINGLE
   

!$OMP BARRIER   
   !!
   !! Save output to monthly averages
   !!
   nsteps = nsteps + 1
   IF (lverbose) THEN
      PRINT*,nsteps,' steps done by thread ',OMP_GET_THREAD_NUM()
   END IF
   
!$OMP SECTIONS
#ifdef psixy
!$OMP SECTION
   psxyu_av = psxyu_av + psxyu
!$OMP SECTION
   psxyv_av = psxyv_av + psxyv
#endif
#ifdef psixz
!$OMP SECTION
   psxzu_av = psxzu_av + psxzu
!$OMP SECTION   
   psxzw_av = psxzw_av + psxzw
#endif
#ifdef psixr
!$OMP SECTION
   psxru_av = psxru_av + psxru
!$OMP SECTION
   psxrw_av = psxrw_av + psxrw
#endif
#ifdef psiyz
!$OMP SECTION
   psiz_av = psiz_av + psiz
!$OMP SECTION
   psiz2_av = psiz2_av + psiz2
!$OMP SECTION
   tsryz_av = tsryz_av + tsryz
#endif
#ifdef psiyr 
!$OMP SECTION  
   psir_av = psir_av + psir
!$OMP SECTION
   psir2_av = psir2_av + psir2
!$OMP SECTION
   tsryr_av = tsryr_av + tsryr
!$OMP SECTION
   tsryr2_av = tsryr2_av + tsryr2
!$OMP SECTION
   volyr_av = volyr_av + volyr
#endif
#ifdef psizr
!$OMP SECTION
   pszr_av = pszr_av + pszr
#endif
#ifdef psirr
!$OMP SECTION
   psrr_av = psrr_av + psrr
!$OMP SECTION
   psrr2_av = psrr2_av + psrr2
!$OMP SECTION
   volrr_av = volrr_av + volrr
#endif
!$OMP END SECTIONS
!$OMP BARRIER

#ifdef psiyr     
!$OMP DO
   DO jk=1,MR2
      DO jj=1,JMT
         DO ji=1,IMT
            IF (nyr(ji,jj,jk) /= 0) THEN
               nyr_av(ji,jj,jk) = nyr_av(ji,jj,jk) + 1
            END IF
         END DO
      END DO
   END DO
!$OMP END DO
#endif
#ifdef psiyz   
!$OMP DO
   DO jk=1,KM
      DO jl=1,LBAS
         DO jj=1,JMT
            IF (nyz(jj,jk,jl) /= 0) THEN
               nyz_av(jj,jk,jl) = nyz_av(jj,jk,jl) + 1
            END IF
         END DO
      END DO
   END DO
!$OMP END DO
#endif
#ifdef psiyr   
!$OMP DO
   DO jzr2=1,NRST
      DO jk=1,MR
         DO jj=1,JMT
            IF (nyr2(jj,jk,jzr2) /= 0) THEN
               nyr2_av(jj,jk,jzr2) = nyr2_av(jj,jk,jzr2) + 1
            END IF
         END DO
      END DO
   END DO 
!$OMP END DO 
#endif
   
   IF (lverbose) THEN
      PRINT*,'Average data? nsteps, ints2, thread = ',&
      & nsteps, ints2, OMP_GET_THREAD_NUM()
      PRINT*,'Average data?',ihour,iday,imon,iyear
   END IF
   IF ((iday == idmax(imon,iyear) .AND. ihour == (24 - hourstep) .AND. tweak_freq == 0) &
   &.OR. ints == intsend &
   &.OR.(nsteps == tweak_freq .AND. tweak_freq /= 0) ) THEN 

!$OMP BARRIER
!$OMP SINGLE
      ints2 = ints2 + 1
!$OMP END SINGLE
!$OMP BARRIER
      IF (lverbose) THEN
         PRINT*,' Averaging and storing data. ints2, thread = ',&
         &      ints2,OMP_GET_THREAD_NUM() 
      END IF
!$OMP SECTIONS
#ifdef psiyz
!$OMP SECTION
      psiz_av = psiz_av / DBLE(nsteps)
!$OMP SECTION
      psiz2_av = psiz2_av / DBLE(nsteps)
#endif
#ifdef psizr       
!$OMP SECTION
      pszr_av = pszr_av / DBLE(nsteps)    
#endif
#ifdef psiyr
!$OMP SECTION
      psir_av = psir_av / DBLE(nsteps)
!$OMP SECTION    
      psir2_av = psir2_av / DBLE(nsteps)
!$OMP SECTION
      volyr_av = volyr_av / DBLE(nsteps) 
#endif
#ifdef psirr
!$OMP SECTION
      psrr_av = psrr_av / DBLE(nsteps)
!$OMP SECTION
      psrr2_av = psrr2_av / DBLE(nsteps)
!$OMP SECTION
      volrr_av = volrr_av / DBLE(nsteps)
#endif
#ifdef psixy
!$OMP SECTION
      psxyu_av = psxyu_av / DBLE(nsteps)
!$OMP SECTION
      psxyv_av = psxyv_av / DBLE(nsteps)
#endif
#ifdef psixz
!$OMP SECTION
      psxzu_av = psxzu_av / DBLE(nsteps)
!$OMP SECTION
      psxzw_av = psxzw_av / DBLE(nsteps)
#endif
#ifdef psixr
!$OMP SECTION
      psxru_av = psxru_av / DBLE(nsteps)
!$OMP SECTION
      psxrw_av = psxrw_av / DBLE(nsteps)
#endif
!$OMP END SECTIONS
!$OMP BARRIER

!$OMP DO
      DO jzr=1,NRST
#ifdef psiyr
         DO jk=1,MR2
            DO jj=1,JMT
               DO ji=1,IMT
                  IF (nyr_av(ji,jj,jk) /= 0) THEN
                     tsryr_av(ji,jj,jk,jzr) = DBLE(tsryr_av(ji,jj,jk,jzr)) / &
                     &                        DBLE(nyr_av(ji,jj,jk))
                  ELSE
                     tsryr_av(ji,jj,jk,jzr) = undef
                  END IF
               END DO
            END DO
         END DO
#endif
#ifdef psiyz
         DO jl=1,LBAS
            DO jk=1,KM
               DO jj=1,JMT
                  IF( nyz_av(jj,jk,jl) /= 0) THEN
                     tsryz_av(jj,jk,jl,jzr) = DBLE(tsryz_av(jj,jk,jl,jzr)) / &
                     &                        DBLE(nyz_av(jj,jk,jl))
                  ELSE
                     tsryz_av(jj,jk,jl,jzr) = undef
                  END IF
               END DO
            END DO
         END DO
#endif
#ifdef psiyr
         DO jzr2=1,NRST
            DO jk=1,MR
               DO jj=1,JMT
                  IF (nyr2_av(jj,jk,jzr) /= 0) THEN
                     tsryr2_av(jj,jk,jzr,jzr2) = DBLE(tsryr2_av(jj,jk,jzr,jzr2)) /&
                     &                           DBLE(nyr2_av(jj,jk,jzr))
                  ELSE
                     tsryr2_av(jj,jk,jzr,jzr2) = undef
                  END IF
               END DO
            END DO
         END DO
#endif
      END DO
!$OMP END DO

!$OMP BARRIER 
!$OMP MASTER
      
      !!
      !! Write all data to netCDF files
      !!
      PRINT*,' Writing data to netCDF files '
      IF (lverbose) THEN
         PRINT*,' year '
      END IF
      ierr = NF90_PUT_VAR(id_nc, id_year, iyear, start=(/ ints2 /) )
      CALL err(ierr)
      
      IF (lverbose) THEN
         PRINT*,' month '
      END IF
      ierr = NF90_PUT_VAR(id_nc, id_mon, imon, start=(/ints2/))
      CALL err(ierr)
#ifdef psiyz      
      IF (lverbose) THEN
         PRINT*,' psi(y,z) '
      END IF
      ierr  = NF90_PUT_VAR(id_nc, id_psiyzv, psiz_av, start=[1,1,1,ints2],      &
      &                                               count=[JMT,KM,LBAS,1])
      CALL err(ierr)
      ierr  = NF90_PUT_VAR(id_nc, id_psiyzw, psiz2_av, start=[1,1,1,ints2],      &
      &                                                count=[JMT,KM,LBAS,1])
      CALL err(ierr)
#endif
#ifdef psizr   
      IF (lverbose) THEN
         PRINT*,' psi(r,z) '
      END IF
      start(1:5) = [  1,  1,      1,    1, ints2]
      count(1:5) = [ MR, KM,   LBAS, NRST,    1]
      ierr  = NF90_PUT_VAR(id_nc, id_psirz, pszr_av, start=[1,1,1,1,ints2],    &
      &                                              count=[MR,KM,LBAS,NRST,1])
      CALL err(ierr)
#endif
#ifdef psiyr      
      IF (lverbose) THEN
         PRINT*,' psi(y,r) '
      END IF
      ierr  = NF90_PUT_VAR(id_nc, id_psiyr, psir_av, start=[1,1,1,1,ints2],   &
      &                                               count=[JMT,MR,LBAS,NRST,1])
      CALL err(ierr)
      
      IF (lverbose) THEN
         PRINT*,' psi2(y,r) '
      END IF
      ierr  = NF90_PUT_VAR(id_nc, id_psiyr2, psir2_av, start=[1,1,1,1,ints2],   &
      &                                               count=[JMT,MR,LBAS2,NRST,1])
      CALL err(ierr)
#endif
#ifdef psirr   
      IF (lverbose) THEN
         PRINT*,' psi(r,r) '
      END IF
      start(1:6) = [  1,  1,      1,    1,    1, ints2]
      count(1:6) = [ MR, MR,   LBAS, NRST, NRST,    1]
      ierr  = NF90_PUT_VAR(id_nc, id_psirr, psrr_av, start=[1,1,1,1,1,ints2],  &
      &                                              count=[MR,MR,LBAS,NRST,NRST,1])
      CALL err(ierr)
   
      IF (lverbose) THEN
         PRINT*,' psi2(r,r) '
      END IF
      ierr  = NF90_PUT_VAR(id_nc, id_psirr2, psrr2_av, start=[1,1,1,1,1,ints2],  &
      &                                                count=[MR,MR,LBAS,NRST,NRST,1])
      CALL err(ierr)
#endif
#ifdef psixy      
      IF (lverbose) THEN
         PRINT*,' psixy(x,y) '
      END IF
      ierr  = NF90_PUT_VAR(id_nc, id_psixyu, psxyu_av, start=[1,1,1,ints2], &
      &                                                count=[IMT,JMT,LBAS,1] )
      CALL err(ierr)
      ierr  = NF90_PUT_VAR(id_nc, id_psixyv, psxyv_av, start=[1,1,1,ints2], &
      &                                                count=[IMT,JMT,LBAS,1] )
      CALL err(ierr)
#endif
#ifdef psixz      
      IF (lverbose) THEN
         PRINT*,' psixz(x,z) '
      END IF
      ierr  = NF90_PUT_VAR(id_nc, id_psixzu, psxzu_av, start=[1,1,1,ints2], &
      &                                                count=[IMT,KM,LBAS,1] )
      CALL err(ierr)
      ierr  = NF90_PUT_VAR(id_nc, id_psixzw, psxzw_av, start=[1,1,1,ints2], &
      &                                                count=[IMT,KM,LBAS,1] )
      CALL err(ierr)
#endif
#ifdef psixr      
      IF (lverbose) THEN
         PRINT*,' psixr(x,r) '
      END IF
      ierr  = NF90_PUT_VAR(id_nc, id_psixru, psxru_av, start=[1,1,1,1,ints2], &
      &                                                count=[IMT,MR,LBAS,NRST,1] )
      CALL err(ierr)
      ierr  = NF90_PUT_VAR(id_nc, id_psixrw, psxrw_av, start=[1,1,1,1,ints2], &
      &                                                count=[IMT,KM,LBAS,NRST,1] )
      CALL err(ierr)
#endif
#ifdef psiyz      
      IF (lverbose) THEN
         PRINT*,' r(y,z) '
      END IF
      start(1:5) = [  1,  1,      1,    1, ints2]
      count(1:5) = [JMT, KM,   LBAS, NRST,    1]
      ierr  = NF90_PUT_VAR(id_nc, id_rtsyz, tsryz_av, start=[1,1,1,1,ints2],   &
      &                                               count=[JMT,KM,LBAS,NRST,1])
      CALL err(ierr)
#endif
#ifdef psiyr   
      IF (lverbose) THEN
         PRINT*,' r(x,y,r) '
      END IF
      start(1:5) = [  1,   1,   1,    1, ints2]
      count(1:5) = [IMT, JMT, MR2, NRST,    1]
      ierr  = NF90_PUT_VAR(id_nc, id_rtsyr, tsryr_av, start=[1,1,1,1,ints2],   &
      &                                               count=[IMT,JMT,MR2,NRST,1])
      CALL err(ierr)
#endif
#ifdef psiyr   
      IF (lverbose) THEN
         PRINT*,' r(y,r) '
      END IF
      start(1:5) = [  1,   1,    1,    1, ints2]
      count(1:5) = [JMT,  MR, NRST, NRST,     1]
      ierr  = NF90_PUT_VAR(id_nc, id_rtsyr2, tsryr2_av, start=[1,1,1,1,ints2], &
      &                                                 count=[JMT,MR,NRST,NRST,1])
      CALL err(ierr)
#endif
#ifdef psiyr   
      IF (lverbose) THEN
         PRINT*,' vol(y,r) '
      END IF
      start(1:5) = [  1,  1,    1,    1, ints2]
      count(1:5) = [JMT, MR, LBAS, NRST,     1]
      ierr  = NF90_PUT_VAR(id_nc, id_volyr, volyr_av, start=[1,1,1,1,ints2],   &
      &                                             count=[JMT,MR,LBAS,NRST,1])
      CALL err(ierr)
#endif
#ifdef psirr      
      IF (lverbose) THEN
         PRINT*,' vol(r,r) '
      END IF
      ierr  = NF90_PUT_VAR(id_nc, id_volrr, volrr_av, start=[1,1,1,1,1,ints2], &
      &                                             count=[MR,MR,LBAS,NRST,NRST,1])
      CALL err(ierr)
#endif
!$OMP END MASTER
!$OMP BARRIER 
      IF (lverbose) THEN
         PRINT*,' Wrote all data to netCDF files', OMP_GET_THREAD_NUM()
      END IF
      
   END IF
   
!$OMP MASTER   
   IF ( (ihour == (24-hourstep) .AND. iday == idmax(imon,iyear) .AND. imon == 12) .OR. &
   &    ints == intsend  ) THEN
      ierr = NF90_CLOSE(id_nc)
      CALL err(ierr)
   END IF
!$OMP END MASTER
!$OMP BARRIER   
   
   END DO !end time loop
   
   !!
   !! End of OpenMP
   !!
   IF (lverbose) THEN
      PRINT*,'End of OpenMP parallelization ',OMP_GET_THREAD_NUM()
   END IF
!$OMP END PARALLEL
   
   
   PRINT*,'--------------------------------------------------------------------'
   PRINT*,'                                                                    '
   PRINT*,' Output will be stored in netCDF file                               '
   PRINT 102, TRIM(ncFile)
   PRINT*,' in directory                                                       '
   PRINT 103, TRIM(outDataDir)

!!------------------------------------------------------------------------------

   CALL system('date')
   PRINT*,'Fin!'
   STOP

!!------------------------------------------------------------------------------


END PROGRAM PSIXX_IFS
!__________________________________________________________________________________________


!!------------------------------------------------------------------------------
!!------------------------------------------------------------------------------
!!------------------------------------------------------------------------------
!!------------------------------------------------------------------------------
!!------------------------------------------------------------------------------
!!------------------------------------------------------------------------------
!!------------------------------------------------------------------------------


!_______________________________________________________________________
