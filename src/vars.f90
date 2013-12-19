MODULE mod_vars
!!------------------------------------------------------------------------------
!!
!!   Module that contains variables and subroutines that should be accessible
!!   by all other processes. 
!!   Has data for grid, code settings, and physical constants
!!
!!-------------------------------------------------------------------------------
   
   IMPLICIT none
   
   integer :: imt, jmt, km
   
   real*4, allocatable, dimension (:)       ::  vlon
   real*4, allocatable, dimension (:)     ::  vlat2
   real*4, allocatable, dimension (:)     ::  vlat
   real*4, allocatable, dimension (:)        ::  vlev
   real*4, allocatable, dimension (:,:)   ::  dxdy, dy
   real*4, allocatable, dimension (:,:) ::  dx
   real*4, allocatable, dimension(:,:,:,:)  ::  dzt
   
   REAL*4, ALLOCATABLE, DIMENSION(:,:,:) :: uflux,vflux,wflux,tem,sal,&
   &                                        rho,geo,vol,rho_i,geo_i
   
   REAL*4                                        ::  undef = -99999., dtstep
   
   REAL*8 :: deg, radian
   
   REAL*8                                        ::  radius = 6371229.,        &
   &                                                 pi = 3.141592653589793238,&
   &                                                 dg = 9.81,                &
   &                                                 Rd = 287.05d0,            &
   &                                                 Rv = 461.5,               &
   &                                                 Lv = 2.5d+6,              &
   &                                                 cp = 1004.d0 
   
   INTEGER :: ihour, iday, imon, iyear, yearstart, monstart,daystart,hourstart,&
   &          hourstep, tweak_tmean, tweak_zmean, tweak_tend, tweak_freq,      &
   &          tweak_entropy
   
   INTEGER, DIMENSION(12,1000:3000)              ::  idmax
   
   CHARACTER :: inDataDir*100, tmpDataDir*100, prefix*100, topoDir*100, rstring*350,     &
   &                                                 lonstring*350,   &
   &                                                 latstring*350
   
   LOGICAL  ::  lverbose, llast, up
   
   
   contains
   
      subroutine init_vars()
      
         allocate ( vlon(imt), vlat2(0:jmt), vlat(1:jmt), vlev(km), &
         &          dxdy(imt,jmt), dy(imt,jmt), dx(imt,0:jmt), dzt(imt,jmt,km,2) )
         
         allocate ( uflux(imt,jmt,km), vflux(imt,0:jmt,km), wflux(imt,jmt,0:km), &
         &          tem(imt,jmt,km), sal(imt,jmt,km), rho(imt,jmt,km), &
         &          geo(imt,jmt,km), vol(imt,jmt,km), rho_i(imt,jmt,0:km), geo_i(imt,jmt,0:km) )
      
      end subroutine
   
END MODULE