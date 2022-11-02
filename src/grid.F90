module mod_grid
   
   use mod_vars
   
   implicit none
   
   integer, private :: ji,ip,im,jj,jp,jm,jk,ku,kb,il
   
contains
   
   
   subroutine calc_dxdy()
      
      real*8 :: dlon, dlat
      
      !! initialise
      dx(:,:) = 0.
      dy(:,:) = 0.
      dxdy(:,:) = 0.
      vlat(:) = 0.
      
      !! 
      !! dx at v-points, dy and u-points, dxdy at t-points 
      !!
      do jj=0,JMT
         
         jp = jj
         jm = jj-1
         
         if (jj /= 0) then
            dlat = vlat2(jp) - vlat2(jm)
            vlat(jj) = 0.5 * ( vlat2(jp) + vlat2(jm) )
         end if
         
         do ji=1,IMT
            
            im = ji-1
            if (im == 0) then
               im = IMT
            end if
            
            dlon = vlon(ji) - vlon(im)
            if (dlon < 0.) then
               dlon = dlon + 360.
            end if
            if (dlon > 360.) then
               dlon = dlon - 360.
            end if
            
            dx(ji,jj) = abs( dlon * deg * dcos( vlat2(jp) * radian ) )
            if (jj /= 0) then
               dy(ji,jj) = abs( dlat * deg )
               dxdy(ji,jp) = 0.5 * (dx(ji,jp)+dx(ji,jm)) * dy(ji,jp)
            end if
            
         end do
         
      end do
      
      if (lverbose) then
         print*,'dx(1,:)',dx(1,:)
         print*,'dy(1,:)',dy(1,:)
         print*,'dxdy(1,:)',dxdy(1,:)
      end if
      
   end subroutine
   
   
   subroutine a2cgrid_t(adata,cdata)
   !!-------------------------------------------------------------------------------
   !!
   !!   Subroutine to re-map A-grid to t-points on a C-grid 
   !!
   !!-------------------------------------------------------------------------------
      
      real*4, dimension(imt,1:jmt)              ::  tmp
      real*4, dimension(imt,0:jmt), intent(in)  ::  adata
      real*4, dimension(imt,1:jmt), intent(out) ::  cdata
      
      !! Interpolate first to u-points
      tmp(1:imt,1:jmt) = 0.5 * (adata(1:imt,1:jmt) + adata(1:imt,0:jmt-1))
      !! Interpolate to t-points
      cdata(2:imt,1:jmt) = 0.5 * (tmp(2:imt,1:jmt) + tmp(1:imt-1,1:jmt))
      cdata(1,1:jmt) = 0.5 * (tmp(1,1:jmt) + tmp(imt,1:jmt))
      
   end subroutine a2cgrid_t
   
   subroutine a2cgrid_u(adata,cdata)
   !!-------------------------------------------------------------------------------
   !!
   !!   Subroutine to re-map A-grid u velocity to u-points on a C-grid
   !!
   !!-------------------------------------------------------------------------------
      
      real*4, dimension(imt,0:jmt), intent(in)  ::  adata
      real*4, dimension(imt,1:jmt), intent(out) ::  cdata
      
      cdata(1:imt,1:jmt) = 0.5 * (adata(1:imt,1:jmt) + adata(1:imt,0:jmt-1))
      
   end subroutine
   
   subroutine a2cgrid_v(adata,cdata)
   !!-------------------------------------------------------------------------------
   !!
   !!   Subroutine to re-map A-grid v velocity to v-points on a C-grid
   !!
   !!-------------------------------------------------------------------------------
      
      real*4, dimension(imt,0:jmt), intent(in)  ::  adata
      real*4, dimension(imt,0:jmt), intent(out) ::  cdata
      
      cdata(2:imt,0:jmt) = 0.5 * (adata(2:imt,0:jmt) + adata(1:imt-1,0:jmt))
      cdata(1,0:jmt) = 0.5 * (adata(1,0:jmt) + adata(imt,0:jmt))
      
   end subroutine
   
   subroutine int_geo()
   !!--------------------------------------------------------------------------------
   !!
   !!   Subroutine to calculate geopotential using the hydrostatic equation
   !!
   !!--------------------------------------------------------------------------------
      
      real*4, allocatable, dimension(:,:) :: tv
      
      allocate (tv(imt,jmt))
      
      tv(:,:) = 0.
      
      do jk=km,1,-1
         
         ku = jk-1
         kb = jk
         
         ! Pressure at top of the atmosphere
         ! Should be zero, but can not for practical reasons
         if (jk == 1) then
            do jj=1,jmt
               do ji=1,imt
                  if (rho_i(ji,jj,ku) == 0.) then
                     print*,'top pressure is zero so geopotential becomes infinite'
                     stop
                  end if
               end do
            end do
         end if
         
         do jj = 1, jmt
             do ji = 1, imt
                 ! Virtual temperature at layer mid-point
                 tv(ji,jj) = ( 1.0 + 0.61 * sal(ji,jj,jk) ) * tem(ji,jj,jk)
                 
                 ! Geopotential at top interface k-1
                 geo_i(ji,jj,ku) = geo_i(ji,jj,kb) + Rd * tv(ji,jj) * log(rho_i(ji,jj,kb)/rho_i(ji,jj,ku))
             end do
         end do
      end do
   
   end subroutine
   
   
   subroutine calc_wflux()
   !!-----------------------------------------------------------------------------
   !!
   !!   Subroutine to calculate vertical mass fluxes by mass conservation
   !!
   !!-----------------------------------------------------------------------------
      
      real*4, dimension(imt,jmt) :: du
      
      wflux(:,:,0) = 0.
      
      do jk=1,KM
         
         du(2:imt,:) = uflux(2:imt,:,jk) - uflux(1:imt-1,:,jk)
         du(1,:) = uflux(1,:,jk) - uflux(imt,:,jk)
            
         wflux(:,:,jk) = wflux(:,:,jk-1) - ( du(:,:) +                         &
         &                           vflux(:,1:jmt,jk) - vflux(:,0:jmt-1,jk) + &
         &                          (dzt(:,:,jk,2) - dzt(:,:,jk,1)) *          &
         &                           dxdy(:,:) / dtstep / dg )
         
      end do
      
   end subroutine
   
   
end module