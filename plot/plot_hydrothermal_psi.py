#!/sw/bin/python2.7
##------------------------------------------------------------------------------
##
##
##    Script to plot stream functions using Python 2.7 with PyNGL and PyNIO
##    Plots hydrothermal stream function with transports and turnover time
##
##    Joakim Kjellsson
##    September 2013
##
##
##------------------------------------------------------------------------------

import os, sys, Nio, Ngl
import numpy as np

from scipy import interpolate
from scipy import ndimage

import sub_psi as sp

##------------------------------------------------------------------------------

ncFile = '/nobackup/vagn2/x_joakj/data/ifs/psi/era_interim_1980-1999_origres_av.nc'

data = sp.read_data(ncFile,ix=0,iy=1,ltend=True)

##
## Open file and open PDF output
##
opt = Nio.options() 
nc  = Nio.open_file(ncFile,'r',opt)
wks = Ngl.open_wks('pdf','era_1980')

##
## Set color map
##
res = Ngl.Resources()
res.wkColorMap = 'BlueWhiteOrangeRed'
Ngl.set_values(wks,res)


##
## Plot psi(r,r)
##
vx = data['rts_lev'][0,:]
vy = data['rts_lev'][1,:]
zplot = data['psirr'][0,:,:]
print np.min(zplot.flatten())
sp.plot_psirr(wks,zplot,vx,vy,xlab='Latent heat [kJ/kg]',ylab='Dry static energy [kJ/kg]')

##------------------------------------------------------------------------------

##
## Overturning in various coordinates
##

zplot = data['psiyr']
vx = data['lat'][:]

zplot = []
vy = []
vsurf = []
vcoord = [4,0,1,2]
flipy = [True,True,False,False]
xon = [False,False,False,True]
ylab = ['Pressure [hPa]','LH [kJ/kg]','DSE [kJ/kg]','MSE [kJ/kg]']
data['psiyr'][0,0,:,:] = data['psiyr'][0,0,:,:]*(-1.)
data['psiyr'][0,4,:,:] = data['psiyr'][0,4,:,:]*(-1.)
for jn in vcoord:
   vy.append(data['rts_lev'][jn,:])
   zplot.append(data['psiyr'][0,jn,:,:])
   vsurf.append(data['rtsyz'][0,jn,-1,:])

sp.plot_psiyr(wks,zplot,vx,vy,vsurf=vsurf,flipy=flipy,xon=xon,ylab=ylab)

Ngl.end()
sys.exit()

##------------------------------------------------------------------------------



##------------------------------------------------------------------------------

tmp = data['rtsyz'][0,2,:,:]
tmp = sp.make_tropical_average( data['lat'], tmp )
vy = [ tmp ]
vx = data['lev'][:]

nc = Nio.open_file('/home/x_joakj/tmp/gfdl_mean.nc','r')
a = nc.variables['a'][::-1]
#a = np.append(nc.variables['a_bnds'][-1,-1],a[::-1])
b = nc.variables['b'][::-1]
#b = np.append(nc.variables['b_bnds'][-1,-1],b[::-1])
p0 = nc.variables['p0'].get_value()
ps = nc.variables['ps'][0,:,:]
ps = np.mean(ps,axis=1)
ps = np.mean(ps[40:50])
nc.close()

p = a*p0+b*ps
p = p/100.
print p

vx = [p]

#p(n,k,j,i) = a(k)*p0 + b(k)*ps(n,j,i)

data = sp.read_data('/nobackup/vagn2/x_joakj/data/ifs/psi/era_does_everything_1980-2009_av.nc')
tmp = sp.make_tropical_average( data['lat'], data['rtsyz'][0,2,:,:] )
vy.append( tmp )


nc = Nio.open_file('/home/x_joakj/tmp/era_mean.nc','r')
a = nc.variables['hyam'][:]
b = nc.variables['hybm'][:]
ps = np.exp(nc.variables['LNSP'][0,0,:,:])
ps = np.mean(ps,axis=1)
ps = np.mean(ps[60:84])
nc.close()

p = a+b*ps
p = p/100.
print p

vx.append(p)

sp.plot_profiles(wks,vy,vx,ylab='Pressure [hPa]',xlab='MSE [kJ/kg]',xmin=300.,xmax=380.)




Ngl.end()
sys.exit()

   
##------------------------------------------------------------------------------

##
## Total mass in (r,r) coordinates
##

print ' '
print '---------------------------------------------'
print '  Plotting total mass of air in (r,r) coordinates'

ilbas   = nc.dimensions['lbas']
volrr   = nc.variables['volrr'][0,:,:,0:ilbas,:,:]

title = ['Tropics','Midlatitudes','Polar','Global']

##
## Log mass
##
volln = np.where(volrr[:,:,:,:,:] > 10**(0), np.log10(volrr[:,:,:,:,:]), 0.)

## Set levels
vol_levels = np.array([0, 2, 4, 6, 8, 10, 12, 14, 15, 16, 17, 18])


del res
res = Ngl.Resources()
res.wkColorMap = 'WhiteBlue'
Ngl.set_values(wks,res)

print '---------------------------------------------'
print 'Max color level: '+repr(np.max(vol_levels))
print 'Max mass(r,r): '+repr(np.max(volln))
print '---------------------------------------------'

colour_labels = []
for clev in vol_levels:
   colour_labels.append(repr(clev))


for jc in range(0,1):
   
   for jb in range(0,1):
      
      vx    =  rts_lev[jb,:]
      vy    =  rts_lev[jc,:]
      
      ## Log of mass
      zplot = volln[jc,jb,:,:,:]
            
      del res
      res = Ngl.Resources()
      
      res.nglMaximize                      =  False
      res.nglPaperOrientation              =  'Portrait'
      res.nglFrame                         =  False
      res.nglDraw                          =  False
      
      res.cnLineLabelDensityF              =  2
      res.cnLineLabelBackgroundColor       =  -1
      res.cnLevelSelectionMode             =  'ExplicitLevels'
   
      res.lbLabelBarOn                     =  True
      res.pmLabelBarSide                   =  'Top'
      res.lbOrientation                    =  'Horizontal'
      res.pmLabelBarDisplayMode            =  'Always'
      res.pmLabelBarWidthF                 =  0.3
      res.pmLabelBarHeightF                =  0.04
      res.lbTitleString                    =  'log~B~10~N~(kg)'
      res.lbLabelFontHeightF               =  0.01
      res.lbTitleFontHeightF               =  0.012
      
      res.sfXArray                         =  vx
      res.sfYArray                         =  vy
      
      res.tiXAxisString                    =  rts_label[jb]
      res.tiYAxisString                    =  rts_label[jc]
      
      res.tmYRLabelsOn                     =  True
      res.tmYLLabelsOn                     =  False
      res.tiYAxisSide                      =  'Right'
      
      if (jc != jb):
         
         psi = psirr[jc,jb,-1,:,:]
                  
         for ja in range(-1,0):
            
            #if(np.mod(ja,2) == 0):
            #   res.tiYAxisSide             =  'Left'
            #   res.tmYLLabelsOn            =  True
            #   res.tmYRLabelsOn            =  False
            #   res.pmLabelBarDisplayMode   =  'Always'
            #else:
            #   res.tiYAxisSide             =  'Right'
            #   res.tmYLLabelsOn            =  False
            #   res.tmYRLabelsOn            =  True
            #   res.pmLabelBarDisplayMode   =  'Never'
            
            res.cnFillOn                   =  True
            res.cnLinesOn                  =  False
            res.cnLineLabelsOn             =  False
            res.cnLevels                   =  vol_levels
            res.vpXF                       =  0.55#vxf[ja]
            res.vpWidthF                   =  0.3
            res.vpYF                       =  0.7#vyf[ja]
            res.vpHeightF                  =  0.3
            #res.tiMainString               =  title[ja]
             
            volplot = Ngl.contour(wks,zplot[ja,:,:],res)
            
            ##
            ## Draw hydrothermal stream lines on top
            ##
            res.pmLabelBarDisplayMode      =  'Never'
            res.cnFillOn                   =  False
            res.cnLinesOn                  =  True
            res.cnLineLabelsOn             =  True
            res.cnInfoLabelOn              =  False
            res.cnLevels                   =  psi_levels
            
            psi2 = psi * 10**(-9)
            psiplot = Ngl.contour(wks,psi2,res)
            
            Ngl.overlay(volplot,psiplot)
            Ngl.draw(volplot)
            
            print '-------------------------------------------'
            print 'mass('+rts_label[jb]+','+rts_label[jc]+')'
            print '    '+title[ja]
            print '     Total mass: '+repr(np.sum(volrr[jb,jc,ja,:,:]))
            
            
            tres = Ngl.Resources()
            fhei = 0.03
            tres.txFontHeightF   =  fhei
            Ngl.text_ndc(wks,'b)',0.82,0.67,tres)
            
         Ngl.frame(wks)



##-----------------------------------------------------------------------------

##
## Zonally integrated quantities
##
print ' '
print '---------------------------------------------'
print '  Plotting zonal-mean quantities'

rtsyz  = nc.variables['rtsyz'][0,:,:,:,:]
vlat   = nc.variables['latitude'][:]
ilev   = nc.dimensions['mod_lev']

rtsyz[0,:,:,:]  = rtsyz[0,:,:,:] * 2.5

vlev = np.arange(1,ilev+1)

del res
res = Ngl.Resources()
res.wkColorMap = 'rainbow'
Ngl.set_values(wks,res)

for jc in range(0,rtsyz.shape[0]):
   
   vx = vlat
   vy = vlev
   zplot = rtsyz[jc,-1,:,:]
   
      
   del res
   res = Ngl.Resources()
   
   res.nglMaximize                      =  False
   res.nglPaperOrientation              =  'Portrait'
   res.nglFrame                         =  False
   res.nglDraw                          =  True
   
   res.cnFillOn                         =  True
   res.cnLinesOn                        =  True
   res.cnLineLabelsOn                   =  True
   res.cnLineLabelDensityF              =  2
   res.cnLineLabelBackgroundColor       =  -1
   
   if( jc != 10):
      res.lbLabelBarOn                  =  True
   else:
      res.lbLabelBarOn                  =  True
   
   res.sfXArray                         =  vlat
   res.sfYArray                         =  vlev
   
   res.tiXAxisString                    =  'Latitude'
   res.tiYAxisString                    =  'Model level index'
   res.tiMainString                     =  rts_label[jc]
   
   res.trYReverse                       =  True
   
   res.vpXF                             =  0.15
   res.vpWidthF                         =  0.8
   res.vpYF                             =  0.7
   res.vpHeightF                        =  0.2
   
   Ngl.contour(wks,zplot,res)
   
   Ngl.frame(wks)


##------------------------------------------------------------------------------


def fill_mask(A):
    '''
    interpolate to fill masked values
    '''
    inds = np.arange(A.shape[0])  #Indices
    nums = np.where(A.mask,1,0)   #1 if mask, 0 if value
    good = np.where(nums == 0)    #Values
    bad  = np.where(nums == 1)    #Mask
    B = np.interp(inds,inds[good],A[good]) #Fill mask with interpolated values
    C = np.where(A.mask,B,A)
    return C

##------------------------------------------------------------------------------

##
## Atmospheric variables in y,r coordinates
##
print ' '
print '---------------------------------------------'
print '  Plotting zonal-mean quantities'

rtsyr    = nc.variables['rtsyr_zm'][0,:,:,:,:]
vlat     = nc.variables['latitude'][:]
mra_lev  = nc.variables['rts_lev'][:,:]

del res
res = Ngl.Resources()
res.wkColorMap = 'rainbow'
Ngl.set_values(wks,res)

for jc in range(0,icoord):  #variable
   
   for jb in range(0,icoord): #vertical coordinate
   
      vx = vlat
      vy = mra_lev[jb,:]
      zplot = rtsyr[jc,jb,:,:]      
      
      ## Interpolate fill values
      for jj in range(0,zplot.shape[1]):
         zplot[:,jj] = fill_mask(zplot[:,jj])
         
      del res
      res = Ngl.Resources()
      
      res.nglMaximize                      =  False
      res.nglPaperOrientation              =  'Portrait'
      res.nglFrame                         =  False
      res.nglDraw                          =  True
      
      res.cnFillOn                         =  False
      res.cnLinesOn                        =  True
      res.cnLineLabelsOn                   =  True
      res.cnLineLabelDensityF              =  3
      res.cnLineLabelBackgroundColor       =  -1
      res.cnInfoLabelOn                    =  True
      
      if( jc != 10):
         res.lbLabelBarOn                  =  True
      else:
         res.lbLabelBarOn                  =  True
   
      res.sfXArray                         =  vx
      res.sfYArray                         =  vy
        
      res.tiXAxisString                    =  'Latitude'
      res.tiYAxisString                    =  rts_label[jb]
      res.tiMainString                     =  rts_label[jc]
      
      res.trYReverse                       =  True
      
      res.vpXF                             =  0.15
      res.vpWidthF                         =  0.8
      res.vpYF                             =  0.7
      res.vpHeightF                        =  0.2
      
      Ngl.contour(wks,zplot,res)
      
      Ngl.frame(wks)


##------------------------------------------------------------------------------

##
## Energy and moisture transports
##

print '---------------------------------------------'
print ' Plotting dry static energy transport'
print ' and water vapor transport ...'


if(1):
   
   ## Select x-coord, y-coord
   ix    =  0
   iy    =  1
   
   psi0  =  -20. #[Sv]
   
   ## Set x-axis
   vxx    =  rts_lev[ix,:]
   
   ## Get the stream function above some reference value
   psi   =  psirr[iy,ix,-1,:,:]  * 10**(-9) #[Sv]
   psi   =  psi - psi0 #psi - psi0
   psi   =  np.where(psi < 0., psi, 0.)     #[Sv]
   
   ## Calculate the integral of (psi-psi0)*dy
   dy    =  rts_lev[iy,1]-rts_lev[iy,0]     #[kJ/kg]
   dpsi  =  psi
   dd    =  dpsi * dy * 10**(-3)            #[PW = 10**15 J/s] Peta-Watt
   vyx    =  np.sum(dd[:,:],axis=0)
   ylabelx = 'Dry static energy transport [PW]'
   
   ##
   ## Print maximum/minimum transport of dry static energy
   ##
   maxval = np.max(np.abs(vyx))
   print ' Max transport of dry static energy: '+repr(np.max(np.abs(vyx)))
   print ' at latent heat: '+repr(vxx[np.where(np.abs(vyx) == maxval)[0]][0])
   
   
   ## Set x-axis
   vxy = rts_lev[iy,:]
   
   dx  =  rts_lev[ix,1] - rts_lev[ix,0]     #[kJ/kg]
   dpsi = psi                               #[Sv]
   dd = (-1) * dpsi * dx / (10**3)          #[PW]
   vyy = np.sum(dd[:,:],axis=1)
   ylabely = 'Latent heat transport [PW]'
   yticks = np.array([0., 1., 2., 3., 4., 5.]) * 2.5
   
   ##
   ## Print maximum/minimum transport of dry static energy
   ##
   maxval = np.max(np.abs(vyy))
   print ' Max transport of latent heat: '+repr(np.max(np.abs(vyy)))
   print ' at dry static energy: '+repr(vxy[np.where(np.abs(vyy) == maxval)[0]][0])
   
   
   del res
   
   res = Ngl.Resources()
   
   res.nglMaximize                       =  False
   res.nglPaperOrientation               =  'Portrait'
   res.nglFrame                          =  False
   res.nglDraw                           =  False
   
   res.vpHeightF                         =  0.2
   res.vpWidthF                          =  0.7
   
   res.tiXAxisFontHeightF                =  0.012
   res.tiYAxisFontHeightF                =  0.012
   res.tmXBLabelFontHeightF              =  0.012
   res.tmYLLabelFontHeightF              =  0.012
   
   res.xyLineThicknessF                  =  3.
   
   plots = []
   
   for ji in range(0,2):
      
      if(ji == 0):
         vx = vxx
         vy = vyx
         
         res.vpXF                        =  0.15
         res.vpYF                        =  0.9
         
         res.tiXAxisString               =  rts_label[ix]
         res.tiYAxisString               =  ylabelx
         #res.tiMainString                =  'Transport of '+rts_label[iy]
         
      elif(ji == 1):
         vx = vxy
         vy = vyy
         
         res.vpXF                        =  0.15
         res.vpYF                        =  0.6
         
         res.tiXAxisString               =  rts_label[iy]
         res.tiYAxisString               =  ylabely
         
         res.trYMaxF      =  yticks[-1]
         res.trYMinF      =  yticks[0]
         res.tmYLMode     =  'Explicit'
         res.tmYLValues   =  yticks
         
         ylabs = []
         for yy in yticks:
            ylabs.append(repr(yy))
         
         res.tmYLLabels   =  ylabs
         
      plots.append(Ngl.xy(wks,vx,vy,res))
      
      if(ji == 1):
         
         res.tmXBOn       = False
         res.tmXBLabelsOn = False
         res.tmXBMinorOn  = False
         res.tmYLOn       = False
         res.tmYLLabelsOn = False
         res.tmYLMinorOn  = False
         res.tmYRLabelsOn = True
         res.tmYROn       = True
         res.tmYUseLeft   = False
         res.tiXAxisOn = False
         res.tiYAxisString = 'Freshwater transport [Sv]'
         res.tiYAxisSide = 'Right'
         
         res.tmYRMode     =  'Explicit'
         res.trYMaxF      =  yticks[-1]/2.5
         res.trYMinF      =  yticks[0]/2.5
         res.tmYRValues   =  yticks / 2.5
         
         ylabs = []
         for yy in yticks:
            ylabs.append(repr(yy/2.5))
         
         res.tmYRLabels   =  ylabs
         
         plots.append(Ngl.xy(wks,vx,vy/2.5,res))
   
   Ngl.draw(plots[0])
   
   #anno = Ngl.add_annotation(plots[1],plots[2])
   
   Ngl.draw(plots[1])
   Ngl.draw(plots[2])
         
         
   
   Ngl.frame(wks)


##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

Ngl.end()

##------------------------------------------------------------------------------
##------------------------------------------------------------------------------