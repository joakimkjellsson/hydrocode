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


##------------------------------------------------------------------------------

ncFile = '/nobackup/vagn2/x_joakj/data/cnrm/psi/cnrm_historical_1980_1980.nc'

##
## Open file and open PDF output
##
opt = Nio.options() 
nc  = Nio.open_file(ncFile,'r',opt)
wks = Ngl.open_wks('pdf','cnrm_1980')

##
## Set color map
##
res = Ngl.Resources()
res.wkColorMap = 'BlueWhiteOrangeRed'
Ngl.set_values(wks,res)



def set_rts_labels():
   """
   """
   rts_label = ['LH [kJ/kg]',\
                'DSE [kJ/kg]',\
                'MSE [kJ/kg]',\
                'Specific volume [m3/kg]',
                'P [hPa]',\
                'T [K]']
   
   flux_label = ['LH flux [PW]',\
                'DSE flux [PW]',\
                'MSE flux [PW]',\
                'Volume flux [m3/s]',
                'Pressure flux [hPa]',\
                'Sensible heat flux [PW]']
   
   labels = {}
   labels['rts_label'] = rts_label
   labels['flux_label'] = flux_label
   
   return labels


def read_data(ncFile,read_list=[],startYear=-1,endYear=-1,ix=0,iy=1,ltend=False):
   """
   Usage: 
      data = read_data(ncFile,startYear,endYear,ltend=False)
   """
   ##
   nc  = Nio.open_file(ncFile,'r')
   var_list = nc.variables.keys()
   if ( len(read_list) == 0 ):
      print ' Reading all variables '
      read_list = var_list
   else: 
      print ' Reading only selected variables '   
   
   print read_list
   
   print '---------------------------------------------------'
   print ' File: '+ncFile
   
   icoord  = nc.dimensions['coords'] 
   ntime   = nc.dimensions['time']
   nlev    = nc.dimensions['mod_lev']
   
   rts_lev = nc.variables['rts_lev'][:,:]
   rts_lev[0,:] = rts_lev[0,:] * 2.5 # sp. hum -> lat. heat
   rts_lev2 = nc.variables['rts_lev2'][:,:]
   rts_lev2[0,:] = rts_lev2[0,:] * 2.5 # sp. hum -> lat. heat
   
   vlon    = nc.variables['longitude'][:]
   vlat    = nc.variables['latitude'][:]
   
   if 'year' in var_list:
      yr      = nc.variables['year'][:]
   elif 'YR' in var_list:
      yr      = nc.variables['YR'][:]
   else:
      print ' I cant find any year array in the file '
      print ' Creating a generic year array (1,2,3... etc) '
      yr      = nc.arange(0,ntime)/12.
   
   if 'month' in var_list:
      mm      = nc.variables['month'][:]
   elif 'MM' in var_list:
      mm      = nc.variables['MM'][:]
   else:
      print ' I cant find any month array in the file '
      print ' Creating a generic month array (1,2,3... etc) '
      mm      = nc.arange(0,ntime)
   
   if (startYear == -1):
      start=0
   else:   
      start = np.min( np.where(yr == startYear) )
   if (endYear == -1):
      stop=-1
   else:
      stop  = np.max( np.where(yr == endYear) )
   
   print ' Start: year '+repr(yr[start])
   print ' End: year '+repr(yr[stop])
   
   print ' Read stream functions '
   data = {}
   data['years'] = yr
   data['months'] = mm
   data['lon'] = vlon
   data['lat'] = vlat
   data['rts_lev'] = rts_lev
   data['rts_lev2'] = rts_lev2
   data['lev'] = np.arange(0,nlev)
   
   if ('psirr' in nc.variables.keys() and 'psirr' in read_list):
      if (ltend):
         psrr = nc.variables['psirr'][start:stop+1,iy,ix,0,:,:] + \
                nc.variables['psirr2'][start:stop+1,iy,ix,0,:,:]
      else: 
         psrr = nc.variables['psirr'][start:stop+1,iy,ix,0,:,:]
      
      for jk in range(1,psrr.shape[2]):
         psrr[:,:,jk] = psrr[:,:,jk] + psrr[:,:,jk-1]
      
      psrr = psrr * 10**(-9)
      data['psirr'] = psrr
      
   if ('psiyr' in nc.variables.keys() and 'psiyr' in read_list):
      psir = nc.variables['psiyr'][start:stop+1,:,0,:,:]
      psir = psir * 10**(-9) 
      data['psiyr'] = psir
       
   if ('rtsyz' in nc.variables.keys() and 'rtsyz' in read_list):
      rtsyz = nc.variables['rtsyz'][start:stop+1,:,0,:,:]
      rtsyz[:,0,:,:] = 2.5 * rtsyz[:,0,:,:]
      data['rtsyz'] = rtsyz
      
   if ('rtsyr' in nc.variables.keys() and 'rtsyr' in read_list):
      rtsyr = nc.variables['rtsyr'][start:stop+1,:,:,:,:]
      data['rtsyr'] = rtsyr
      
   if ('rtsyr_zm' in nc.variables.keys() and 'rtsyr_zm' in read_list):
      rtsyr_zm = nc.variables['rtsyr_zm'][start:stop+1,:,:,:,:]
      data['rtsyr_zm'] = rtsyr_zm
   
   nc.close()
   
   return data

def cclap(temp,press,rh=1.,temp0=273.,e_s0=611.):
    """
    cclap - Module to calculate saturation specific humidity from
            the Clausius-Clapeyron equation.
    
    The equation: es = rh * e_s0 * exp( Lv / (Rv*temp0) - Lv / (Rv*temp) )
    
    
    Specific humidity: q = w / (w+1), w = epsilon*es/(press-es)
    
    Usage:
    
    qs, es = cclap(temp,press,rh=1.,temp0=273.,e_s0=611)
    
    In:
      temp     -  Temperature in Kelvins
      press    -  Pressure in Pascals.
      rh       -  Relative humidity (optional)
      temp0    -  Temperature in Kelvins with a measured es (optional)
      e_s0     -  Saturation vapor pressure at temp0 in Pascals (optional)
      
    Out:
      qs       -  Specific humidity
      es       -  Vapor pressure
    """
    
    import numpy as np
    
    ##
    ## Constants
    ##
    Rv  =  461.5      #Gas constant for water vapor
    Lv  =  2.5*10**6  #Latent heat in 1 kg of water vapor
    eps =  0.622
    
    es = rh * e_s0 * np.exp( Lv / (Rv*temp0) - Lv / (Rv*temp) )
    
    ws = eps*es/(press-es)
    
    qs = ws / (ws+1.)
    
    return qs,ws


def plot_cc(wks,cont,x='lh',y='dse'):
   """
   
   """
   ##
   if (x == 'lh' and y == 'dse'):
      vy = np.linspace(240,360,101)
      press = 101300.
      vx,es = cclap(vy*1000/1004.,press)
      vx = 2.5 * 10**3 * vx
      ccstring = 'C-C, RH=100%'
      ## Line of constant moist static energy
      vy2 = 345. - vx
      lres = Ngl.Resources()
      lres.gsLineDashPattern         =  1  
      lres.gsLineLabelFontHeightF    =  0.009
      lres.gsLineLabelString         =  'MSE = 345 kJ/kg'
      lres.gsLineDashSegLenF         =  0.09
      line = Ngl.add_polyline(wks,cont,vx,vy2,lres)
      
   if (x == 'lh' and y == 'mse'):
      vy = np.linspace(240,360,101)
      press = 101300.
      vx,es = 2.5 * 10**3 * cclap(vy*1000/1004.,press)
      vy = vy + vx
      ccstring = 'C-C, RH=100%'
      
      ## Line of constant moist static energy
      vy2 = 345. * np.ones(vx.shape[0])
      lres = Ngl.Resources()
      lres.gsLineDashPattern         =  1  
      lres.gsLineLabelFontHeightF    =  0.009
      lres.gsLineLabelString         =  'MSE = 345 kJ/kg'
      lres.gsLineDashSegLenF         =  0.09
      line = Ngl.add_polyline(wks,cont,vx,vy2,lres)
   
   
   lres = Ngl.Resources()
   lres.gsLineDashPattern         =  15  
   lres.gsLineLabelFontHeightF    =  0.009
   lres.gsLineLabelString         =  ccstring
   line = Ngl.add_polyline(wks,cont,vx,vy,lres)
   
   return wks,cont
         

def plot_psirr(wks, zplot, vx, vy, xlab = 'x', ylab = 'y', vpx=0.2, vpy=0.8, \
               hei=0.4, cc='dse',lsmooth=False, frame=True):
   """
   """
   
   
   ##
   ## Colors for hydrothermal stream functions
   ##
   ld   = 450.
   dl   = 60
   
   res = Ngl.Resources()
   res.wkColorMap = 'WhiteBlue'
   #res.wkColorMap = 'BlueWhiteOrangeRed'
   Ngl.set_values(wks,res)
   
   ##
   ## Reverse colormap
   ##
   del res
   res = Ngl.Resources()
   cmap = Ngl.retrieve_colormap(wks)
   cmap[2:,:] = cmap[-1:1:-1,:]
   res.wkColorMap = cmap
   Ngl.set_values(wks,res)
         
   psi_levels = np.arange(-ld,0.,dl)
   
   ##
   ## Positions of four subplots
   ##
   
   print '  Plotting stream functions in (r,r) coordinates'
   
   if(1):
      
      if(1):
         
         #vu = np.zeros(zplot[-1,:,:].shape)
         #vv = np.zeros(zplot[-1,:,:].shape)
         
         #vu[1:,:] = zplot[-1,0:-1,:] - zplot[-1,1:,:]
         #vv[:,1:] = zplot[-1,:,1:] - zplot[-1,:,0:-1]
               
         title = 'psi('+xlab+','+ylab+')'
                  
         print ' Amplitude, max. err : '+repr(np.max(np.abs(zplot[:,:]))),\
                                         repr(np.max(np.abs(zplot[:,-1])))
         
         del res
         res = Ngl.Resources()
         
         res.nglMaximize                      =  False
         res.nglPaperOrientation              =  'Portrait'
         res.nglFrame                         =  False
         res.nglDraw                          =  False
         
         res.cnFillOn                         =  True
         res.cnLinesOn                        =  True
         res.cnLineLabelsOn                   =  False
         res.cnLineLabelDensityF              =  2
         res.cnLineLabelBackgroundColor       =  -1
         res.cnLevelSelectionMode             =  'ExplicitLevels'
         res.cnLevels                         =  psi_levels
      
         res.pmLabelBarSide                   =  'Right'
         res.lbLabelBarOn                     =  True
         res.lbOrientation                    =  'Vertical'
         res.pmLabelBarDisplayMode            =  'Always'
         res.pmLabelBarWidthF                 =  0.03
         res.pmLabelBarHeightF                =  hei
         res.lbTitleString                    =  'Sv (10~S~9~N~kg/s)'
         res.lbLabelFontHeightF               =  0.01
         res.lbTitleFontHeightF               =  0.01
         
         res.sfXArray                         =  vx
         res.sfYArray                         =  vy
         
         res.tiXAxisString                    =  xlab
         res.tiYAxisString                    =  ylab
         
         res.vpWidthF                         =  hei
         res.vpHeightF                        =  hei
         
         res.vpYF                             =  vpy
         res.vpXF                             =  vpx
         
         if (lsmooth):
            sigma = 2
            order = 0
            for ji in range(0,2):
               ndimage.filters.gaussian_filter(zplot[:,:],sigma,order=order,\
                                               output=zplot[:,:],\
                                               mode='reflect', cval=0.0)
         
         cont = Ngl.contour(wks,zplot[:,:],res)
         
         wks,cont = plot_cc(wks,cont,x='lh',y=cc)
         Ngl.draw(cont)
         
         if(frame):
            Ngl.frame(wks)


def fill_blanks(ylab,nn):
   """
   """
   ## If not a list, make it a one-element list
   if ( isinstance(ylab,list) ):
      tmp = len(ylab)
      ylab2 = ylab
   else:
      tmp = 1
      ylab2 = [ylab]
   
   ## If not enough elements
   if (tmp != nn and nn > 1):
      print ' Not enough list elements. Filling in the blanks ...'
      ylab3 = []
      for n in range(0,nn):
         ylab3.append(ylab2[0])
      
      ylab2 = ylab3
   
   return ylab2



def plot_psiyr(wks,data,vx,vy,xlab='Latitude [deg N]',ylab = 'y',\
               vsurf=[],hei=0.1,ylon=True,xon=False,lbar=False,flipy=False,\
               vpx=0.15,vpy=0.8,lsmooth=False):
   """
   Plots meridional overturning stream functions
   
   Usage: 
      plot_psiyr(wks,data,vx,vy,xlab='Latitude',ylab='y',vsurf=np.array([]),\
                 hei=0.15,ylon=True,xon=True,lbar=True,flipy=False)
   
   Input:
      wks    -  Workstation ID
     data    -  2D matrix. Can be several matrices in a list
       vx    -  Latitude vector
       vy    -  y-vector. Can be several in a list
       
   Optional:
     xlab    -  Title for x-axis
     ylab    -  Title for y-axis
    vsurf    -  Vector of surface values to draw in
      hei    -  Height of plot(s)
     ylon    -  Turn on y labels on left side (else right)
      xon    -  Turn on x labels on bottom
     lbar    -  Turn on colorbar
   
   Output:
      None
      
   """
   ##
   res = Ngl.Resources()
   res.wkColorMap = 'ViBlGrWhYeOrRe'
   Ngl.set_values(wks,res)
   
   print '  Plotting overturning stream functions in (y,r) coordinates'
   print lbar
   if ( isinstance(data,list) ):
      nplots = len(data)
   else:
      nplots = 1
      data = [data]
   
   vy = fill_blanks(vy,nplots)
   ylab = fill_blanks(ylab,nplots)
   vsurf = fill_blanks(vsurf,nplots)
   ylon = fill_blanks(ylon,nplots)
   xon = fill_blanks(xon,nplots)
   lbar = fill_blanks(lbar,nplots)
   lbar[0] = True
   flipy = fill_blanks(flipy,nplots)
   
   hei = max([0.8/(nplots+1),hei])
   wth = min([0.7,hei * 4])
   vpy = 0.9 - hei * np.arange(0,nplots+1)
   
   ##
   ## Color levels for y-r stream functions
   ##
   colour_levels = np.arange(-190,210,20)
   
   for jn in range(0,nplots):
      
      zplot = data[jn]
      
      res = Ngl.Resources()
      
      res.nglMaximize                      =  False
      res.nglPaperOrientation              =  'Portrait'
      res.nglFrame                         =  False
      res.nglDraw                          =  False
      
      res.cnFillOn                         =  True
      res.cnLinesOn                        =  True
      res.cnLineLabelsOn                   =  False
      res.cnLineLabelDensityF              =  2
      res.cnLineLabelBackgroundColor       =  -1
      res.cnLevelSelectionMode             =  'ExplicitLevels'
      res.cnLevels                         =  colour_levels
      
      res.sfXArray                         =  vx
      res.sfYArray                         =  vy[jn]
      res.vpXF                             =  vpx
      res.vpWidthF                         =  wth
      res.vpYF                             =  vpy[jn]
      res.vpHeightF                        =  hei
      
      ##
      ## Axis
      ##
      res.tiYAxisFontHeightF               =  0.012
      res.tiYAxisString                    =  ylab[jn]
      res.tiXAxisString                    =  xlab
      
      if( ylon[jn] == True ):
         res.tiYAxisSide                   =  'Left'
         res.tmYLLabelsOn                  =  True
         res.tmYRLabelsOn                  =  False
         res.tmYLLabelFontHeightF          =  0.012
      else:
         res.tiYAxisSide                   =  'Right'
         res.tmYLLabelsOn                  =  False
         res.tmYRLabelsOn                  =  True
         res.tmYRLabelFontHeightF          =  0.012
      
      if( xon[jn] == True ):
         res.tiXAxisOn                     =  True
         res.tmXBLabelsOn                  =  True
         res.tiXAxisFontHeightF            =  0.012
         res.tmXBLabelFontHeightF          =  0.012
      else:
         res.tiXAxisOn                     =  False
         res.tmXBLabelsOn                  =  False
      
      ##
      ## Labelbar
      ##
      if( lbar[jn] == True ):
         res.lbTitleString         =  'Sverdrup (10~S~9~N~kg/s)'
         res.lbLabelFontHeightF    =  0.012
         res.lbTitleFontHeightF    =  0.015
         res.lbLabelBarOn          =  True
         res.lbOrientation         =  'Horizontal'
         
         res.pmLabelBarSide        =  'Top'
         res.pmLabelBarDisplayMode =  'Always'
         res.pmLabelBarWidthF      =  0.65
         res.pmLabelBarHeightF     =  0.07
      else:
         res.lbLabelBarOn            =  False
      
      ##
      ## Flip y axis?
      ##
      if( flipy[jn] == True ):
         res.trYReverse                    =  True
      else:
         res.trYReverse                    =  False
      
      
      ## Smooth
      if (lsmooth==True):
         sigma = 2
         order = 0
         for ji in range(0,5):
            ndimage.filters.gaussian_filter(zplot[:,:],sigma,order=order,\
                                            output=zplot[:,:],\
                                            mode='reflect', cval=0.0)
      
      cont = Ngl.contour(wks,zplot,res)
      
      ##
      ## Also draw the surface values
      ##
      if (vsurf[jn].shape[0] == vx.shape[0]):
         lres = Ngl.Resources()
         lres.gsLineThicknessF = 1.5
         lres.gsLineColor = 'black'
         lres.gsLineDashPattern = 1
         line = Ngl.add_polyline(wks,cont,vx,vsurf[jn],lres)
      else:
         print ' Length of vx and vsurf do not agree ' 
      
      ## And show the max/min values
      zmax = np.max( zplot.flatten() )
      zmin = np.min( zplot.flatten() )
      max_point = np.where(zplot == zmax)
      min_point = np.where(zplot == zmin)
      
      lmax = vx[max_point[1]][0]
      ymax = vy[jn][max_point[0]][0]
      print ' psimax = '+repr(zmax)
      print ' at lat = '+repr(lmax)+', y = '+repr(ymax)
      
      lmin = vx[min_point[1]][0]
      ymin = vy[jn][min_point[0]][0]
      print ' psimin = '+repr(zmin)
      print ' at lat = '+repr(lmin)+', y = '+repr(ymin)
      
      pres = Ngl.Resources()
      pres.gsMarkerColor = 'black'
      pres.gsMarkerThicknessF = 3.
      mark = Ngl.add_polymarker(wks,cont,lmax,ymax,pres)
      mark = Ngl.add_polymarker(wks,cont,lmin,ymin,pres)
      
      Ngl.draw(cont)
   
   #Ngl.frame(wks)
   
   flux_label = ['LH flux [PW]',\
                'DSE flux [PW]',\
                'MSE flux [PW]',\
                'Volume flux [m3/s]',
                'Pressure flux [hPa]',\
                'Sensible heat flux [PW]']
   
   jj = 0
   if(1):
      
      vy2 = np.zeros((nplots,vx.shape[0]))
      for jc in range(0,nplots):
         zplot = data[jc][:,:]
         chi = 0.5 * (vy[jc][1:] + vy[jc][:-1])
         for jj in range(0,vx.shape[0]):
            dpsi = zplot[1:,jj] - zplot[:-1,jj]
            vy2[jc,jj] = np.sum(dpsi * chi)
      
      vy2 = 10**(-3) * (-1) * vy2[:,:] # W->PW and change sign
      
      linecolors = ['blue','red','black','green','yellow','purple']
      
      del res
      res = Ngl.Resources()
      
      res.nglMaximize                      =  False
      res.nglPaperOrientation              =  'Portrait'
      res.nglFrame                         =  False
      res.nglDraw                          =  True
      
      res.xyLineColors                     =  linecolors
      #res.xyLineThicknesses                = [2,2,2]
      #res.xyLabelMode                      = 'Custom'
      #res.xyExplicitLabels                 = ['LH','DSE','MSE']
      #
      #if( np.mod(jc,2) == 0):
      #   res.tiYAxisSide                   =  'Left'
      #   res.tmYLLabelsOn                  =  True
      #   res.tmYRLabelsOn                  =  False
      #else:
      #   res.tiYAxisSide                   =  'Right'
      #   res.tmYLLabelsOn                  =  False
      #   res.tmYRLabelsOn                  =  True
      # 
      #if( jc == icoord-1):
      #   res.tiXAxisOn                     =  True
      #   res.tmXBLabelsOn                  =  True
      #else:
      #   res.tiXAxisOn                     =  False
      #   res.tmXBLabelsOn                  =  False
      
      ## Axis
      res.tiXAxisString                    =  'Latitude'
      res.tiYAxisString                    =  '[PW]'
      res.tiXAxisFontHeightF               =  0.012
      res.tiYAxisFontHeightF               =  0.012
      res.tmXBLabelFontHeightF             =  0.012
      res.tmYLLabelFontHeightF             =  0.012            
      res.trYMaxF                          =  6.
      res.trYMinF                          =  -6.
      
      res.vpXF                             =  vpx
      res.vpWidthF                         =  wth
      res.vpYF                             =  vpy[-1]
      res.vpHeightF                        =  hei
      
      cont = Ngl.xy(wks,vx,vy2[0:3],res)
      
      
   Ngl.frame(wks)


def make_tropical_average(vlat,data2d):
   """
   """
   
   j1 = np.min( np.where(vlat >= 15.)[0] )
   j2 = np.max( np.where(vlat <= -15.)[0] )
   
   data1d = np.mean( data2d[:,j2:j1+1], axis=1 )
   
   return data1d

def plot_xy(wks,vx,data_vy,vpx=0.2,vpy=0.8,hei=0.3,xlab='',ylab='',\
                  prof=False,yrev=True,\
                  linecolors = ['red','blue','black','green','yellow'],\
                  linethicknesses=[1,1,1,1,1],
                  linepatterns=[0,0,0,0,0],\
                  titles=[''],\
                  xmin=0.,xmax=1000.,ymin=0.,ymax=100.,\
                  frame=True):
   """
   """
   
   if ( isinstance(data_vy,list) ):
      nprof = len(data_vy)
   else:
      data_vy = [data_vy]
      nprof = 1
   
   vx = fill_blanks(vx,nprof)
   linecolors = fill_blanks(linecolors,nprof)
   linepatterns = fill_blanks(linepatterns,nprof)
   linethicknesses = fill_blanks(linethicknesses,nprof)
   
   if (prof):
      tmp = vx
      vx = data_vy
      data_vy = tmp
   
   
   
   res = Ngl.Resources()
   res.nglFrame = False
   res.nglPaperOrientation = 'Portrait' 
   res.nglMaximize = False
   res.nglDraw = False
   
   res.xyLineColors = linecolors
   
   res.tiYAxisString = ylab
   res.tiXAxisString = xlab
   
   res.trYReverse = yrev
   res.trXMinF = xmin
   res.trXMaxF = xmax
   res.trYMinF = ymin
   res.trYMaxF = ymax
   
   if(prof):
      wth = hei * 1./np.sqrt(2.)
   else:
      wth = hei * np.sqrt(2.)
   res.vpXF = vpx
   res.vpYF = vpy
   res.vpHeightF = hei
   res.vpWidthF = wth
   
   xy = Ngl.xy(wks,vx[0],data_vy[0],res)
   
   for jn in range(1,nprof):
      res = Ngl.Resources()
      res.gsLineColor                           =  linecolors[jn]
      res.gsLineDashPattern                     =  linepatterns[jn]
      res.gsLineThicknessF                      =  linethicknesses[jn]
      
      line = Ngl.add_polyline(wks,xy,vx[jn],data_vy[jn],res)
       
   Ngl.draw(xy)
   
   wth2 = 0.15
   hei2 = 0.08
   for jn in range(0,2):
      if (jn == 0):
         x0 = vpx+wth-(wth2)
         y0 = vpy+hei2
         st = 0
         sp = 4
      elif (jn == 1):
         x0 = vpx+wth-(2*wth2+0.05)
         y0 = vpy+hei2
         st = 4
         sp = 7
      lres = Ngl.Resources()
      lres.vpWidthF            =  wth2
      lres.vpHeightF           =  hei2
      lres.lgLineColors        =  linecolors[st:sp]
      lres.lgLineThicknesses   =  linethicknesses[st:sp] 
      lres.lgDashIndexes       =  linepatterns[st:sp]
      lres.lgLineLabelsOn      =  False
      lres.lgLabelFontHeightF  =  0.01
      lg = Ngl.legend_ndc(wks,sp-st,titles[st:sp],x0,y0,lres)
   #tres = Ngl.Resources()
   #tres.txFontHeightF = 0.015
   #tres.txJust = 'CenterLeft'
   #y = vpy + 0.03
   #x0 = vpx
   #dx = 0.15
   #print nprof
   #print titles
   #for jn in range(0,nprof):
   #   tres.txFontColor = linecolors[jn]
   #   x = x0 + dx*jn
   #   Ngl.text_ndc(wks,titles[jn],x,y,tres)
   
   if (frame):
      Ngl.frame(wks)


def plot_dots(wks,vx,vy,vpx=0.2,vpy=0.8,hei=0.3,xlab='',ylab='',\
              prof=False,lreg=True,\
              dotcolors = ['red','blue','black','green','yellow'],\
              dotsizes=[1,1,1,1,1],
              dotindices=[9,6,2,3,4,5],\
              titles=[''],\
              xmin=0.,xmax=1000.,ymin=0.,ymax=100.,\
              frame=True):
   """
   """
   
   ndots = len(vy)
   
   print dotcolors
   print dotindices
   print dotsizes
   
   dotcolors = fill_blanks(dotcolors,ndots)
   dotindices = fill_blanks(dotindices,ndots)
   dotsizes = fill_blanks(dotsizes,ndots)
   
   
   res = Ngl.Resources()
   res.nglFrame = False
   res.nglPaperOrientation = 'Portrait' 
   res.nglMaximize = False
   res.nglDraw = False
   
   res.xyMarkLineMode = 'Markers'
   res.xyMonoMarkerColor = False
   res.xyMarkerColors = dotcolors
   res.xyMarkerThicknesses = 2
   res.xyMarkerSizes = np.array(dotsizes) * 0.01
   res.xyMarkers = dotindices
   
   res.tiYAxisString = ylab
   res.tiXAxisString = xlab
   
   res.trXMinF = xmin
   res.trXMaxF = xmax
   res.trYMinF = ymin
   res.trYMaxF = ymax
   
   if(prof):
      wth = hei * 1./np.sqrt(2.)
   else:
      wth = hei * np.sqrt(2.)
   res.vpXF = vpx
   res.vpYF = vpy
   res.vpHeightF = hei
   res.vpWidthF = wth
   
   print vx
   print vy
   print xmin,xmax,ymin,ymax
   xy = Ngl.xy(wks,np.transpose(np.vstack((vx,vx))),np.transpose(np.vstack((vy,vy))),res)
   
   if (lreg):
      ## Fit line
      k,m,cor,sig,err = stats.linregress(vx,vy)
      vx = np.linspace(xmin,xmax,10)
      vy = vx * k + m
      
      # Plot line
      lres = Ngl.Resources()
      lres.gsLineColor = 'black'
      lres.gsLineDashPattern = 1
      lres.gsLineThicknessF = 2
      line = Ngl.add_polyline(wks,xy,vx,vy,lres)
      
      # Print correlations
      text = 'R = '+repr(cor)
      if (sig*100. < 1.):
         text = text + ' (99% sign.)'
      elif (sig*100. < 5.):
         text = text + ' (95% sign.)'
      else:
         text = text + ' (not sign.)'
      tres = Ngl.Resources()
      tres.txFontHeightF = 0.01
      tres.txJust = 'BottomLeft'
      tres.txFontColor = 'black'
      Ngl.text_ndc(wks,text,vpx+wth,vpy-hei-0.03,tres)
      
   Ngl.draw(xy)
   
   #wth2 = 0.15
   #hei2 = 0.08
   #for jn in range(0,2):
   #   if (jn == 0):
   #      x0 = vpx+wth-(wth2)
   #      y0 = vpy+hei2
   #      st = 0
   #      sp = 4
   #   elif (jn == 1):
   #      x0 = vpx+wth-(2*wth2+0.05)
   #      y0 = vpy+hei2
   #      st = 4
   #      sp = 7
   #   lres = Ngl.Resources()
   #   lres.vpWidthF            =  wth2
   #   lres.vpHeightF           =  hei2
   #   lres.lgLineColors        =  linecolors[st:sp]
   #   lres.lgLineThicknesses   =  linethicknesses[st:sp] 
   #   lres.lgDashIndexes       =  linepatterns[st:sp]
   #   lres.lgLineLabelsOn      =  False
   #   lres.lgLabelFontHeightF  =  0.01
   #   lg = Ngl.legend_ndc(wks,sp-st,titles[st:sp],x0,y0,lres)
   if (frame):
      Ngl.frame(wks)



data = read_data(ncFile,ix=0,iy=1,ltend=True)


##
## Plot psi(r,r)
##
vx = data['rts_lev'][0,:]
vy = data['rts_lev'][1,:]
zplot = data['psirr'][0,:,:]
print np.min(zplot.flatten())
plot_psirr(wks,zplot,vx,vy,xlab='Latent heat [kJ/kg]',ylab='Dry static energy [kJ/kg]')

##------------------------------------------------------------------------------

##
## Overturning in various coordinates
##

zplot = data['psiyr']
vx = data['lat'][:]

zplot = []
vy = []
vsurf = []
vcoord = [0,1,2]
flipy = [True,False,False]
xon = [False,False,True]
ylab = ['LH [kJ/kg]','DSE [kJ/kg]','MSE [kJ/kg]']
data['psiyr'][0,0,:,:] = data['psiyr'][0,0,:,:]*(-1.)
#data['psiyr'][0,4,:,:] = data['psiyr'][0,4,:,:]*(-1.)
for jn in vcoord:
   vy.append(data['rts_lev'][jn,:])
   zplot.append(data['psiyr'][0,jn,:,:])
   vsurf.append(data['rtsyz'][0,jn,-1,:])

plot_psiyr(wks,zplot,vx,vy,vsurf=vsurf,flipy=flipy,xon=xon,ylab=ylab)

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