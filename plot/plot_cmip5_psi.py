##
##   Python script to plot stream functions from CMIP5 runs
##

import matplotlib.pyplot as plt

import numpy as np
import os, sys, Nio, Ngl
from scipy import stats

def read_data(ncFile,read_list=[],startYear=-1,endYear=-1,ix=0,iy=1,ltend=True):
   """
   Usage: 
      data = read_data(ncFile,startYear,endYear,ltend=False)
   """
   ##
   nc  = Nio.open_file(ncFile,'r')
   var_list = nc.variables.keys()
   if ( len(read_list) == 0 ):
      read_list = var_list   
   
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
    

def plot_psirr(vpx=0.2, vpy=0.8, hei=0.4, lbar=True, xlab=1,ylab=1):
   
   ##
   ## Colors for hydrothermal stream functions
   ##
   ld   = 550.
   dl   = 50
   
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
   
   if(1):
      
      if(1):
         
         
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
         
         if (lbar):
            res.pmLabelBarSide                   =  'Top'
            res.lbLabelBarOn                     =  True
            res.lbOrientation                    =  'Horizontal'
            res.pmLabelBarDisplayMode            =  'Always'
            res.pmLabelBarWidthF                 =  0.6
            res.pmLabelBarHeightF                =  0.03
            res.pmLabelBarOrthogonalPosF         =  0.07
            #res.pmLabelBarParallelPosF           =  1.1
            res.lbTitleString                    =  'Sv (10~S~9~N~kg/s)'
            res.lbLabelFontHeightF               =  0.012
            res.lbTitleFontHeightF               =  0.015
         else:
            res.pmLabelBarDisplayMode            =  'Never'
         
         res.sfXArray                         =  vx
         res.sfYArray                         =  vy
         
         res.trXMinF                          =  -3
         res.trXMaxF                          =  55
         res.trYMinF                          =  250
         res.trYMaxF                          =  360
         
         res.tiXAxisString                    =  'Latent heat [kJ/kg]'
         res.tiYAxisString                    =  'Dry static energy [kJ/kg]'
         res.tiMainString                     =  title
         res.tiMainOffsetYF                   =  -0.045
         res.tiMainFontHeightF                =  0.015
         
         if (xlab == 1):
            res.tiXAxisOn                     =  True
            res.tmXBLabelsOn                  =  True
            res.tiXAxisSide                   =  'Bottom'
         elif (xlab == -1):
            res.tiXAxisOn                      =  True
            res.tmXBLabelsOn                  =  False
            res.tmXTLabelsOn                  =  True
            res.tiXAxisSide                   =  'Top'
         elif (xlab == 0):
            res.tiXAxisOn                     =  False
            res.tmXBLabelsOn                  =  False
            res.tmXTLabelsOn                  =  False
            
         if (ylab == 1):
            res.tiYAxisOn                     =  True
            res.tmYLLabelsOn                  =  True
            res.tmYRLabelsOn                  =  False
            res.tiYAxisSide                   =  'Left'
         elif (ylab == -1):
            res.tiYAxisOn                     =  True
            res.tmYLLabelsOn                  =  False
            res.tmYRLabelsOn                  =  True
            res.tiYAxisSide                   =  'Right'
         elif (ylab == 0):
            res.tiYAxisOn                     =  False
            res.tmYLLabelsOn                  =  False
            res.tmYRLabelsOn                  =  False
            
         res.vpWidthF                         =  hei
         res.vpHeightF                        =  hei
         
         res.vpYF                             =  vpy
         res.vpXF                             =  vpx
         
         lsmooth = False
         if (lsmooth):
            sigma = 2
            order = 0
            for ji in range(0,2):
               ndimage.filters.gaussian_filter(zplot[:,:],sigma,order=order,\
                                               output=zplot[:,:],\
                                               mode='reflect', cval=0.0)
         
         cont = Ngl.contour(wks,zplot[:,:],res)
         
         ## Clausius-Clapeyron line
         vy_cc = np.linspace(240,360,101)
         press = 101300.
         vx_cc,es = cclap(vy_cc*1000/1004.,press)
         vx_cc = 2.5 * 10**3 * vx_cc
         #ccstring = 'C-C, RH=100%'
         lres = Ngl.Resources()
         lres.gsLineDashPattern         =  15  
         lres.gsLineLabelFontHeightF    =  0.009
         lres.gsLineThicknessF          =  2
         #lres.gsLineLabelString         =  ccstring
         line = Ngl.add_polyline(wks,cont,vx_cc,vy_cc,lres)
         
         ## MSE profile
         lres = Ngl.Resources()
         lres.gsLineDashPattern         =  1
         lres.gsLineThicknessF          =  3  
         line = Ngl.add_polyline(wks,cont,vlh,vdse,lres)
         
         Ngl.draw(cont)


def plot_psirr_lines(vpx=0.2, vpy=0.8, hei=0.4, xlab=1,ylab=1):
   
   psi_levels = np.array([-500,-50])
   
   if(1):
      
      if(1):
         
         res = Ngl.Resources()
         
         res.nglMaximize                      =  False
         res.nglPaperOrientation              =  'Portrait'
         res.nglFrame                         =  False
         res.nglDraw                          =  False
         
         res.cnFillOn                         =  False
         res.cnLinesOn                        =  True
         res.cnLineLabelsOn                   =  True
         res.cnLineLabelDensityF              =  2
         res.cnLineLabelBackgroundColor       =  -1
         res.cnLevelSelectionMode             =  'ExplicitLevels'
         res.cnLevels                         =  psi_levels
         res.cnMonoLineColor                  =  True
         res.cnLineColor                      =  'blue'
         res.cnLineThicknessF                 =  2
         res.cnInfoLabelOn                    =  False
         
         res.sfXArray                         =  vx
         res.sfYArray                         =  vy
         
         res.trXMinF                          =  -3
         res.trXMaxF                          =  55
         res.trYMinF                          =  250
         res.trYMaxF                          =  360
         
         res.tiXAxisString                    =  'Latent heat [kJ/kg]'
         res.tiYAxisString                    =  'Dry static energy [kJ/kg]'
         res.tiMainString                     =  title
         
         if (xlab == 1):
            res.tiXAxisOn                     =  True
            res.tmXBLabelsOn                  =  True
            res.tiXAxisSide                   =  'Bottom'
         elif (xlab == -1):
            res.tiXAxisOn                      =  True
            res.tmXBLabelsOn                  =  False
            res.tmXTLabelsOn                  =  True
            res.tiXAxisSide                   =  'Top'
         elif (xlab == 0):
            res.tiXAxisOn                     =  False
            res.tmXBLabelsOn                  =  False
            res.tmXTLabelsOn                  =  False
            
         if (ylab == 1):
            res.tiYAxisOn                     =  True
            res.tmYLLabelsOn                  =  True
            res.tmYRLabelsOn                  =  False
            res.tiYAxisSide                   =  'Left'
         elif (ylab == -1):
            res.tiYAxisOn                     =  True
            res.tmYLLabelsOn                  =  False
            res.tmYRLabelsOn                  =  True
            res.tiYAxisSide                   =  'Right'
         elif (ylab == 0):
            res.tiYAxisOn                     =  False
            res.tmYLLabelsOn                  =  False
            res.tmYRLabelsOn                  =  False
            
         res.vpWidthF                         =  hei
         res.vpHeightF                        =  hei
         
         res.vpYF                             =  vpy
         res.vpXF                             =  vpx
         
         lsmooth = False
         if (lsmooth):
            sigma = 2
            order = 0
            for ji in range(0,2):
               ndimage.filters.gaussian_filter(zplot1[:,:],sigma,order=order,\
                                               output=zplot[:,:],\
                                               mode='reflect', cval=0.0)
         
         cont1 = Ngl.contour(wks,zplot1[:,:],res)
         
         ## Clausius-Clapeyron line
         vy_cc = np.linspace(240,360,101)
         press = 101300.
         vx_cc,es = cclap(vy_cc*1000/1004.,press)
         vx_cc = 2.5 * 10**3 * vx_cc
         #ccstring = 'C-C, RH=100%'
         lres = Ngl.Resources()
         lres.gsLineDashPattern         =  15  
         lres.gsLineLabelFontHeightF    =  0.009
         lres.gsLineThicknessF          =  3
         #lres.gsLineLabelString         =  ccstring
         line = Ngl.add_polyline(wks,cont1,vx_cc,vy_cc,lres)
         
         ## MSE profiles
         lres = Ngl.Resources()
         lres.gsLineDashPattern         =  1
         lres.gsLineThicknessF          =  7  
         lres.gsLineColor               =  'blue'
         line = Ngl.add_polyline(wks,cont1,vlh1,vdse1,lres)
         
         
         res.cnLineColor = 'red'
         
         cont2 = Ngl.contour(wks,zplot2[:,:],res)
         
         ## MSE profiles
         lres = Ngl.Resources()
         lres.gsLineDashPattern         =  1
         lres.gsLineThicknessF          =  7  
         lres.gsLineColor               =  'red'
         line = Ngl.add_polyline(wks,cont2,vlh2,vdse2,lres)
         
         Ngl.overlay(cont2,cont1)
         Ngl.draw(cont2)
         

def mon2ann(mondata):
   """
   """
   
   nmon = mondata.shape[0]
   anndata = np.array([])
   n = 0
   tmp = 0
   for jn in range(0,nmon):
      tmp = tmp + mondata[jn]
      n = n + 1
      if (n == 12 or n == nmon):
         tmp = tmp / 12.
         anndata = np.append(anndata,tmp)
         tmp = 0
         n = 0
   
   return anndata
   

def plot_xy(vpx=0.2,vpy=0.8,hei=0.3,wth=0.4,xlab='',ylab='',legend=True, yrev=False, lab=''):
   
   nprof = len(vy)
   
   res = Ngl.Resources()
   res.nglFrame = False
   res.nglPaperOrientation = 'Portrait' 
   res.nglMaximize = False
   res.nglDraw = False
   
   res.xyLineColors = colors
   res.xyLineThicknesses = thicknesses
   
   res.tiYAxisString = ylab
   res.tiXAxisString = xlab
   res.tiMainString = lab
   res.tiMainFontHeightF = 0.015
   res.tiMainOffsetYF = -hei/4.7
   res.tiMainOffsetXF = -wth * 2./5.
   
   res.trXMinF = xmin
   res.trXMaxF = xmax
   res.trYMinF = ymin
   res.trYMaxF = ymax
   res.trYReverse = yrev
   
   res.vpXF = vpx
   res.vpYF = vpy
   res.vpHeightF = hei
   res.vpWidthF = wth
   
   xy = Ngl.xy(wks,vx[0],vy[0],res)
   
   for jn in range(1,nprof):
      res = Ngl.Resources()
      res.gsLineColor                           =  colors[jn]
      res.gsLineDashPattern                     =  patterns[jn]
      res.gsLineThicknessF                      =  thicknesses[jn]
      
      if ( isinstance(vx,list) and len(vx) >= len(vy) ):
         tmpx = vx[jn]
      elif (isinstance(vx,list) and len(vx)-1 < jn ):
         tmpx = vx[0]
      else:
         tmpx = vx
         
      line = Ngl.add_polyline(wks,xy,tmpx,vy[jn],res)
       
   Ngl.draw(xy)
   
   wth2 = 0.08
   hei2 = 0.08
   for jn in range(0,2):
      if (jn == 0):
         x0 = vpx+wth-(wth2)
         y0 = vpy+hei2
         st = 0
         sp = 4
      elif (jn == 1 and nprof > 4):
         x0 = vpx+wth-(2*wth2+0.05)
         y0 = vpy+hei2
         st = 4
         sp = nprof
      lres = Ngl.Resources()
      lres.vpWidthF            =  wth2
      lres.vpHeightF           =  hei2
      lres.lgLineColors        =  colors[st:sp]
      lres.lgLineThicknesses   =  thicknesses[st:sp] 
      lres.lgDashIndexes       =  patterns[st:sp]
      lres.lgLineLabelsOn      =  False
      lres.lgLabelFontHeightF  =  0.008
      if(legend):
         lg = Ngl.legend_ndc(wks,sp-st,titles1[st:sp],x0,y0,lres)
   


def plot_dots(vpx=0.2,vpy=0.8,hei=0.3,wth=0.4,xlab='',ylab='',legend=True, yrev=False, lab=''):
   
   
   
   res = Ngl.Resources()
   res.nglFrame = False
   res.nglPaperOrientation = 'Portrait' 
   res.nglMaximize = False
   res.nglDraw = False
   
   res.xyMarkLineMode = 'Markers'
   res.xyMonoMarkerColor = True
   res.xyMarkerColor = 'black'
   res.xyMarkerThicknessF = 2
   res.xyMarkerSizeF = 0.01
   
   res.tiYAxisString = ylab
   res.tiXAxisString = xlab
   res.tiMainString = lab
   res.tiMainFontHeightF = 0.015
   res.tiMainOffsetYF = -hei/4.7
   res.tiMainOffsetXF = -wth * 2./5.
   
   res.trXMinF = xmin
   res.trXMaxF = xmax
   res.trYMinF = ymin
   res.trYMaxF = ymax
   res.trYReverse = yrev
   
   res.vpXF = vpx
   res.vpYF = vpy
   res.vpHeightF = hei
   res.vpWidthF = wth
   
   xy = Ngl.xy(wks,np.transpose(np.vstack((vx,vx))),np.transpose(np.vstack((vy,vy))),res)
   
   if(1):
      ## Fit line
      k,m,cor,sig,err = stats.linregress(vx,vy)
      vx2 = np.linspace(xmin,xmax,10)
      vy2 = vx2 * k + m
      
      # Plot line
      lres = Ngl.Resources()
      lres.gsLineColor = 'black'
      lres.gsLineDashPattern = 1
      lres.gsLineThicknessF = 2
      line = Ngl.add_polyline(wks,xy,vx2,vy2,lres)
      
      # Print correlations
      text = 'k = %.2f  R = %4.2f' % (k,cor)
      if (sig*100. < 1.):
         text = text + ' (99% sign.)'
      elif (sig*100. < 5.):
         text = text + ' (95% sign.)'
      else:
         text = text + ' (not sign.)'
      tres = Ngl.Resources()
      tres.txFontHeightF = 0.012
      tres.txJust = 'CenterLeft'
      tres.txFontColor = 'black'
      Ngl.text_ndc(wks,text,vpx,vpy+0.02,tres)
       
   Ngl.draw(xy)
   
   


##
##---------------------------------------------------------------------------
##

ncFiles1 = ['/nobackup/vagn2/x_joakj/data/canesm/psi/canesm_historical_1980-1999',\
            '/nobackup/vagn2/x_joakj/data/ccsm/psi/ccsm_historical_1980-1999',\
            '/nobackup/vagn2/x_joakj/data/ifs/psi/ecearth_historical_1980-1999_new',\
            '/nobackup/vagn2/x_joakj/data/gfdl/psi/gfdl_historical_1980-1999',\
            '/nobackup/vagn2/x_joakj/data/ipsl/psi/ipsl_historical_1980-1999',\
            '/nobackup/vagn2/x_joakj/data/noresm/psi/noresm_historical_1980-1999',\
            '/nobackup/vagn2/x_joakj/data/ifs/psi/era_interim_1980-1999']

ncFiles2 = ['/nobackup/vagn2/x_joakj/data/canesm/psi/canesm_rcp85_2080-2099',\
            '/nobackup/vagn2/x_joakj/data/ccsm/psi/ccsm_rcp85_2080-2099',\
            '/nobackup/vagn2/x_joakj/data/ifs/psi/ecearth_rcp85_2080-2099',\
            '/nobackup/vagn2/x_joakj/data/gfdl/psi/gfdl_rcp85_2080-2099',\
            '/nobackup/vagn2/x_joakj/data/ipsl/psi/ipsl_rcp85_2080-2099',
            '/nobackup/vagn2/x_joakj/data/noresm/psi/noresm_rcp85_2080-2099']
           
titles1 = ['CanESM2','CCSM4','EC-Earth2.3','GFDL-CM3','IPSL-CM5A-LR','NorESM1-M','ERA-Interim']
titles2 = ['CanESM2','CCSM4','EC-Earth2.3','GFDL-CM3','IPSL-CM5A-LR','NorESM1-M']


colors      = [ 2,3,4,7,13,23,1 ]
thicknesses = [ 2,2,2,2,2,2,5 ]
patterns    = [ 0,0,0,0,0,0,0 ]

##
## Open file and open PDF output
##
wks = Ngl.open_wks('pdf','cmip5')


hei = 0.22
x0 = 0.1 ; y0 = 0.9
dx = hei * 1.1 ; dy = hei * 1.1
vpx = [] ; vpy = []
xlab = [] ; ylab = []
for jn in range(0,len(ncFiles1)):
   if (jn == len(ncFiles1)-1):
      x = x0 + dx 
   else:
      x = x0 + np.mod(jn,3) * dx
      
   if (np.mod(jn,3) == 0):
      ylab.append(1)
   elif (np.mod(jn,3) == 2): 
      ylab.append(-1)
   else:
      ylab.append(0)
   
   if (jn == len(ncFiles1)-1 or jn == len(ncFiles1)-2 or jn == len(ncFiles1)-4):
      xlab.append(1)
   else:
      xlab.append(0)
   
   y = y0 - np.floor(jn/3) * dy
   vpx.append(x)
   vpy.append(y)


##
## Read data and calculate some stuff
##
amp_h = []
for jn in range(0,len(ncFiles1)):
   ##
   ## Read data
   ##
   ncFile = ncFiles1[jn]+'_av.nc'
   data = read_data(ncFile,ix=0,iy=1,ltend=True)
   
   vx = data['rts_lev'][0,:]
   vy = data['rts_lev'][1,:]
   zplot = data['psirr'][0,:,:]
   vlh = np.mean(data['rtsyz'][0,0,:,30:42],axis=1)
   vdse = np.mean(data['rtsyz'][0,1,:,30:42],axis=1)
   
   title = titles1[jn]
   amp_h.append( np.max(np.abs( zplot[:,:].flatten() )) )
   
   print title,amp_h[-1],repr(np.max(np.abs(zplot[:,-1])))
   
   
   ## Plot hydrothermal stream function
   if (jn == 1):
      plot_psirr(vpx=vpx[jn],vpy=vpy[jn],hei=hei,xlab=xlab[jn],ylab=ylab[jn],lbar=True)
   else:
      plot_psirr(vpx=vpx[jn],vpy=vpy[jn],hei=hei,xlab=xlab[jn],ylab=ylab[jn],lbar=False)
   

Ngl.frame(wks)

amp_r = []
xlab[4] = 1
for jn in range(0,len(ncFiles2)):
   ncFile = ncFiles2[jn]+'_av.nc' 
   data2 = read_data(ncFile,ix=0,iy=1,ltend=True)
      
   zplot = data2['psirr'][0,:,:]
   vlh = np.mean(data2['rtsyz'][0,0,:,30:42],axis=1)
   vdse = np.mean(data2['rtsyz'][0,1,:,30:42],axis=1)
   
   title = titles2[jn]
   amp_r.append( np.max(np.abs( zplot[:,:].flatten() )) )
   print title,amp_r[-1],repr(np.max(np.abs(zplot[:,-1])))
   
   ## Plot hydrothermal stream function
   if (jn == 1):
      plot_psirr(vpx=vpx[jn],vpy=vpy[jn],hei=hei,xlab=xlab[jn],ylab=ylab[jn],lbar=True)
   else:
      plot_psirr(vpx=vpx[jn],vpy=vpy[jn],hei=hei,xlab=xlab[jn],ylab=ylab[jn],lbar=False)


Ngl.frame(wks)

damp = np.array([amp_r])-np.array([amp_h[:-1]])
print ' Mean amp diff ',np.mean(damp)

##
## Colorful colormap
##
res = Ngl.Resources()
res.wkColorMap = 'hlu_default'
Ngl.set_values(wks,res)


##
## An overlay plot
##
ncFile = ncFiles1[0]+'_av.nc'
data = read_data(ncFile,ix=0,iy=1,ltend=True)
   
vx = data['rts_lev'][0,:]
vy = data['rts_lev'][1,:]
zplot1 = data['psirr'][0,:,:]
vlh1 = np.mean(data['rtsyz'][0,0,:,30:42],axis=1)
vdse1 = np.mean(data['rtsyz'][0,1,:,30:42],axis=1)

ncFile = ncFiles2[0]+'_av.nc'
data = read_data(ncFile,ix=0,iy=1,ltend=True)
zplot2 = data['psirr'][0,:,:]
vlh2 = np.mean(data['rtsyz'][0,0,:,30:42],axis=1)
vdse2 = np.mean(data['rtsyz'][0,1,:,30:42],axis=1)

title = titles1[0]
plot_psirr_lines(vpx=0.2,vpy=0.8,hei=0.3,xlab=1,ylab=1)

Ngl.frame(wks)


##
## Read in surface temperatures
##
tsFiles1 = ['/nobackup/vagn2/x_joakj/data/CanESM2_histr_tropical.txt_mon.txt',
            '/nobackup/vagn2/x_joakj/data/CCSM4_histr_tropical.txt_mon.txt',
            '/nobackup/vagn2/x_joakj/data/EC-Earth_histr_tropical.txt_mon.txt',
            '/nobackup/vagn2/x_joakj/data/GFDL-CM3_histr_tropical.txt_mon.txt',
            '/nobackup/vagn2/x_joakj/data/IPSL-CM5A-LR_histr_tropical.txt_mon.txt',
            '/nobackup/vagn2/x_joakj/data/NorESM1-M_histr_tropical.txt_mon.txt',
            '/nobackup/vagn2/x_joakj/data/ERA-Interim_histr_tropical.txt_mon.txt']
            
tsFiles2 = ['/nobackup/vagn2/x_joakj/data/CanESM2_rcp85_tropical.txt_mon.txt',
            '/nobackup/vagn2/x_joakj/data/CCSM4_rcp85_tropical.txt_mon.txt',
            '/nobackup/vagn2/x_joakj/data/EC-Earth_rcp85_tropical.txt_mon.txt',
            '/nobackup/vagn2/x_joakj/data/GFDL-CM3_rcp85_tropical.txt_mon.txt',
            '/nobackup/vagn2/x_joakj/data/IPSL-CM5A-LR_rcp85_tropical.txt_mon.txt',
            '/nobackup/vagn2/x_joakj/data/NorESM1-M_rcp85_tropical.txt_mon.txt']

ts_h = []
for jn in range(0,len(tsFiles1)):
   f = open(tsFiles1[jn])
   lines = f.readlines()
   ts_mean = 0
   n = 0
   for line in lines:
      year = eval(line[0:4])
      if (year >= 1980 and year <= 1999):
         ts_mean = ts_mean + eval(line[16:])
         n = n + 1
   ts_mean = ts_mean / float(n)
   ts_h.append(ts_mean)
   f.close()

ts_r = []
diff_ts = []
for jn in range(0,len(tsFiles2)):
   f = open(tsFiles2[jn])
   lines = f.readlines()
   ts_mean = 0
   n = 0
   for line in lines:
      year = eval(line[0:4])
      if (year >= 2080 and year <= 2099):
         ts_mean = ts_mean + eval(line[16:])
         n = n + 1
   ts_mean = ts_mean / float(n)
   ts_r.append(ts_mean)
   diff_ts.append( ts_mean - ts_h[jn] )
   f.close()

##
## Now plot tropical MSE profiles
##
mse_h = []
dse_h = []
lh_h = []
vp_h = []
for jn in range(0,len(ncFiles1)):
   ##
   ## Read data
   ##
   ncFile = ncFiles1[jn]+'_av.nc'
   data = read_data(ncFile,ix=0,iy=1,ltend=True)
   
   vlat = data['lat'][30:42] * np.pi/180.
   mse = data['rtsyr'][0,2,:,30:42,:]
   dse = data['rtsyr'][0,1,:,30:42,:]
   lh = data['rtsyr'][0,0,:,30:42,:] * 2.5
   
   mse = np.ma.mean(mse,axis=2)
   dse = np.ma.mean(dse,axis=2)
   lh  = np.ma.mean(lh,axis=2)
   
   for k in range(0,mse.shape[0]):
      mse[k,0] = np.ma.sum( mse[k,:] * np.cos(vlat) )/np.sum(np.cos(vlat)) 
      dse[k,0] = np.ma.sum( dse[k,:] * np.cos(vlat) )/np.sum(np.cos(vlat))
      lh[k,0]  = np.ma.sum( lh[k,:]  * np.cos(vlat) )/np.sum(np.cos(vlat))
   mse = mse[:,0]
   dse = dse[:,0]
   lh = lh[:,0]
   
   p = np.linspace(0,1100,mse.shape[0])
   
   vp_h.append(p)
   mse_h.append(mse)
   dse_h.append(dse)
   lh_h.append(lh)

vp_r = []
mse_r = []   
dse_r = []
lh_r = []

diff_mse = []
diff_dse = []
diff_lh = []
diff_mse_m = []
diff_dse_m = []
diff_lh_m = []
for jn in range(0,len(ncFiles2)):
   ##
   ## Read data
   ##
   ncFile = ncFiles2[jn]+'_av.nc'
   data = read_data(ncFile,ix=0,iy=1,ltend=True)
   
   vlat = data['lat'][30:42] * np.pi/180.
   mse = data['rtsyr'][0,2,:,30:42,:]
   dse = data['rtsyr'][0,1,:,30:42,:]
   lh = data['rtsyr'][0,0,:,30:42,:] * 2.5
   mse = np.ma.mean(mse,axis=2)
   dse = np.ma.mean(dse,axis=2)
   lh = np.ma.mean(lh,axis=2)
   for k in range(0,mse.shape[0]):
      mse[k,0] =  np.ma.sum( mse[k,:] * np.cos(vlat) )/np.sum(np.cos(vlat)) 
      dse[k,0] =  np.ma.sum( dse[k,:] * np.cos(vlat) )/np.sum(np.cos(vlat)) 
      lh[k,0] =  np.ma.sum( lh[k,:] * np.cos(vlat) )/np.sum(np.cos(vlat)) 
   mse = mse[:,0]
   dse = dse[:,0]
   lh = lh[:,0]
   p = np.linspace(0,1100,mse.shape[0])
   
   vp_r.append(p)
   mse_r.append(mse)
   dse_r.append(dse)
   lh_r.append(lh)
   
   diff_mse.append( (mse - mse_h[jn]) /mse_h[jn] * 100. )
   diff_dse.append( (dse - dse_h[jn]) /dse_h[jn] * 100. * dse_h[jn]/mse_h[jn] )
   diff_lh.append(  (lh  - lh_h[jn] ) /lh_h[jn]  * 100. * lh_h[jn] /mse_h[jn] )
   
   kmin = np.min( np.where( p >= 900 ) )
   kmax = np.max( np.where( p <= 1100 ) )
   print kmin,kmax
   diff_mse_m.append( np.ma.mean(diff_mse[-1][kmin:kmax+1]) )
   diff_dse_m.append( np.ma.mean(diff_dse[-1][kmin:kmax+1]) )
   diff_lh_m.append( np.ma.mean(diff_lh[-1][kmin:kmax+1]) )
   
   print titles2[jn],np.mean(diff_mse[-1][kmin:kmax+1])
   

vx = diff_lh
vy = vp_r
xmin = -5 ; xmax = 5.
ymin = 0. ; ymax = 1100.
plot_xy(vpx=0.2,vpy=0.9,hei=0.3,wth=0.2,yrev=True,legend=True,xlab='LH [kJ/kg]',lab='a)')

vx = diff_dse
vy = vp_r
xmin = -5 ; xmax = 5.
plot_xy(vpx=0.45,vpy=0.9,hei=0.3,wth=0.2,yrev=True,legend=False,xlab='DSE [kJ/kg]',lab='b)')

vx = diff_mse
vy = vp_r
xmin = -5 ; xmax = 5.
plot_xy(vpx=0.7,vpy=0.9,hei=0.3,wth=0.2,yrev=True,legend=False,xlab='MSE [kJ/kg]',lab='c)')

Ngl.frame(wks)

xmin = 0. ; xmax = 5.
ymin = 0. ; ymax = 5.
vx = np.array(diff_ts)
vy = np.array(diff_lh_m)
plot_dots(vpx=0.2,vpy=0.9,hei=0.2,wth=0.3,yrev=False,legend=True,xlab='~F33~D~F21~T',ylab='~F33~D~F21~LH',lab='a)')
vy = np.array(diff_dse_m)
plot_dots(vpx=0.58,vpy=0.9,hei=0.2,wth=0.3,yrev=False,legend=True,xlab='~F33~D~F21~T',ylab='~F33~D~F21~DSE',lab='b)')
vy = np.array(diff_mse_m)
plot_dots(vpx=0.2,vpy=0.6,hei=0.2,wth=0.3,yrev=False,legend=True,xlab='~F33~D~F21~T',ylab='~F33~D~F21~MSE',lab='c)')

Ngl.frame(wks)

##
## Meridional mass and energy fluxes
##
vlats = []

lh_mflux_h = []
lh_flux_h = []
lh_mflux_r = []
lh_flux_r = []

dse_mflux_h = []
dse_flux_h = []
dse_mflux_r = []
dse_flux_r = []

mse_mflux_h = []
mse_flux_h = []
mse_mflux_r = []
mse_flux_r = []

diff_lh_m = []
diff_dse_m = []
diff_mse_m = []
diff_lh_e = []
diff_dse_e = []
diff_mse_e = []

for jn in range(0,len(ncFiles1)):
   ncFile = ncFiles1[jn]+'_av.nc'
   data = read_data(ncFile,ix=0,iy=1,ltend=True)
   vlats.append(data['lat'][:])
   
   dl = data['rts_lev'][0,1]-data['rts_lev'][0,0]
   ds = data['rts_lev'][1,1]-data['rts_lev'][1,0]
   dm = data['rts_lev'][2,1]-data['rts_lev'][2,0]
   
   psi_l = data['psiyr'][0,0,:,:]
   psi_d = data['psiyr'][0,1,:,:]
   psi_m = data['psiyr'][0,2,:,:]
   
   lh_mflux_h.append( np.max(psi_l,axis=0) - np.min(psi_l,axis=0) )
   dse_mflux_h.append( np.max(psi_d,axis=0) - np.min(psi_d,axis=0) )
   mse_mflux_h.append( np.max(psi_m,axis=0) - np.min(psi_m,axis=0) )
   
   lh_flux_h.append( np.sum(psi_l * dl,axis=0) / 10**3 )
   dse_flux_h.append( np.sum(psi_d * ds,axis=0) / 10**3 )
   mse_flux_h.append( np.sum(psi_m * dm,axis=0) / 10**3 )
   

for jn in range(0,len(ncFiles2)):
   ncFile = ncFiles2[jn]+'_av.nc'
   data2 = read_data(ncFile,ix=0,iy=1,ltend=True)
   
   psi_l = data2['psiyr'][0,0,:,:]
   psi_d = data2['psiyr'][0,1,:,:]
   psi_m = data2['psiyr'][0,2,:,:]
   
   lh_mflux_r.append( np.max(psi_l,axis=0) - np.min(psi_l,axis=0) )
   dse_mflux_r.append( np.max(psi_d,axis=0) - np.min(psi_d,axis=0) )
   mse_mflux_r.append( np.max(psi_m,axis=0) - np.min(psi_m,axis=0) )
      
   lh_flux_r.append( np.sum(psi_l * dl,axis=0) / 10**3  )
   dse_flux_r.append( np.sum(psi_d * ds,axis=0) / 10**3  )
   mse_flux_r.append( np.sum(psi_m * dm,axis=0) / 10**3  )
      
   #diff_lh_m.append(  (lh_mflux_r[jn]  - lh_mflux_h[jn])  / lh_mflux_h[jn] *100.)
   #diff_dse_m.append( (dse_mflux_r[jn] - dse_mflux_h[jn]) / dse_mflux_h[jn]*100. )
   #diff_mse_m.append( (mse_mflux_r[jn] - mse_mflux_h[jn]) / mse_mflux_h[jn]*100. )
   diff_lh_m.append(  (lh_mflux_r[jn]  - lh_mflux_h[jn]) )
   diff_dse_m.append( (dse_mflux_r[jn] - dse_mflux_h[jn]))
   diff_mse_m.append( (mse_mflux_r[jn] - mse_mflux_h[jn]))
   
   diff_lh_e.append( lh_flux_r[jn] - lh_flux_h[jn] )
   diff_dse_e.append( dse_flux_r[jn] - dse_flux_h[jn] )
   diff_mse_e.append( mse_flux_r[jn] - mse_flux_h[jn] )


vx = vlats
vy = lh_mflux_h
xmin = -90 ; xmax = 90
ymin = 0 ; ymax = 200
plot_xy(vpx=0.15,vpy=0.9,hei=0.2,wth=0.3,ylab='Mass flux [Sv]',lab='a)')

vy = dse_mflux_h
plot_xy(vpx=0.15,vpy=0.65,hei=0.2,wth=0.3,ylab='Mass flux [Sv]',legend=False,lab='b)')

vy = mse_mflux_h
plot_xy(vpx=0.15,vpy=0.4,hei=0.2,wth=0.3,ylab='Mass flux [Sv]',legend=False,lab='c)',xlab='Latitude')

#vy = lh_mflux_r
vy = diff_lh_m
ymin = -50 ; ymax = 50
plot_xy(vpx=0.55,vpy=0.9,hei=0.2,wth=0.3,ylab='~F33~D~F21~ Mass flux [Sv]',legend=False,lab='d)')

#vy = dse_mflux_r
vy = diff_dse_m
plot_xy(vpx=0.55,vpy=0.65,hei=0.2,wth=0.3,ylab='~F33~D~F21~ Mass flux [Sv]',legend=False,lab='e)')

#vy = mse_mflux_r
vy = diff_mse_m
plot_xy(vpx=0.55,vpy=0.4,hei=0.2,wth=0.3,ylab='~F33~D~F21~ Mass flux [Sv]',legend=False,lab='f)',xlab='Latitude')


for jn in range(0,len(ncFiles2)):
   fig = plt.figure()
   ax = fig.add_subplot(111)
   ax.plot(vx[0],np.array(lh_mflux_h[jn]),'k-')
   ax.plot(vx[0],np.array(lh_mflux_r[jn]),'k--')
   ax.plot(vx[0],np.array(dse_mflux_h[jn]),'g-')
   ax.plot(vx[0],np.array(dse_mflux_r[jn]),'g--')
   ax.plot(vx[0],np.array(mse_mflux_h[jn]),'r-')
   ax.plot(vx[0],np.array(mse_mflux_r[jn]),'r--')
   
   #print titles1[jn]
   #print 'lh'
   #print np.max( lh_mflux_h[jn] ), np.max( lh_mflux_r[jn] )
   #print 'dse'
   #print np.max( dse_mflux_h[jn] ), np.max( dse_mflux_r[jn] )
   #print 'mse'
   #print np.max( mse_mflux_h[jn] ), np.max( mse_mflux_r[jn] )
   
Ngl.frame(wks)

#plt.show()


vx = vlats
vy = lh_flux_h
xmin = -90 ; xmax = 90
ymin = -6 ; ymax = 6
plot_xy(vpx=0.15,vpy=0.9,hei=0.2,wth=0.3,ylab='LH flux [PW]',lab='a)')

vy = dse_flux_h
ymin = -6 ; ymax = 6
plot_xy(vpx=0.15,vpy=0.65,hei=0.2,wth=0.3,ylab='DSE flux [PW]',legend=False,lab='b)')

vy = mse_flux_h
ymin = -6 ; ymax = 6
plot_xy(vpx=0.15,vpy=0.4,hei=0.2,wth=0.3,xlab='Latitude',ylab='MSE flux [PW]',legend=False,lab='c)')

#vy = lh_flux_r
vy = diff_lh_e
ymin = -2 ; ymax = 2
plot_xy(vpx=0.55,vpy=0.9,hei=0.2,wth=0.3,ylab='~F33~D~F21~ LH flux [PW]',legend=False,lab='d)')

#vy = dse_flux_r
vy = diff_dse_e
ymin = -2 ; ymax = 2
plot_xy(vpx=0.55,vpy=0.65,hei=0.2,wth=0.3,ylab='~F33~D~F21~ DSE flux [PW]',legend=False,lab='e)')

#vy = mse_flux_r
vy = diff_mse_e
ymin = -2 ; ymax = 2
plot_xy(vpx=0.55,vpy=0.4,hei=0.2,wth=0.3,xlab='Latitude',ylab='~F33~D~F21~ MSE flux [PW]',legend=False,lab='f)')


Ngl.frame(wks)


##
## Hydrothermal mass and energy fluxes
##
dse_mflux_h = []
lh_mflux_h = []

dse_flux_h = []
lh_flux_h = []
lh_flux1_h = []
lh_flux2_h = []

for jn in range(0,len(ncFiles1)):
   data = read_data(ncFiles1[jn]+'_av.nc',ix=0,iy=1,ltend=True)
   psirr = data['psirr'][0,:,:]
   
   psi0 = 50
   psi = np.where( psirr < -50., psirr, 0. )
   dpsi = psi[:,1:] - psi[:,0:-1]
   psi = 0.5 * (psi[:,1:] + psi[:,0:-1])
   vdse = data['rts_lev'][1,:]
   vlh  = (data['rts_lev'][0,1:]+data['rts_lev'][0,0:-1])/2.
   
   i = np.min( np.where( vlh >= 10. )[0] )
   j = np.min( np.where( vdse >= 320. )[0] )
   
   lh_mflux_h.append( psi[j,:] )
   dse_mflux_h.append( psi[:,i] )
   
   lh_flux = np.zeros((psi.shape[0]))
   for jj in range(0,psi.shape[0]):
      lh_flux[jj] = np.sum( dpsi[jj,:] * vlh / 10**(3) )
   
   dse_flux = np.zeros((psi.shape[1]))
   for jj in range(0,psi.shape[1]):
      dse_flux[jj] = np.sum( dpsi[:,jj] * vdse / 10**3 )
   
   lh_flux1  = []
   lh_flux2  = []
   for jj in range(0,psi.shape[0]):
      if (np.min(psi[jj,:]) < 0.):
         i = np.min( np.where( psi[jj,:] - np.min(psi[jj,:]) <= 0. )[0] )
         #print data['rts_lev'][0,i]
         lh_flux1.append( np.sum(dpsi[jj,i:] * vlh[i:]  / 10**3) )
         lh_flux2.append( np.sum(dpsi[jj,0:i]* vlh[0:i] / 10**3) )
         if (jj == j):
            tmp1 = np.sum( (dpsi[jj,i:]) * (vdse[1]-vdse[0]) )
            tmp2 = np.sum( (dpsi[jj,0:i]) * (vdse[1]-vdse[0]) )
            print tmp1,tmp2,i,titles1[jn]
      else:
         lh_flux1.append(0.)
         lh_flux2.append(0.)
   lh_flux1 = np.array(lh_flux1)
   lh_flux2 = np.array(lh_flux2)
   #print data['rts_lev'][1,j],lh_flux[j],lh_flux1[j],lh_flux2[j]
   lh_flux1_h.append(lh_flux1[j])
   lh_flux2_h.append(lh_flux2[j])
   
   lh_flux_h.append( lh_flux )
   dse_flux_h.append( dse_flux )
   

dse_mflux_r = []
lh_mflux_r = []

dse_flux_r = []
lh_flux_r = []
lh_flux1_r = []
lh_flux2_r = []

diff_lh_mflux = []
diff_dse_mflux = []

diff_lh_flux = []
diff_dse_flux = []

for jn in range(0,len(ncFiles2)):
   data = read_data(ncFiles2[jn]+'_av.nc',ix=0,iy=1,ltend=True)
   psirr = data['psirr'][0,:,:]
   
   psi0 = 50
   psi = np.where( psirr < -50., psirr, 0. )
   dpsi = psi[:,1:] - psi[:,0:-1]
   psi = 0.5 * (psi[:,1:] + psi[:,0:-1])
   vdse = data['rts_lev'][1,:]
   vlh  = (data['rts_lev'][0,1:]+data['rts_lev'][0,0:-1])/2.
   
   i = np.min( np.where( vlh >= 10. )[0] )
   j = np.min( np.where( vdse >= 320. )[0] )
   lh_mflux_r.append( psirr[j,:] )
   dse_mflux_r.append( psirr[:,i] )
   lh_flux = np.zeros((psi.shape[0]))
   for jj in range(0,psi.shape[0]):
      lh_flux[jj] = np.sum( dpsi[jj,:] * vlh / 10**(3) )
   dse_flux = np.zeros((psi.shape[1]))
   for jj in range(0,psi.shape[1]):
      dse_flux[jj] = np.sum( dpsi[:,jj] * vdse / 10**3 )
   
   lh_flux1  = []
   lh_flux2  = []
   for jj in range(0,psi.shape[0]):
      if (np.min(psi[jj,:]) < 0.):
         i = np.min( np.where( psi[jj,:] - np.min(psi[jj,:]) <= 0. )[0] )
         lh_flux1.append( np.sum(dpsi[jj,i:] * vlh[i:]  / 10**3) )
         lh_flux2.append( np.sum(dpsi[jj,0:i]* vlh[0:i] / 10**3) )
         if (jj == j):
            tmp1 = np.sum( (dpsi[jj,i:]) * (vdse[1]-vdse[0]) )
            tmp2 = np.sum( (dpsi[jj,0:i]) * (vdse[1]-vdse[0]) )
            print tmp1,tmp2,i,titles2[jn]
      else:
         lh_flux1.append(0.)
         lh_flux2.append(0.)
   lh_flux1 = np.array(lh_flux1)
   lh_flux2 = np.array(lh_flux2)
   lh_flux1_r.append(lh_flux1[j])
   lh_flux2_r.append(lh_flux2[j])
   
   dse_flux_r.append( dse_flux )
   lh_flux_r.append( lh_flux )
   
   #tmp = np.where(lh_flux_h[jn] != 0.,(lh_flux - lh_flux_h[jn])/lh_flux_h[jn] * 100.,0.)
   diff_lh_flux.append( lh_flux - lh_flux_h[jn] )

   
vx = [vdse]
vy = lh_flux_h

xmin = 240.
xmax = 360.
ymin = 0.
ymax = 15.
plot_xy(vpx=0.1,vpy=0.9,hei=0.2,wth=0.3,ylab='LH flux [PW]',lab='a)')

ymin = -5.
ymax = 5.
vy = diff_lh_flux
plot_xy(vpx=0.1,vpy=0.65,hei=0.2,wth=0.3,legend=False,xlab='DSE [kJ/kg]',ylab='~F33~D~F21~F',lab='b)')

#print dse_flux_h[0]

vx = [vlh]
vy = dse_flux_h
xmin = 0.
xmax = 50.
ymin = -25.
ymax = 0.
#plot_xy(vpx=0.5,vpy=0.9,hei=0.2,wth=0.3)

vy = dse_flux_r
#plot_xy(vpx=0.5,vpy=0.6,hei=0.2,wth=0.3)

Ngl.frame(wks)

Ngl.end()
