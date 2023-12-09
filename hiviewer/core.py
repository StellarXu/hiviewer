#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
core.py
Created on someday Jan, 2020
@author: Xu Chen, Jilin University / NAOC / UCAS.

This code is used to: Inspect and show the fits image,  Overlap different images (or contours), Process the data cube preliminarily.

"""
#import matplotlib as mpl
#mpl.use('Agg')

from astropy.io import fits
from matplotlib import pyplot as plt
import numpy as np
from astropy import units as u
from astropy import wcs, log
from spectral_cube import SpectralCube # 3D or 4D fits need
import os 

from .utils import line_set

#####################################################################################################
#####################################################################################################     
"""
%%% Using 'self.ndim' to check its shape.
%%% for (x,y,z),(x,y)is the 2D position and z is the third axis (usually freq or velocity).
%%% for (x,y,z,u),(x,y)is the 2D position and z is the third axis (usually freq or velocity),
u is the fourth (maybe stokes or something...)
"""
class FitsPic(object):
    """ 
    
    Attributes & Properties
    ----------
    filename: str   Name of FITS file
    ndim:  dimension number of the data
    ctype: ['coor--pro']    Coordinate type  and  projection.
    crval: Coordinate value at reference point (one per axis).
    crpix: Array location of the reference point in pixels (one per axis)
    cuint: Units for stated coordinate axes
    cdelt: Coordinate increment at reference point
    velo:  spectral axis of the cube
    fill_nans_data: fill Nans of the data with zero.
    beam:  resolution information of the telescope
    
    Methods
    ----------
    deg2pix:         ra,dec coordinates to pixel position 
    pix2deg:         pixel position to ra,dec coordinates
    velo2sel:        velocity(km/s) to selected slice number
    get/plot_slice:  get a 2D slice of the cube / show the  slice image.
    get/plot_spec:   get/plot the spectrum of a certain point.
    plot_side_slice: show a cube slice from one side.
    slab_moment:     compute 3 moment maps.
    plot_contour:    plot the contour on a certain slice
    add_beam:        add a radio beam.
    reproj_convolve_2D: astropy.convolve and reproject the 2D fits to target fits
    reproj_convolve_cube: SpectralCube.convolve and reproject the fits cube to target fits.
    center_levels :  Expirical contour levels for a 2D FitsPic central part.
    
    """
    def __init__(self, filename):
        self.filename = filename
        self.hdul = fits.open(filename)
        self.hd0 = self.hdul[0].header
        if self.hd0['NAXIS'] != 0:
            self.data = self.hdul[0].data #  DATA
        else:
            self.data = self.hdul[1].data
            self.hd0 = self.hdul[1].header
        self.ndim = len(self.data.shape)
        self.wcs_obj = wcs.WCS(self.hd0)
        self.wcs_cel = wcs.WCS(self.hd0).celestial
        # delete the z(or z,u) axis, and only reserve the celestial body.
        if (self.ndim == 3) or (self.ndim == 4) :
            self.cube = SpectralCube.read(self.hdul).with_spectral_unit(u.km / u.s) 
            #Initiate a SpectralCube.Use this to inspect the cube.
        self.hdul.close()
        
        
    @property    
    def ctype(self):
        return self.wcs_obj.wcs.ctype
    @property
    def crval(self):
        return self.wcs_obj.wcs.crval
    @property
    def crpix(self):
        return self.wcs_obj.wcs.crpix
    @property
    def cunit(self):
        return self.wcs_obj.wcs.cunit
    @property
    def cdelt(self):
        return self.wcs_obj.wcs.cdelt
    @property
    def velo(self):
        """
        spectral axis of the cube
        """ 
        if self.ndim == 2:
            log.warning("It doesn't have the velocity parameter!")
            v = " Nothing"
        else:
            v = self.cube.spectral_axis
        return v
    
      
    def fill_nans_data(self,fillnum = 0):
        """
        fill Nans of the data with a constant,default = 0.
        """
        from copy import deepcopy
        nan_loc = np.isnan(self.data)
        data_nonans = deepcopy(self.data)
        data_nonans[nan_loc] = fillnum
        return data_nonans
    
    
    def deg2pix(self,ra,dec,around=True):
        """
        ra,dec coordinates to pixel position 
        pixels start from (0,0)
        """
        px,py=self.wcs_cel.all_world2pix(ra,dec,0)
        if around:
            x,y=np.around([px,py]).astype(int)
            return x,y
        else:
            return px,py

    
    def pix2deg(self,px,py):
        """ 
        pixel position to ra,dec coordinates
        pixels start from (0,0)
        """
        ra,dec=self.wcs_cel.all_pix2world(px,py,0)
        return ra,dec

    def velo2sel(self,v):
        """ 
        velocity(km/s) to selected slice number,just like z in (z,y,x) of a cube
        """
        velo = self.velo.value

        if self.cunit[2] == 'm/s':
            vdelt=self.cdelt[2]/1e3
        else:
            vdelt=self.cdelt[2]
        
        vdelt = np.abs(vdelt)
        sel=np.where((velo>v-vdelt/2)&(velo<v+vdelt/2))[0][0]
        return sel
    
    def get_slice(self, sel = 0,axis = 0 ):
        """
        for (x,y,z),(x,y)is the 2D position and z is the third parameter 
        (usually freq or velocity).select which slice(z) you want.
        for dimension is 3: axis = 0 means you fixed the velocity, then get an image;
            axis = 1 or 2 means you fixed dec or ra
        """
        if self.ndim == 3:
            if axis == 0:
                s=self.data[sel,:,:]
            elif axis == 1:
                s=self.data[:,sel,:]
            elif axis == 2:
                s=self.data[:,:,sel]
        elif self.ndim == 4:
            #usually , its shape is (0,0,dec,ra)
            s=self.data[0,0,:,:]
        elif self.ndim == 2:
            s=self.data

        return s
     

    def plot_slice(self ,sel =0, cmap = 'gray', figsize=(6,5),picname = './',per_vmin_max = None,
                   coord_format = None, save = False, xylim = None, vmin_max = None, cbar_label = '$\mathrm{Jy \ beam^{-1}}$'):
        """
        select which slice you want.It will show a slice image.
        Or use #self.cube[sel,:,:].quicklook()#
        The syntax for the format string is the following:
            format	result
            'dd'	'15d'
            'dd:mm'	'15d24m'
            'dd:mm:ss'	'15d23m32s'
            'dd:mm:ss.s'	'15d23m32.0s'
            'dd:mm:ss.ssss'	'15d23m32.0316s'
            'hh'	'1h'
            'hh:mm'	'1h02m'
            'hh:mm:ss'	'1h01m34s'
            'hh:mm:ss.s'	'1h01m34.1s'
            'hh:mm:ss.ssss'	'1h01m34.1354s'
            'd'	'15'
            'd.d'	'15.4'
            'd.dd'	'15.39'
            'd.ddd'	'15.392'
            'm'	'924'
            'm.m'	'923.5'
            'm.mm'	'923.53'
            's'	'55412'
            's.s'	'55412.0'
            's.ss'	'55412.03'
            'x.xxxx'	'15.3922'
            '%.2f'	'15.39'
            '%.3f'	'15.392'
            '%d'	'15'

        """
        if self.ndim != 2:
            print(f"The 3rd parameter is {self.velo[sel]}")
        else:
            print('Showing image.')
        fig = plt.figure(figsize=figsize)
        ax=fig.add_subplot(111, projection=self.wcs_cel)
        #print(f'max={self.get_slice(sel).max()},min={self.get_slice(sel).min()}')
        if vmin_max != None:
            im=ax.imshow(self.get_slice(sel),vmin=vmin_max[0],vmax=vmin_max[1],origin='lower', cmap = cmap)
        else:
            if per_vmin_max == None:   
                im=ax.imshow(self.get_slice(sel),origin='lower', cmap = cmap )
            else:
                from .utils import percent_vminmax
                data = self.get_slice(sel)
                vmin,vmax = percent_vminmax(data,percent = per_vmin_max)
                im=ax.imshow(self.get_slice(sel),vmin=vmin,vmax=vmax,origin='lower', cmap = cmap)

   
        cbar = plt.colorbar(im,pad=.01)
        cbar.set_label(f'{cbar_label}', size=15)
        if xylim is not None:
            xlim = xylim[:2]; ylim = xylim[-2:]
        else:
            xlim = None; ylim = None

        line_set(ax, xlabel = 'RA (J2000)',ylabel = 'Dec (J2000)', xlim=xlim,ylim=ylim,)

        if coord_format is not None:
            ra=ax.coords['ra'];dec=ax.coords['dec'];ra.set_major_formatter(coord_format);dec.set_major_formatter(coord_format)
        plt.tight_layout()
        if save:
            plt.savefig(picname+'_slice.png',dpi=300,bbox_inches='tight')
        return ax 
        #self.cube[sel,:,:].quicklook()
    

    def plot_side_slice(self,sel =0, axis=1,cmap = 'gray',figsize=(6,5),per_vmin_max = None,
                        xylim=None,picname = './',save = False,vmin_max = None,velo = None,trans = False):
        """
        select which slice you want.It will show a slice image along ra axis(2) or dec axis(1)
        It looks like a waterfall figure at axis = 1 or 2 ,means you fixed dec or ra
        print an approximate fixed ra or dec.
        """

        if self.ndim != 3:
            raise(ValueError("It must be a 3D cube!")) 
        else:
            
            if axis == 1:
                dec = self.pix2deg(px=self.crpix[0],py=sel)[1]
                print(f"The fixed dec is nearly {dec} deg")
            elif axis==2:
                ra = self.pix2deg(py=self.crpix[1],px=sel)[0]
                print(f"The fixed ra is nearly {ra} deg")
            else:
                raise(ValueError("axis should be 1 or 2 !"))

        fig,ax = plt.subplots(figsize=figsize)
        #print(f'max={self.get_slice(sel).max()},min={self.get_slice(sel).min()}')

        if not isinstance(velo,np.ndarray):
            velo = self.velo.value
             
        data = self.get_slice(sel,axis)
        if trans:
            extent=(velo[0],velo[-1],0,data.shape[1])
            data = data.T
        else:
            extent=(0,data.shape[1],velo[0],velo[-1])

        if vmin_max != None:
            im=ax.imshow(data,origin='lower', cmap = cmap,aspect='auto',
                        vmin = vmin_max[0],vmax = vmin_max[1],extent = extent)
        else:
            if per_vmin_max == None:   
                im=ax.imshow(data,origin='lower', cmap = cmap,aspect='auto',extent = extent)

            else:
                from .utils import percent_vminmax
                vmin,vmax = percent_vminmax(data,percent = per_vmin_max)
                im=ax.imshow(data,origin='lower', cmap = cmap,aspect='auto',
                            vmin = vmin,vmax = vmax, extent = extent)
        
        if axis == 1:
            title = f"dec nearly {dec:.4f} deg"
        elif axis==2:
            title = f"ra nearly {ra:.4f} deg"
        
        if xylim is not None:
            xlim = xylim[:2]; ylim = xylim[-2:]
        else:
            xlim = None; ylim = None

        if trans == False:
            if axis==1:
                xlabel = 'RA [pixel]'
            else:
                xlabel = 'DEC [pixel]'
            ylabel = 'Velocity [km/s]' 
        else:
            if axis==1:
                ylabel = 'RA [pixel]'
            else:
                ylabel = 'DEC [pixel]'
            xlabel = 'Velocity [km/s]' 
                
        line_set(ax, xlabel = xlabel,ylabel = ylabel, xlim=xlim,ylim=ylim, title = title)        
        
        plt.tight_layout()
        plt.colorbar(im,pad=.01)
        if save:
            plt.savefig(picname+'_side.png',dpi=300,bbox_inches='tight')
        plt.show() 


    
    def get_spec(self,x,y,coord_type = 'radec'):
        """ 
        Only 3D
        Using this to get the spectrum of each point.
        x,y : coord_type = 'radec': ra,dec unit: deg
              coord_type = 'pixel': px,py
        """
        if coord_type == 'radec':
            px,py = self.deg2pix(ra = x,dec = y)
            return self.data[:,py,px] 
            #self.cube[:,py,px].quicklook()
        elif coord_type == 'pixel':
            return self.data[:,y,x]

    
    def plot_spec(self,x,y,coord_type = 'radec' ,picname = './',xlim = None, ylim = None,
                  velo = None, save = False):
        """ 
        Only 3D
        Using this to get the spectrum of each point.
        x,y : coord_type = 'radec': ra,dec unit: deg
              coord_type = 'pixel': px,py
        Or use #self.cube[:,py,px].quicklook()
        """
        if coord_type == 'radec':
            ra,dec = x,y
            px,py = self.deg2pix(ra,dec)
        elif coord_type == 'pixel':
            px,py = x,y
            ra,dec = self.pix2deg(px,py)
        spec = self.data[:,py,px]
        print(f'ra={ra},dec={dec}, px,py={px,py}')
        #log.info(' If you are uncertain to the unit of x axis, use self.cube[:,py,px].quicklook() to check .')
        #velo ,_,_= self.cube.world[:,py,px]#velo,dec,ra
        if not isinstance(velo,np.ndarray):
            velo = self.velo.value

        fig,ax=plt.subplots(figsize=(6,5))
        ax.plot(velo,spec)
        ax.grid()
        
        line_set(ax, xlabel = f"{self.ctype[2]}"+"$\mathrm{km \ s^{-1}}$", 
                 ylabel = 'intensity', xlim=xlim,ylim=ylim,)
        plt.tight_layout()
        if save:
            plt.savefig(picname+'_spec.png',dpi=300,bbox_inches='tight')
        plt.show() 
       
        #self.cube[:,y,x].quicklook()
        
        
    def slab_moment(self,v1,v2,filepath='./data/fast',
                    m0=True,m1=False,m2=False,m2_type='sigma',NAME_WITH_V = True):
        """ 
        Only for 3D cubes
        slab velocity from v1 km/s to v2 km/s.
            eg. For M33, v1 = -300km/s and v2 = -50km/s , like this
            
        You can choose to compute  3 moment maps. m0, m1 and m2.
        These new 2D maps will be saved as new FITS files.
        m2_type：
                'variance': M2
                'sigma'：sigma = np.sqrt(M2)
                'fwhm':fwhm = np.sqrt(8ln2)*sigma

        """
        cubeKMS_HI = self.cube.with_spectral_unit(u.km / u.s)
        slabHI = cubeKMS_HI.spectral_slab(v1*u.km/u.s, v2*u.km/u.s)

        if NAME_WITH_V:
            outname = filepath + '%+d'%v1 + '%+d'%v2
        else:
            outname = filepath

        # For M33, v1 = -300km/s and v2 = -50km/s , like this
        #slabHI = cubeKMS_HI.spectral_slab(-300*u.km/u.s, -50*u.km/u.s)
        log.info('Producing moment maps.'
              f'Slab velocity :{v1}km/s and {v2}km/s. ')
        if m0:
            m0_HI = slabHI.moment0()
            m0_HI.write(outname+'_HI-moment0.fits', overwrite=True)
            log.info(f" saved {outname}_HI-moment0.fits ")
        if m1:
            m1_HI = slabHI.moment1()
            m1_HI.write(outname+'_HI-moment1.fits', overwrite=True)
            log.info(f" saved {outname}_HI-moment1.fits ")
        if m2:
            if m2_type=='variance':
                m2_HI = slabHI.moment2()
            elif m2_type=='sigma':
                ## sigma = np.sqrt(M2)
                m2_HI = slabHI.linewidth_sigma() 
            elif m2_type=='fwhm':
                ## fwhm = np.sqrt(8ln2)*sigma
                m2_HI = slabHI.linewidth_fwhm() 
            m2_HI.write(f"{outname}_{m2_type}_HI-moment2.fits", overwrite=True)
            log.info(f" saved {outname}_{m2_type}_HI-moment2.fits ")
        
    
        # Convert Moment_0 to a Column Density assuming optically thin media
        #hi_column_density = m0_HI * 1.82 * 10**18 / (u.cm * u.cm) * u.s / u.K / u.km
        
    
    def center_levels(self,data_unit = 'K'):
        """
        Expirical contour levels for a 2D FitsPic central part.
        But I'm not sure about the Suitabe contour levels!!
        At present it = (the Max data - np.log(a,b,23))[::-1] emmmm
        Sometimes it's wrong because of the bad data. 
        """
        m = np.max(self.fill_nans_data(fillnum = 0))
        if data_unit == 'K':
            level = (m-np.logspace(0.1,2.9,23))[::-1]
        else:
            level = (m-np.logspace(0.1,1.9,23))[::-1]
        log.info("It's just a expirical level created automaticly!")
        return level
        
        
    def plot_contour(self,levels,sel=0 ,vmin_max = None, per_vmin_max = None,
                     alpha=1,clabel=False,figsize=(6,5),xylim=None,
                     cbar_label='label', picname ='./',save = False,cmap='gray',coord_format = None):
        """ 
        plot contour on your selected slice . Levels decides the contour number.
        alpha belongs to the backgrond image.
        """
        
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection=self.wcs_cel)
        
        if vmin_max != None:
            im=ax.imshow(self.get_slice(sel),vmin=vmin_max[0],vmax=vmin_max[1],origin='lower', cmap = cmap)
        else:
            if per_vmin_max == None:   
                im=ax.imshow(self.get_slice(sel),origin='lower', cmap = cmap )
            else:
                from .utils import percent_vminmax
                data = self.get_slice(sel)
                vmin,vmax = percent_vminmax(data,percent = per_vmin_max)
                im=ax.imshow(self.get_slice(sel),vmin=vmin,vmax=vmax,origin='lower', cmap = cmap)   
        
        
        cbar = plt.colorbar(im,pad=.01)
        cbar.set_label(f'{cbar_label}', size=15)

        #ax.grid(color = 'white', ls = 'dotted', lw = 2)
        [ax.coords[ci].display_minor_ticks(True) for ci in range(2)]
        bx=ax.contour(self.get_slice(sel),levels = levels,colors='orange')
        if xylim is not None:
            xlim = xylim[:2]; ylim = xylim[-2:]
        else:
            xlim = None; ylim = None
        line_set(ax, xlabel = 'RA (J2000)',ylabel = 'Dec (J2000)', xlim=xlim,ylim=ylim,)    
        
        if clabel:
            ax.clabel(bx,fmt='%.1f')
        
        if coord_format is not None:
        	ra=ax.coords['ra'];dec=ax.coords['dec'];ra.set_major_formatter(coord_format);dec.set_major_formatter(coord_format)
        plt.tight_layout()
        if save:
            plt.savefig(picname+'_contour.png',dpi=300,bbox_inches='tight') 
        
        return ax
   
     
    def add_beam(self,bmaj=3/60,bmin=3/60,bpa=0,newname=None):
        """
        If there's no beam is defined for this SpectralCube or the beam information,
        add major/minor/position angle(deg). 
        """
        hdl=fits.open(self.filename)
        hd=hdl[0].header
        #hd['CTYPE3']='VELOCITY'#spectralcube 不可识别的要先转化一下
        hd['BMAJ']=bmaj
        hd['BMIN']=bmin
        hd['BPA']=bpa
        if newname == None:
            newname = self.filename
        fits.writeto(newname,hdl[0].data,hd,overwrite=True)
        hdl.close()
        
    def write_column_density(self, z = 0, beam = [3, 3], newname=None):
        """
        Change Moment 0 to log NHI map
        """
        def flux2col_density(S, z, beam):
            beam = np.array(beam) * 60
            logN = np.log10(1.1e24 * (1 + z)**2 * S / (beam[0] * beam[1]))
            return logN
        
        hdl=fits.open(self.filename)
        hd=hdl[0].header
        hd['BUNIT'] = 'log cm-2'
        
        data = flux2col_density(self.data, z, beam)
        
        if newname == None:
            newname = self.filename.replace('moment0', 'column_density')
        fits.writeto(newname, data, hd, overwrite=True)
        hdl.close()


    @property
    def beam(self):
        """
        resolution information of the telescope
        """
        from radio_beam import Beam
        self_beam = Beam.from_fits_header(self.hd0)
        return self_beam

    def reproj_convolve_2D(self,reproj_target,beam_target,fitsname,method='montage'):
        """       
        We use package montage/kapteyn/reproject to reproject self to reproj_target shape. 
        In this process, the pixel size will be the same as reproject target.
        Then use  convolve_fft in astropy to convolve beam to target beam.

        Read the fits header to get beam infomation.If it doesn't has,there will be a warning. Just 
        add a beam on the header.


        Parameters
        ------------
        self,reproj_target,beam_target: FitsPic Class Object.
        method: montage/kapteyn/reproject. Choose one of these packages to reproject.

        fitsname:str.The result fits will be named as [fitsname+'_convolve_reproj_result.fits']
        REMEMBER the 'self' object must have 2 dimensions.
        reproj_target should be 2D FitsPics.
        beam_target can be 2D/3D/4D FitsPics.
        
        Result
        --------
        A 2D fits        """
        from .utils import reproject_fits
        from astropy.convolution import convolve, convolve_fft
        result = fitsname+'_convolve_reproj_result.fits'
        reprojFits = fitsname+'_reproj.fits'
        if self.ndim != 2:
            raise(ValueError("The self object must have 2 dimensions."))
        else:
            if os.path.exists(result):
                print(f"{result} already exsist.Delete it...")
                os.remove(result)
            if os.path.exists(reprojFits):
                print(f"{reprojFits} already exsist.Delete it...")
                os.remove(reprojFits)

        reproject_fits(self.filename,reproj_target,beam_target,result=reprojFits,method=method,slabm0 = True)

        beam_self = self.beam
        beam_targ = beam_target.beam

        if beam_self == beam_targ:
            log.warning("The target beam is identical to the current beam. "
                          "Skipping convolution.")
            log.info(f' Create a new fits file named {reprojFits}')
            return 
        else:
            beam_self_to_targ = beam_targ.deconvolve(beam_self)
        
        rj_target_hd = reproj_target.hd0
        pix_scale = rj_target_hd['CDELT2']*3600*u.arcsec
        gauss_kern = beam_self_to_targ.as_kernel(pix_scale)

        img = fits.getdata(reprojFits)
        imgSmooth = convolve_fft(img, gauss_kern, allow_huge=True, normalize_kernel=True)
        fits.writeto(result,imgSmooth,fits.getheader(reprojFits),overwrite=True)

        if os.path.exists(result):
            log.info(' Create a new fits file '+result)
        os.remove(reprojFits)


    def reproj_convolve_cube(self,reproj_target,beam_target,fitsname,method='montage',
                             v1=None,v2=None,slabm0=True):
        """       
        We use package montage/kapteyn/reproject to reproject self to reproj_target shape. 
        In this process, the pixel size will be the same as reproject target.
        Then use convolve_to in SpectralCube to convolve beam to the target beam. 

        Read the fits header to get beam infomation.If it doesn't has,there will be a warning. Just 
        add a beam on the header.


        Parameters
        ------------
        self,reproj_target,beam_target: FitsPic Class Object.
        method: montage/kapteyn/reproject. Choose one of these packages to reproject.
        v1,v2: slab moment 0 , from v1 to v2. the same as slab_moment if slabm0=True
        fitsname:str.The result fits will be named as [fitsname+'_convolve_reproj_result.fits']

        REMEMBER the 'self' object must have 3 or 4 dimensions.In another words, it is a cube!
        reproj_target should be 2D FitsPics.
        beam_target can be 2D/3D/4D FitsPics.
        """ 
        from .utils import reproject_fits
        
        result = fitsname+'_convolve_reproj_result.fits'  
        if self.ndim == 2:
            raise(ValueError("The self object must have 3 or 4 dimensions.In another words, it is a cube!"))
        else:
            if os.path.exists(result):
                print(f"{result} already exsist.Delete it...")
                os.remove(result)

        beam_self = self.beam
        beam_targ = beam_target.beam

        if beam_self == beam_targ:
            log.warning("The target beam is identical to the current beam. "
                          "Skipping convolution.")
            new_cube= self.cube
        else:
            new_cube = self.cube.convolve_to(beam_targ)
            log.warning("It seems that SpectralCube.convolve_to will change the original cube resolution,\
                        so you should redefine the original cube if needed.")

        if slabm0:
            if self.data.shape[-3] != 1: #real cube
                cubeKMS_HI = new_cube.with_spectral_unit(u.km / u.s)
                #v1=-300
                #v2=-50
                slabHI = cubeKMS_HI.spectral_slab(v1*u.km/u.s, v2*u.km/u.s)
                print('Producing convolved moment 0 from the cube.'
                      f'Slab velocity :{v1}km/s and {v2}km/s. ')

                m0_HI = slabHI.moment0()
                m0_name = fitsname+'_convolved-moment0.fits'
                m0_HI.write(m0_name, overwrite=True)
                name = m0_name

            else: #don't slab moment 
                print(f"Its data shape is {self.data.shape}.So don't slab moment 0. ")
                new_cube_name = fitsname+'_convolved_cube.fits'
                new_cube.write(new_cube_name,overwrite=True)
                name = new_cube_name
            log.info("Reprojecting the moment0 of the cube...")
        else:
            new_cube_name = fitsname+'_convolved_cube.fits'
            new_cube.write(new_cube_name,overwrite=True)
            name = new_cube_name
            log.info("Reprojecting the whole cube...")

        reproject_fits(name,reproj_target,beam_target,result,method,slabm0)  

        if os.path.exists(fitsname+'_convolved-moment0.fits'):
            os.remove(fitsname+'_convolved-moment0.fits') 
        if os.path.exists(fitsname+'_convolved_cube.fits'):
            os.remove(fitsname+'_convolved_cube.fits')
        if os.path.exists(result):
            log.info(' Create a new fits file '+result)


