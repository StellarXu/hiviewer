#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
utils.py
Created on someday Jan, 2020
@author: Xu Chen, Jilin University / NAOC / UCAS.

"""
#import matplotlib as mpl
#mpl.use('Agg')

from astropy.io import fits
from matplotlib import pyplot as plt
import numpy as np
from astropy import units as u
from astropy import wcs, log
import os

def percent_vminmax(data,percent = None):
    if percent != None:
        vmin = np.nanpercentile(data,q = (1-percent)/2* 100,interpolation='nearest')
        vmax = np.nanpercentile(data,q = (1+percent)/2* 100,interpolation='nearest')
        return vmin,vmax
    else:
        return np.nanmin(data),np.nanmax(data)


def coord2str(ra,dec,output='decimal',frame='icrs'):
    """ 
    ra,dec coordinate to string
    
    eg. input ra=215.53478433,dec=-0.21847319 whose type is a number, default unit is deg.
    input '14h22m08.3482s' '-00d13m06.5035s' whose type is string.
    
    input string like:  '00h42m30s', '+41d12m00s' , '00h42.5m', '+41d12m' are all fine
    output: 'decimal','dms','hmsdms'
    """
    from astropy.coordinates import SkyCoord  
    if isinstance(ra,str):
        c=SkyCoord(ra,dec,frame=frame)
    else:
        from astropy import units as u
        c=SkyCoord(ra=ra*u.degree,dec=dec*u.degree,frame=frame)

    return c.to_string(output) 

def line_set(ax,xlabel,ylabel,xlim=None,ylim=None,legend=False,
             title=None,loc='best', size = None,frameon=True):
    """
    plot settings
    """
    if size is None:
        size = {'mz': 1,   # Set thickness of the tick marks
                'lz': 3,   # Set length of the tick marks
                'lbz': 14,  # Set label size
                'tkz': 12,  # Set tick size
                }
    mz = size['mz']; lz = size['lz']
    lbz = size['lbz']; tkz = size['tkz']
    # Make tick lines thicker
    for l in ax.get_xticklines():
        l.set_markersize(lz)
        l.set_markeredgewidth(mz)
    for l in ax.get_yticklines():
        l.set_markersize(lz)
        l.set_markeredgewidth(mz)

    # Make figure box thicker
    for s in ax.spines.values():
        s.set_linewidth(mz)
    ax.minorticks_on()
    ax.tick_params(labelsize=tkz)
    if xlim is not None: ax.set_xlim(xlim)
    if ylim is not None: ax.set_ylim(ylim)
    if xlabel is not None: ax.set_xlabel(xlabel,fontsize=lbz)
    if ylabel is not None: ax.set_ylabel(ylabel,fontsize=lbz)
    if legend: ax.legend(loc = loc,fontsize=tkz, frameon=frameon)
    if title is not None: ax.set_title(title,fontsize=lbz)

def reproject_fits(name,reproj_target,beam_target,result,method,slabm0):
    """
    see reproj_convolve_cube, reproj_convolve_2D
    """
    rj_target_hd = reproj_target.hd0
    new_header = rj_target_hd.copy() 
    
    beam = beam_target.beam
    new_header['BMAJ'] = beam.major.to_value(u.deg)                                                 
    new_header['BMIN'] = beam.minor.to_value(u.deg)                                                  
    new_header['BPA']  = beam.pa.to_value(u.deg)
    new_header['BEAM'] = f'Beam: BMAJ={beam.major.to_value(u.arcsec)} arcsec BMIN={beam.minor.to_value(u.arcsec)} arcsec BPA ={beam.pa.to_value(u.deg)} deg'
    

    if method == 'montage':
        """
        https://montage-wrapper.readthedocs.io/en/latest/
        """
        import montage_wrapper as mt #[sudo apt get montage] together
        new_header.totextfile('reproj_target_header.hdr', endcard=True, overwrite=True)
        if slabm0:    
            mt.reproject(name, result, header='reproj_target_header.hdr')
        else:
            mt.reproject_cube(name, result, header='reproj_target_header.hdr')
        os.remove('reproj_target_header.hdr')

    elif method == 'kapteyn':
        """
        https://www.astro.rug.nl/software/kapteyn/maputilstutorial.html#re-projections-and-image-overlays
        """
        if slabm0:
            from kapteyn import maputils
            hdl = fits.open(name)
            hdr = hdl[0].header
            if hdr['NAXIS'] != 2:
                if hdr['NAXIS'] == 3:
                    data = hdl[0].data[0,:,:]
                    del hdr['*3'];del new_header['*3']
                elif hdr['NAXIS'] == 4:
                    data = hdl[0].data[0,0,:,:]
                    del hdr['*3'];del hdr['*4'];del new_header['*3'];del new_header['*4']
                hdr['NAXIS'] = 2 ; new_header['NAXIS'] = 2

                map_self=maputils.FITSimage(externaldata=data,externalheader=hdr)
            else:
                map_self=maputils.FITSimage(name)
            hdl.close() 

            rj_target_data = reproj_target.data
            map_rj_target = maputils.FITSimage(externaldata=rj_target_data,externalheader=rj_target_hd)

            map_reprojected = map_self.reproject_to(rj_target_hd)
            fits.writeto(result,map_reprojected.dat,new_header,overwrite=True)

        else:
            raise ValueError("kapteyn method doesn't support cube reproject now, please use montage.")
    elif method == 'reproject':
        """
        https://reproject.readthedocs.io/en/stable/celestial.html#adaptive-resampling
        """
        if slabm0:
            from reproject import reproject_adaptive
            hdl = fits.open(name)
            hdr = hdl[0].header

            if hdr['NAXIS'] != 2:
                if hdr['NAXIS'] == 3:
                    data = hdl[0].data[0,:,:]
                    del new_header['*3']
                elif hdr['NAXIS'] == 4:
                    data = hdl[0].data[0,0,:,:]
                    del new_header['*3'];del new_header['*4']
                wcs_ = wcs.WCS(hdr).celestial
                reprojected_data,footprint = reproject_adaptive((data,wcs_),rj_target_hd)
                new_header['NAXIS'] = 2
            else:    
                reprojected_data,footprint = reproject_adaptive(hdl,rj_target_hd)

            fits.writeto(result,reprojected_data,new_header,overwrite=True)
            hdl.close()
        else:
            raise ValueError("reproject method doesn't support cube reproject now, please use montage.")
    else:
        raise ValueError("A method not support!")
        

################################ ellipse_fit function 1 ###################################
def ellipse_fit(x,y):
    '''
    Reference:
    [1] Halir R., Flusser J. 'Numerically Stable Direct Least Squares
    Fitting of Ellipses'
    [2] Weisstein, Eric W. "Ellipse." From MathWorld--A Wolfram Web Resource.
    http://mathworld.wolfram.com/Ellipse.html
    [3]Ben Hammel, & Nick Sullivan-Molina. (2020, March 21). bdhammel/
    least-squares-ellipse-fitting: v2.0.0 (Version v2.0.0). Zenodo. 
    http://doi.org/10.5281/zenodo.3723294 
    https://github.com/bdhammel/least-squares-ellipse-fitting
    '''
    #[1]Algorithm
    D1 = np.mat(np.vstack([x**2, x * y, y**2]).T)
    D2 = np.mat(np.vstack([x,y,np.ones(len(x))]).T)

    S1 = D1.T * D1
    S2 = D1.T * D2
    S3 = D2.T * D2

    C1=np.mat(np.zeros((3,3)))
    C1[0,2]=2;C1[1,1]=-1;C1[2,0]=2

    M = C1.I * (S1 - S2*S3.I*S2.T)
    eig_val,eig_arr = np.linalg.eig(M)

    DELTA=4*np.multiply(eig_arr[0, :], eig_arr[2, :])- np.power(eig_arr[1, :], 2)
    a1 = eig_arr[:, np.nonzero(DELTA > 0)[0]]
    a2 = -S3.I * S2.T * a1

    para=np.asarray((np.vstack((a1,a2)))).flatten() 
    
    center,major,minor,phi = get_param(para)

    
    return center,major,minor,phi,para

def get_param(para):
    #[2] ellipse equation
    a,b,c,d,f,g = para
    b = b/2 ;d = d/2; f=f/2

    x0 = (c*d - b*f)/(b**2 - a*c)
    y0 = (a*f - b*d)/(b**2 - a*c)
    center = [x0, y0]

    n = 2*(a*f**2 + c*d**2 + g*b**2 - 2*b*d*f - a*c*g)
    d1 = (b*b - a*c)*((c-a)*np.sqrt(1+ 4*b*b/((a-c)*(a-c)))-(c+a))
    d2 = (b*b - a*c)*((a-c)*np.sqrt(1 + 4*b*b/((a-c)*(a-c)))-(c+a))
    major = np.sqrt(n/d1)
    minor = np.sqrt(n/d2)
    phi = .5*np.arctan((2*b)/(a-c))
    
    return center,major,minor,phi

################################ ellipse_fit function 2 ###################################
def kapteyn_ellipse_fit(x,y):
    """
    https://www.astro.rug.nl/software/kapteyn/kmpfittutorial.html#special-topics
    """
    from kapteyn import kmpfit
    from math import sin, cos, radians, sqrt, atan

    def getestimates( x, y ):
        """
        Method described in http://en.wikipedia.org/wiki/Image_moments
        in section 'Raw moments' and 'central moments'.
        Note that we work with scalars and not with arrays. Therefore
        we use some functions from the math module because the are 
        faster for scalars
        """
        m00 = len(x)
        m10 = np.add.reduce(x)
        m01 = np.add.reduce(y) 
        m20 = np.add.reduce(x*x) 
        m02 = np.add.reduce(y*y) 
        m11 = np.add.reduce(x*y) 

        Xav = m10/m00
        Yav = m01/m00

        mu20 = m20/m00 - Xav*Xav
        mu02 = m02/m00 - Yav*Yav
        mu11 = m11/m00 - Xav*Yav

        theta = (180.0/np.pi) * (0.5 * atan(-2.0*mu11/(mu02-mu20)))
        if (mu20 < mu02):                   # mu20 must be maximum
            (mu20,mu02) = (mu02,mu20)        # Swap these values
            theta += 90.0

        d1 = 0.5 * (mu20+mu02)
        d2 = 0.5 * sqrt( 4.0*mu11*mu11 + (mu20-mu02)**2.0 )
        maj = sqrt(d1+d2)
        minm = sqrt(d1-d2)
        return (Xav, Yav, maj, minm, theta)


    def func(p, data):
        """
        Calculate z = (x/maj)**2 + (y/min)**2
        Note that z = 1 is an ellipse
        """
        x0, y0, major, minor, pa = p
        x, y = data
        pa   = radians(pa)   
        sinP = sin(-pa)
        cosP = cos(-pa)    
        xt = x - x0
        yt = y - y0    
        xr = xt * cosP - yt * sinP
        yr = xt * sinP + yt * cosP
        return (xr/major)**2 + (yr/minor)**2

    def residuals(p, data):
        """
        Note that the function calculates the height z of the 3d landscape
        of the function z = (x/maj)**2 + (y/min)**2.
        An ellipse is defined where z=1. The residuals we use
        for the least squares fit to minimize, is then 1-z
        """
        x, y = data
        return 1.0 - func(p,data)
    
    def get_parameter(xcf, ycf, majorf, minorf, angf):
        x0, y0, a, b, t =xcf, ycf, majorf, minorf, -np.deg2rad(angf-180) 
        
        A = a**2*sin(t)**2+b**2*cos(t)**2
        B = 2*(a**2-b**2)*sin(t)*cos(t)
        C = a**2*cos(t)**2+b**2*sin(t)**2
        D = A*(-2*x0)-y0*B
        F = C*(-2*y0) -x0*B
        G = -a**2*b**2 + A*x0**2 + C* y0**2 + B*x0*y0
        
        para = np.array([A,B,C,D,F,G])
        return para

    p0 = getestimates(x, y)

    fitter = kmpfit.Fitter(residuals=residuals, data=(x,y))
    fitter.fit(params0=p0)
    """
    print("\n========= Fit results ellipse model ==========")
    print("Initial params:", fitter.params0)
    print("Params:        ", fitter.params)
    print("Iterations:    ", fitter.niter)
    print("Function ev:   ", fitter.nfev)
    print("Uncertainties: ", fitter.xerror)
    print("dof:           ", fitter.dof)
    print("chi^2, rchi2:  ", fitter.chi2_min, fitter.rchi2_min)
    print("stderr:        ", fitter.stderr)
    print("Status:        ", fitter.status)
    """
    xcf, ycf, majorf, minorf, angf = fitter.params 
    
    para = get_parameter(xcf, ycf, majorf, minorf, angf)
    phi = np.deg2rad(angf-180)
    
    return [xcf, ycf], majorf, minorf,  phi , para

    
        



