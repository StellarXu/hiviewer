#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
overlap.py
Created on someday Jan, 2020
@author: Xu Chen, Jilin University / NAOC / UCAS.

This code is used to: Overlap different images (or contours).

"""
#import matplotlib as mpl
#mpl.use('Agg')


from matplotlib import pyplot as plt
import numpy as np
from astropy import  log

from .utils import line_set

"""
%%% Overlap Functions
%%% overlap two contours or overlap the contour on the image
"""        
        
def overlap_two_contours(bk,fo,bklevels,folevels,c='orange',alphabk=0.6,alphafo=0.91,coord_format = None,
                         radec_range = None,figsize = (9, 6),save = False,picname ='./'):
    """
    overlap two contours
    
    Parameters
    ------------
    bk: background FitsPic example:FAST HI
    fo: foreground FitsPic
    levels and alphas belongs to the two contours
    c: the color of fore contour
    picname: save filepath or picname , like [picname+'_overlap_img_contour.png'].
             Also titled with it.
    radec_range: unit deg. eg.When 'radec_range = [22.7,24.1,30,31.5]', it only shows 
                 area ra ranges from 22.7degto 24.1deg, dec ranges from 30deg to 31.5deg. 
    save: if False, it won't save the pic. Default to True.
    """
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection=bk.wcs_cel)
    #a background pic only to limit the fig size 
    im_bk=ax.imshow(bk.get_slice() ,alpha=0)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    fo_transform = ax.get_transform(fo.wcs_cel)  
    # extract axes Transform information for the HI data
    con_bk=ax.contour(bk.get_slice(), levels=bklevels,alpha=alphabk,linestyles='solid',
                      colors='k',linewidths=1.5)
    con_fo=ax.contour(fo.get_slice(), levels=folevels,alpha=alphafo,linestyles='solid',
                      transform = fo_transform,colors=c,linewidths=1.5) 

        
    if coord_format is not None:
        ra=ax.coords['ra'];dec=ax.coords['dec'];ra.set_major_formatter(coord_format);dec.set_major_formatter(coord_format)
 
    ax.grid()
    if radec_range is not None:
        xlim1,ylim1 = bk.deg2pix(radec_range[1],radec_range[2])
        xlim2,ylim2 = bk.deg2pix(radec_range[0],radec_range[3])
        xlim = (xlim1,xlim2)
        ylim = (ylim1,ylim2)
    line_set(ax, xlabel = 'RA (J2000)',ylabel = 'Dec (J2000)', xlim=xlim,ylim=ylim, title = picname,)
    plt.tight_layout()
    if save:
        plt.savefig(picname+'_overlap_contours.png',dpi=300,bbox_inches='tight')
    #plt.show()

    
def overlap_img_contour(bk,fo,levels,c='orange',vmin_max=None,alpha=0.8,bar_label='label',coord_format = None,
                        clabel=False,radec_range = None,figsize = (9, 6),save = False,picname ='./'):
    """
    put fo contour on bk image.
    
    Parameters
    ------------
    bk: background.FitsPic, example: image of FAST HI
    fo: foreground.FitsPic
    levels: belongs to the  contour
    alpha: of the image
    c: the color of contour
    vmin_vmax: eg.(a,b) means that img color ranges a from b.
    bar_label: label of the colorbar
    clabel: if True, it shows contour with number.
    picname: save filepath or picname , like [picname+'_overlap_img_contour.png']
             Also titled with it.
    radec_range: unit deg. eg.When 'radec_range = [22.7,24.1,30,31.5]', it only shows 
                 area ra ranges from 22.7degto 24.1deg, dec ranges from 30deg to 31.5deg. 
    save: if False, it won't save the pic. Default to True.
    """
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection=bk.wcs_cel)

    if vmin_max is None:
        im=ax.imshow(bk.get_slice())
    else:
        im=ax.imshow(bk.get_slice(),vmin=vmin_max[0],vmax=vmin_max[1])
    im_bar=plt.colorbar(im,pad=.01)
    im_bar.set_label(bar_label,size=15)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    hi_transform = ax.get_transform(fo.wcs_cel)  
    con_fo=ax.contour(fo.get_slice(), alpha=alpha, levels=levels,linestyles='solid',
                   transform = hi_transform,colors=c,linewidths=1.3) 
    if clabel:
        ax.clabel(con_fo,fmt='%.f')

    ax.grid(color = 'white', ls = 'dotted', lw = 2,alpha=0.5)   
    
    if coord_format is not None:
        	ra=ax.coords['ra'];dec=ax.coords['dec'];ra.set_major_formatter(coord_format);dec.set_major_formatter(coord_format)
    
    if radec_range is not None:
        xlim1,ylim1 = bk.deg2pix(radec_range[1],radec_range[2])
        xlim2,ylim2 = bk.deg2pix(radec_range[0],radec_range[3])
        xlim = (xlim1,xlim2)
        ylim = (ylim1,ylim2)
    line_set(ax, xlabel = 'RA (J2000)',ylabel = 'Dec (J2000)', xlim=xlim,ylim=ylim, title = picname,)
    plt.tight_layout()
    if save:
        plt.savefig(picname+'_overlap_img_contour.png',dpi=300,bbox_inches='tight')
    #plt.show()
    
    
def two_contours_with_peaks(bk,fo,bklevels,folevels,min_distance = 2,threshold_rel=0.2, c='m',alphabk=0.3,alphafo=0.5,
                         radec_range = None,figsize = (9, 6),save = False,picname ='./',coord_format = None,):
    """
    overlap two contours and mark the local peaks and valleys.
    
    Parameters
    ------------
    bk: background 2D FitsPic example:FAST HI moment0
    fo: foreground 2D FitsPic
    levels and alphas belongs to the two contours
    min_distance: int, Minimum number of pixels separating peaks in a region of 2 * min_distance + 1 
                 (i.e. peaks are separated by at least min_distance). To find the maximum number of 
                 peaks, use min_distance=1.
    threshold_rel: float, Minimum intensity of peaks, calculated as max(image) * threshold_rel.
                   Set to None if don't need.
    c: the color of fore contour
    picname: save filepath or picname , like [picname+'_overlap_img_contour.png'].
             Also titled with it.
    radec_range: unit deg. eg.When 'radec_range = [22.7,24.1,30,31.5]', it only shows 
                 area ra ranges from 22.7degto 24.1deg, dec ranges from 30deg to 31.5deg. 
    save: if False, it won't save the pic. Default to True.
    """
    from skimage.feature import peak_local_max
    
    if (bk.ndim,fo.ndim) != (2,2) :
        raise(ValueError("Should input 2D FitsPic."))
    
    def pix2pix(bk,fo,coord_fo):
        """Transforms foreground pixel coordinates to background pixel coordinates."""
        ra,dec = fo.wcs_obj.all_pix2world(coord_fo[:, 1],coord_fo[:, 0],0)
        #.all_pix2world(y,x,0),0 means Numpy and C standard index begins with 0
        px1,py1 = bk.wcs_obj.all_world2pix(ra,dec,0)
        px = np.rint(px1).astype(int)
        py = np.rint(py1).astype(int)
        return px,py


    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection=bk.wcs_cel)
    #a background pic only to limit the fig size 
    img_bk=ax.imshow(bk.get_slice() ,alpha=0)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    [ax.coords[ci].display_minor_ticks(True) for ci in range(2)]
    #extract axes Transform information for the HI data
    fo_transform = ax.get_transform(fo.wcs_cel)  
    con_bk=ax.contour(bk.get_slice(), levels=bklevels,alpha=alphabk,linestyles='solid',
                      colors='k',linewidths=1.5)
    con_fo=ax.contour(fo.get_slice(), levels=folevels,alpha=alphafo,linestyles='solid',
                      transform = fo_transform,colors=c,linewidths=1.5) 
    #fill nans
    im_bk = bk.fill_nans_data()
    im_fo = fo.fill_nans_data()
    im_bk2 = -im_bk
    im_fo2 = -im_fo
    #find peaks and valleys
    coord_bk = peak_local_max(im_bk, min_distance=min_distance,threshold_rel=threshold_rel)
    coord_fo = peak_local_max(im_fo, min_distance=min_distance,threshold_rel=threshold_rel)
    coord_bk2 = peak_local_max(im_bk2, min_distance=min_distance)
    coord_fo2 = peak_local_max(im_fo2, min_distance=min_distance)
    #transform pixel coordinates
    px_fo,py_fo = pix2pix(bk,fo,coord_fo)
    px_fo2,py_fo2 = pix2pix(bk,fo,coord_fo2)
    
    px_bk,py_bk = coord_bk[:, 1],coord_bk[:, 0]
    px_bk2,py_bk2 = coord_bk2[:, 1],coord_bk2[:, 0]
    #mark points
    ax.scatter(px_bk,py_bk,marker='o',c='r',s=150,label='background peak',zorder=2)
    ax.scatter(px_fo,py_fo,marker='o',c='g',label='foreground peak',zorder=2)
    ax.scatter(px_bk2,py_bk2,marker='o',c='b',s=150,label='background valley',zorder=2)
    ax.scatter(px_fo2,py_fo2,marker='o',c='orange',label='foreground valley',zorder=2)
        
    if coord_format is not None:
        	ra=ax.coords['ra'];dec=ax.coords['dec'];ra.set_major_formatter(coord_format);dec.set_major_formatter(coord_format)
    
    ax.grid()
    if radec_range is not None:
        xlim1,ylim1 = bk.deg2pix(radec_range[1],radec_range[2])
        xlim2,ylim2 = bk.deg2pix(radec_range[0],radec_range[3])
        xlim = (xlim1,xlim2)
        ylim = (ylim1,ylim2)
    line_set(ax, xlabel = 'RA (J2000)',ylabel = 'Dec (J2000)', xlim=xlim,ylim=ylim,
             title = picname, legend = True)
    plt.tight_layout()
    if save:
        plt.savefig(picname+'_peaks.png',dpi=300,bbox_inches='tight')
    log.info(" The first you input(bk =) is background.The second you input(fo =) is foreround. "
         "When comparing these two, the valley floor positions are less accurate than the peaks, \
         because their gradient is smaller.")
    #plt.show()
    
