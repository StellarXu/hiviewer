#!/usr/bin/env python
# coding: utf-8

"""
HIsource_fio.py
Created on someday November, 2020
@author: Xu Chen, Jilin University / NAOC / UCAS.

A simple of version of HIsource.py for hiviewer, is used to fit contour and plot pv diagram.

"""


import numpy as np
from matplotlib import pyplot as plt
import hiviewer as FIO
from astropy import log
import heapq

from .utils import line_set

class HIsource(object):
    def __init__(self,subcube,subm0):
        self.cube = subcube
        self.m0 = subm0
    
    def find_contour_fit(self,con,center_thld=4,function = 'default'):
        """
        find a contour to fit with an ellipse
        """
        from numpy.linalg import LinAlgError
        
        len_list=list(map(len,con.allsegs[0]))
        if len(len_list) < 4: # try 4 times at most
            contour_num_max = len(len_list)
        else:
            contour_num_max = 4
        
        len_list_sort=heapq.nlargest(contour_num_max,len_list)
        for contour_num in range(contour_num_max):
            print(f"fitting contour {contour_num}")
            loc = np.where(np.array(len_list)==len_list_sort[contour_num])[0][0]

            con_loc=con.allsegs[0][loc]#[np.argmax(len_list)]
            conx_,cony_=con_loc[:,0],con_loc[:,1]
            ################## use ellipse to fit the contour ####################
            try:
                if function == 'default':
                    from .utils import ellipse_fit
                    center_,major_,minor_,phi_,para_=ellipse_fit(conx_,cony_)
                elif function == 'kapteyn':
                    from .utils import kapteyn_ellipse_fit
                    center_,major_,minor_,phi_,para_=kapteyn_ellipse_fit(conx_,cony_)

                if np.sqrt(sum((center_-self.m0.crpix)**2)) < center_thld:
                    center,major,minor,phi,para,conx,cony = center_,major_,minor_,phi_,para_,conx_,cony_
                    print("Found a near-center ellipse :D  hhh")# guess the target should stay in central area
                    
                    return center,major,minor,phi,para,conx,cony
                
                else:
                    print("Far from center.Next try...")

            except LinAlgError as L:
                print(f"{repr(L)}. Can't fit this ellipse. Next try...")
            except ValueError as V:
                print(f"{repr(V)}. Can't fit this ellipse. Next try...")
                
        log.warn("failed  ellipse fitting !")
        return []
            

    def fit_contour(self ,factor=None,levels=None, line_length = 10,
                    vmin_max = None,clabel=False,figsize=(6, 5),
                    cbar_label='Jy/beam', picname ='./',save = False,cmap='gray',alpha = 0.9,
                    function = 'default',**kwargs):
        """ 
        plot the image, contours, ellipse, and get fit paraments and major line. 
        """
        
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection=self.m0.wcs_cel)
                                                                 #plt.rcParams['ytick.direction']='in'
        if vmin_max is None:
            im=ax.imshow(self.m0.get_slice(),alpha=alpha,cmap=cmap)
        else:
            im=ax.imshow(self.m0.get_slice(),vmin=vmin_max[0],vmax=vmin_max[1],alpha=alpha,cmap=cmap)        
        cbar = plt.colorbar(im,pad=.01)
        cbar.set_label(f'{cbar_label}', size=14)
        ax.set_xlabel("RA", fontsize=14)
        ax.set_ylabel("DEC", fontsize=14)
        [ax.coords[ci].display_minor_ticks(True) for ci in range(2)]
        if levels != None:
            contours=ax.contour(self.m0.get_slice(),levels = levels,colors='yellow')
        ################## definite a level and plot a contour ##############
        peak=np.nanmax(self.m0.data)
        if factor == None:
            factor = 0.5
        print(f"contour level to fit: {peak:.3f} * {factor}")
        con = ax.contour(self.m0.get_slice(),levels = [peak*factor] ,colors='orange')# contour to fit
        
        # plot other contours
        con_half=ax.contour(self.m0.get_slice(),levels = [peak*0.5] ,colors='r')
        con_quad=ax.contour(self.m0.get_slice(),levels = [peak*0.25] ,colors='b')
        if clabel:
            [ax.clabel(c,fmt='%.1f') for c in [con,con_half,con_quad]]
        ax.axis('equal')
           
        ################# find a contour to fit with an ellipse #################
        try:
            center,major,minor,phi,para,conx,cony = self.find_contour_fit(con,function=function)
            # coefficients 
            self.coef = [center,major,minor,phi]
            self.center = center
            self.phi = phi
            print("\n========= Fit results ellipse model ==========")
            print("Center pixel:",center)
            print("SemiMajor pixel:",major)
            print("SemiMinor pixel:",minor)
            print("Phi rad:",phi)
            print("=============== Fit Finished ===================")
        except ValueError as V:
            self.success_fit = False
            plt.show()
            log.warn(f"{repr(V)}")
            return self
        ################## get pixels position inside the ellipse ###############
#        inside_px = self.pix_inside_ellipse(para,**kwargs)
        # plot pixels inside the ellipse
#        if plot_inside_pix:
#            ax.plot(inside_px[:,0],inside_px[:,1],'m.',label='inside pixels')
        ######################## plot major line ###############################
        ps0, ps1 = self.maj_line(line_length)
        ax.plot([ps0[0], ps1[0]], [ps0[1], ps1[1]], c='r', linewidth=1)
        
        # plot FAST center
        ax.plot(center[0],center[1],'bx',label='FAST HI')
            
        ########################### plot the ellipse  ############################
        from matplotlib.patches import Ellipse
        ellipse = Ellipse(xy=center, width=2*major, height=2*minor, angle=np.rad2deg(phi),
            edgecolor='g', fc='None',label='ellipse fit',lw=2, zorder=2)
        ax.add_patch(ellipse)
        
        line_set(ax, xlabel = 'RA (J2000)',ylabel = 'Dec (J2000)', )
        ax.legend(loc=2,bbox_to_anchor=(1.07,1.07),borderaxespad=2.5)
        plt.tight_layout()
        if save:
            plt.savefig(picname+'_contour.png',dpi=300,bbox_inches='tight')
        plt.show()   
        
        self.success_fit = True
        return self
    
    def maj_line(self, line_length = 10):
        data = self.m0.data
        ny,nx = data.shape
        center,_,_,pa = self.coef
        ix0, iy0 = center
        
        ps0 = [ix0 - (- line_length * np.sin(pa)), iy0 - line_length * np.cos(pa)]
        ps1 = [ix0 + (- line_length * np.sin(pa)), iy0 + line_length * np.cos(pa)]

        return ps0, ps1
    
    
    def pv_cut(self, center, pa, line_length):
        """
        code modified from Ren Zhiyuan, 2022 multiband astronomy courses.

        center: pixel center
        pa: angle in rad
        line_length: pv line length in pixel
        """
        rsig = 3      # pixel numbers
        # rmax = 1         # searching range to find the int-pixel around the input center.

        cube_data = self.cube.data

        ix0, iy0 = center
        posi = np.arange(-line_length, line_length, 0.1)
        npos = len(posi)

        nv, ny, nx = np.shape(cube_data)
        pv_slice = np.zeros((nv, npos))   # array to restore pv plot
        # Tb_near = m0_data[int(iy0)-rmax:int(iy0)+rmax+1, int(ix0):int(ix0)+rmax+1]
        # [iy1, ix1] = np.where(Tb_near == Tb_near.max())   
        # find the real maximum position:
        # iy0, ix0 = int(iy0) + (iy1[0] - rmax), int(ix0) + (ix1[0] - rmax)
        ix0, iy0 = int(ix0), int(iy0)
        for i in np.arange(npos):
            weight_sum = 0
            # i-c along the slice:
            ixc = ix0 - posi[i] * np.sin(pa)
            iyc = iy0 + posi[i] * np.cos(pa)
            for j in np.arange(-2, 3):  # weighted average of the surrounding pixels
                for k in np.arange(-2, 3):
                    try:
                        ixs, iys = int(ixc) + j, int(iyc) + k
                        rs = np.sqrt((ixs - ixc)**2 + (iys - iyc)**2) # distance weight
                        weight = np.exp(-(rs)**2 / (2 * rsig **2))
                        pv_slice[:, i] += weight * cube_data[:, iys, ixs]
                        weight_sum += weight
                    except IndexError:
                        continue
            pv_slice[:, i] /= weight_sum
        return posi, pv_slice

    def pv_plot(self, center, pa, line_length, xylim = None, ovlim = None, save_name = None):
        """
        code modified from Ren Zhiyuan, 2022 multiband astronomy courses.

        center: pixel center
        pa: angle in rad
        line_length: pv line length in pixel
        xylim: xlim and ylim, default [0, nx, 0, ny]
        save_name: str
        """
        ix0, iy0 = center

        posi, pv_slice = self.pv_cut(center, pa, line_length)

        m0_data = self.m0.data
        ny, nx = np.shape(m0_data)

        ps0 = [ix0 - (-line_length * np.sin(pa)), iy0 - line_length * np.cos(pa)]
        ps1 = [ix0 + (-line_length * np.sin(pa)), iy0 + line_length * np.cos(pa)]

        vel = self.cube.velo.value
        p_mesh, v_mesh = np.meshgrid(posi, vel)

        lvs = m0_data.max()*np.arange(0.15, 1.1, 0.15)
        fig = plt.figure(figsize=(16, 7))  # figure window
        ax = fig.add_subplot(121)  # plot area with projection assigned
        # false-color image:
        vms = [m0_data.min(), m0_data.max()*0.9] 
        im = ax.imshow(m0_data, origin='lower', aspect='equal', cmap='terrain_r',
                         vmin= vms[0], vmax=vms[1], alpha=0.7)
        # line
        ax.plot([ps0[0], ps1[0]], [ps0[1], ps1[1]], c='black', linewidth=1)
        # contour levels:
        ax.contour(m0_data, levels=lvs, colors='white', origin='lower', linewidths=0.9)
        # center
        ax.plot(ix0, iy0,'rx')
        if xylim is not None:
            xlim = xylim[:2]; ylim = xylim[-2:]
        else:
            xlim = None; ylim = None
        line_set(ax, xlabel = 'RA [pixel]',ylabel = 'Dec [pixel]', xlim=xlim,ylim=ylim,)
        plt.colorbar(im, pad = 0.01)    

        ax = fig.add_subplot(122)
        # ax.imshow(pv_slice, origin='lower', aspect='auto')
        im = ax.pcolor(p_mesh, v_mesh, pv_slice, shading='auto')
        lvs = pv_slice.max()*np.arange(0.05, 1.1, 0.2)
        ax.contour(p_mesh, v_mesh, pv_slice, levels=lvs, colors='white', origin='lower', linewidths=0.9)
        
        if ovlim is not None:
            olim = ovlim[:2]; vlim = ovlim[-2:]
        else:
            olim = None; vlim = None
        line_set(ax, xlabel = 'Offset [pixel]',ylabel = 'Velocity [km/s]',xlim=olim,ylim=vlim,)
        plt.tight_layout()
        plt.colorbar(im, pad = 0.01)
        if save_name is not None:
            plt.savefig(save_name +'.png',bbox_inches='tight')
        plt.show()