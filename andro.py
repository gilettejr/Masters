#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 15:31:10 2020

@author: cameronrobertson
"""

from crossmatching_utils import topcatcross,crossed_data
from HESSCMD import plotHess
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from iso_utils import import_isos

#takes in m31 file as input, applies cls and extinction cuts, saves to csv file
def clsextcut():
    
    f=topcatcross('lot_m31.unique',0)
        
    f.loadascii()
    f.ciscuts()
    f.sbsextinction()
    f.create_tangent_coords(10.6847083,41.26875)
    f.make_dataframe()
    
    data=f.wfdat
    
    data.to_parquet('m31_int_files/int1')

def topcrossed():
    rmatch=crossed_data()
    
    n147topcross = topcatcross('lot_m31.unique',0)
    rmatch.topmatch(n147topcross,'crossedm31.csv')
    n147topcross.delete_crossed_points()
    
    n147topcross.create_tangent_coords(10.6847083,41.26875)
    n147topcross.make_dataframe()
    
    data=n147topcross.wfdat
    
    print(data)
    
    data.to_parquet('m31_int_files/int2')
    

    
    
    


#index descriptors
#1 - data just with cls cut
#2 - data with cls and topcat cuts
#3 - bulge
#4 - bulge of top right stream, 1.46,2.1, 0.43, 0.73
#5 - sample of bottom left disc -1.60, -1.30,-0.80,-0.60
#6 - g1 clump, -1.92,-1.44,-1.7,-1.4


    

class m31data:
    
    #class for andromeda methods
    def __init__(self,intno):
        
        data=pd.read_parquet('m31_int_files/int' + str(intno))
        self.data=data
        

        
     
    
    def plot_cmd(self):
        
        data=self.data
        plt.figure()
        plt.rc('axes',labelsize=15)
        plt.scatter(data.jmag-data.kmag,data.kmag,color='black',marker='.')
        plt.gca().invert_yaxis()
        plt.show()
        plt.ylabel('K$_0$')
        plt.xlabel('J$_0$-K$_0$')
    
    def plot_hess(self):
        data=self.data
        plotHess(data.kmag.dropna(),data.jmag.dropna()-data.kmag.dropna())
        
    def plot_tan(self):
        
        data=self.data
        #plt.scatter(data.xi,data.eta,marker='.')
        #plotHess(data.eta.dropna(),data.xi.dropna())
        plt.figure()
        plt.rc('axes',labelsize=15)
        plt.plot(data.eta,data.xi,color='black',linestyle='none',marker=',')
        plt.gca().set_ylabel(r'$\eta$/degrees')
        plt.gca().set_xlabel(r'$\xi$/degrees')
        #plt.xlim(-2.1,2.1)
        #plt.ylim(-2.1,2.1)
        
        #plt.xlim(-3,3)
        #plt.ylim(-5,3)
        plt.gca().invert_xaxis()
        
    def plot_spatial(self):
        
        data=self.data
        plt.figure()
        plt.plot(data.ra,data.dec,color='blue',linestyle='none',marker=',')
        plt.gca().set_xlabel('Ra')
        plt.gca().set_ylabel('Dec')
        
    def plot_cc(self):
        
        data=self.data
        
        jherr=np.sqrt(data.jerr**2 + data.herr**2)
        hkerr=np.sqrt(data.herr**2 + data.kerr**2)
        
        plt.figure()
        plt.rc('axes',labelsize=15)
        plt.scatter(data.hmag-data.kmag,data.jmag-data.hmag,color='black',s=3,marker='o')
        
        plt.ylabel('$J_0$-$H_0$')
        plt.xlabel('$H_0$-$K_0$')
        
    def plot_ccerr(self):
        
        data=self.data
        
        jherr=np.sqrt(data.jerr**2 + data.herr**2)
        hkerr=np.sqrt(data.herr**2 + data.kerr**2)
        
        plt.figure()
        
        plt.errorbar(data.hmag-data.kmag,data.jmag-data.hmag,xerr=hkerr,yerr=jherr,ecolor='red',color='black',linestyle='none',marker='.')
        
        plt.ylabel('$J_0$-$H_0$')
        plt.xlabel('$H_0$-$K_0$')
        
    def make_anuli(self,slice_size):
            
        data=self.data
        
        fig,ax=plt.subplots(1)
        
        #ellipse data
        
        xheight=0.93
        yheight=1.00
        
        xwidth=0.095
        ywidth=0.665
        
        ell_center=(0,0)
        ell_width=np.sqrt(xwidth**2+ywidth**2)
        ell_height=np.sqrt(xheight**2 + yheight**2)
        angle=np.arctan(1/0.93)
        angle=360-np.degrees(angle)
        
        ellipse=patches.Ellipse(ell_center,ell_width,ell_height,angle)
        #ax.add_patch(ellipse)
        
        cos_angle = np.cos(np.radians(180.-angle))
        sin_angle = np.sin(np.radians(180.-angle))
        
        x=data.xi.dropna()
        y=data.eta.dropna()
        
        xc = x - ell_center[0]
        yc = y - ell_center[1]
        
        xct = xc * cos_angle - yc * sin_angle
        yct = xc * sin_angle + yc * cos_angle 
        
        rad_cc = (xct**2/(ell_width/2.)**2) + (yct**2/(ell_height/2.)**2)
        
        colors_array = []
        
        for r in rad_cc:
            if r <= 1.:
                # point in ellipse
                colors_array.append('green')
            else:
                # point not in ellipse
                colors_array.append('blue')
        
        ax.scatter(x,y,c=colors_array,linewidths=0.3,alpha=1)
        
        plt.show()
    #bulge region, identified with index 3 on initiaisation
    def bulge_cut(self,radiusarcsec):
        
        data=self.data
        
        for i in range(len(data.ra)):
            
            if np.isnan(data.ra[i])==True:
                
                continue
            
            elif np.sqrt((data.ra[i]-10.6847083)**2 + (data.dec[i]-41.26875)**2) > radiusarcsec/3600:
                
                data.loc[i]=np.nan
                
        data.to_parquet('m31_int_files/int3')
    
    #box, take from tangent distribution
    def box_cut(self,x1,x2,y1,y2):
        
        data=self.data
        
        x=data.xi
        y=data.eta
        
        for i in range(len(data.ra)):
            
            if np.isnan(data.ra[i])==True:
                
                continue
            
            elif x[i]<x1 or x[i]>x2 or y[i]<y1 or y[i]>y2:
                
                data.loc[i]=np.nan
                
        data.to_parquet('m31_int_files/int6')
        
class andro_utils(m31data):
    
    def make_c_m(self,fore,jklim,hklim,jhlim,trgb):
        
        data=self.data
        mdata=data.copy()
        cdata=data.copy()
        
        jk=data.jmag-data.kmag
        hk=data.hmag-data.kmag
        jh=data.jmag-data.hmag
        
        for i in range(len(mdata.kmag)):
            
            if np.isnan(mdata.kmag[i])==True:
                
                continue
            
            elif jk[i]>jklim or jk[i]<fore or (hk[i]>hklim and jh[i]>jhlim) or mdata.kmag[i]>trgb:
                
                mdata.loc[i]=np.nan
                
        for i in range(len(cdata.kmag)):
            
            if np.isnan(cdata.kmag[i])==True:
                
                continue
            
            elif hk[i]<hklim or jk[i]<jklim or jh[i]<jhlim or cdata.kmag[i]>trgb:
                cdata.loc[i]=np.nan
                
        self.cdata=cdata
        self.mdata=mdata
        
    def find_FEH(self):
            
        def CM_to_FEH(CM):
        
            FEH=-1.39 -0.47*np.log10(CM)
        
            return(FEH)
        
        mdata=self.mdata
        cdata=self.cdata
        
        m_no_top=[]
        for i in mdata.kmag:
            if np.isnan(i)==False:
                m_no_top.append(0)
        c_no_top=[]   
        for j in cdata.kmag:
            if np.isnan(j)==False:
                c_no_top.append(0)
        
        C_top=len(c_no_top)
        M_top=len(m_no_top)
        
        ratio_top = C_top/M_top
        FEH_top=CM_to_FEH(ratio_top)
        
        
        print('Gaia crossmatch: C/M = ' + str(ratio_top) + ' and [Fe/H] = ' + str(FEH_top))
        
            
        
def main():
    #m31.box_cut(-1.92,-1.44,-1.7,-1.4)
    
    #m31.find_FEH()
    
    #m31=andro_utils(2)
    #m31.plot_cmd()
    
    bulge=andro_utils(3)
    #bulge.plot_tan()
    #bulge.find_FEH()
    #bulge.plot_cmd()
    bulge.plot_cmd()
    
    
    #tr=andro_utils(4)
    #tr.make_c_m(1.03,1.37,0.50,0.82,17.8)
    #tr.find_FEH()
    #tr.plot_cc()
    
    bl=andro_utils(5)
    bl.plot_cmd()
    #bl.plot_cc()
    #bl.make_c_m(1.02,1.37,0.47,0.78,17.8)
    #bl.find_FEH()
    #bl.plot_cc()
    #m31.box_cut(1.46,2.1,0.43,0.73)
    
    g1=andro_utils(6)
    g1.plot_cmd()
    
    #i=import_isos('andromeda','Isochrones/disc_app.dat')
    #i.overlay_agb_tip()
    
    #g1.plot_cc()
    #g1.make_c_m(0.97,1.43,0.5,0.77,17.8)
    ##g1.find_FEH()
    
    
    
main()

        