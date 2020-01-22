#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 14:43:33 2020

@author: cameronrobertson
"""

from astropy.io import ascii
from astropy.visualization import quantity_support
quantity_support()
from astropy.coordinates import SkyCoord
from astropy import units as u
from dustmaps.sfd import SFDQuery
from asymphr import asymphr
from graphing_class import graphs
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class crossall:
    
    def __init__(self,filename):
        self.filename=filename
    
    def panda_cls(self):
        
        data = pd.read_csv(self.filename)
        

        
        for i in range(len(data.g)):
            if (data.gcls[i] !=-1.0 and data.gcls[i] !=-2.0) or (data.icls[i] !=-1.0 and data.icls[i] !=-2.0):
                
                data.loc[i]=np.nan
                

                
        print(data)
                
                
        data.to_csv('NGC185_pandcut.csv')
    
    def make_cmd(self,infile):
        
        data=pd.read_csv(infile)
        
        self.ra=data.ra_1
        self.dec=data.dec_1
        self.index=data.id
        self.gmag=data.g
        self.imag=data.i
        self.irbluemag=data.threesix
        self.irredmag=data.fourfive
        
        hmag=[]
        jmag=[]
        kmag=[]
        
            
        
        
        a=asymphr('lot_n147.unique',0)
        a.loadascii()
        a.sbsextinction()
        
        hmag=[]
        jmag=[]
        kmag=[]
        
        for i in self.index:
            
            hmag.append(a.hmag[i])
            jmag.append(a.jmag[i])
            kmag.append(a.kmag[i])
        
        self.hmag=np.array(hmag)
        self.jmag=np.array(jmag)
        self.kmag=np.array(kmag)
        
        #plotter=graphs()
        #plotter.kj_cmd(self)
        
    def pandaviz_extinction(self):
        
        #ra and dec variables set from class attributes for ease
        
        ra = self.ra
        dec = self.dec
            
        #astropy skycoord called, units and frame set
            
        coords=SkyCoord(ra,dec,unit='deg',frame='icrs')
            
        #sfd map chosen
            
        sfd=SFDQuery()
            
        #E(B-V) for each star loaded into array from map
            
        sfdred=sfd(coords)
               
        #corrections from table 6 of Schlafly and Finkbeiner(2011) made
            
            
        gext = sfdred * 3.172
        iext = sfdred * 1.682
            
            
        #extinction corrections carried out on class attributes
            
        
            
        self.gmag=self.gmag - gext
        self.imag=self.imag - iext
        
        #method to carry out star by star extinction on object using SFD dustmap. Takes object target as input
        
        print('Pandas Extinction corrections done')
        
    def choose_stars(self,xmin,xmax,ymin,ymax):
        
        #NGC147 C cluster
        
        #1.8-2
        #16.6-16.9
        
        #NGC147 M cluster
        
        #1.110-1.135
        #17.50-17.56
        kmag=self.kmag
        jmag=self.jmag
        j_k=jmag-kmag
        
        indices=[]
        
        for i in range(len(self.index)):
            if xmin<j_k[i]<xmax and ymin<kmag[i]<ymax:
                indices.append(self.index[i])
                
        locind=[]
        
        for i in indices:
            
            locind.append((np.where(self.index==i)))
        
        self.sedindices=np.array(locind)
        print(self.sedindices)
        
    def plot_sed(self):
        
        gwave=4876.7
        iwave=7520.8
        jwave=12482.9
        hwave=16588.4
        kwave=21897.7
        irbluewave=36000
        irredwave=45000
        
        
        i=106
        
        lam=gwave * u.AA
        
        ydatmag=np.array([self.gmag[i],self.imag[i],self.jmag[i],self.hmag[i],self.kmag[i],self.irbluemag[i],self.irredmag[i]]) * u.mag
        gflux=self.gmag[i]*u.mag.to('Jy', u.spectral_density(lam))
        xdat=np.array([gwave,iwave,jwave,hwave,kwave,irbluewave,irredwave]) * u.AA
        

        
        #plt.gca().invert_yaxis()
        #plt.plot(xdat,ydatflux,label='Star 1')
        
        
        #plt.legend()
        #plt.show()
    
    
c = crossall('NGC147_pand.csv')
c.make_cmd('ngc147_sed.csv')
c.pandaviz_extinction()
c.choose_stars(1.8,2.0,16.6,16.9)
c.plot_sed()