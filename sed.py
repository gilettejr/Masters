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
        
        plotter=graphs()
        plotter.kj_cmd(self)
        
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
            
            locind.append(np.where(self.index==i)[0][0])
        self.sedindices=locind
        print(self.sedindices)
        
        plotter=graphs()
        plotter.kj_cmd_select(self,xmin,xmax,ymin,ymax)
        
        
    def plot_sed(self):
        
        def mag_to_flux(m,zpflux):
            
            f=zpflux*10**(-m/2.5)
            
            return f
        
        def angwave_to_HZ(lam):
            
            c=3*10**8
            lam=lam*10**-10
            
            freq=c/lam
            
            return freq
        
        gwave=4876.7
        iwave=7520.8
        jwave=12482.9
        hwave=16588.4
        kwave=21897.7
        irbluewave=36000
        irredwave=45000
        
        gzpoint=3893.0
        izpoint=2577.0
        jzpoint=1556.8
        hzpoint=1038.3
        kzpoint=644.1
        irbluezpoint=277.2
        irredzpoint=179.0
        
        zpoints=np.array([gzpoint,izpoint,jzpoint,hzpoint,kzpoint,irbluezpoint,irredzpoint])
        angwaves=np.array([gwave,iwave,jwave,hwave,kwave,irbluewave,irredwave])
        freqs= angwave_to_HZ(angwaves)
        micwaves=angwaves * 10 ** -4
        
        xticks=np.concatenate((micwaves,np.array([1])))
        
        plt.figure()
        
        j=self.sedindices
        
        for i in range(len(j)):
            
            print(i)
            
            mags=np.array([self.gmag[j[i]],self.imag[j[i]],self.jmag[j[i]],self.hmag[j[i]],self.kmag[j[i]],self.irbluemag[j[i]],self.irredmag[j[i]]])  
        
            fluxes=mag_to_flux(mags,zpoints)
        
            ydata=fluxes * freqs
        

            plt.plot(micwaves,ydata,marker='o',label = 'SED for NGC147 Star, SED Index ' + str(j[i]))
            
            if (i%11==0 and i!=0) or i==len(self.sedindices)-1:
                
                print('yeet')
            
                plt.yscale('log')
                plt.xscale('log')
                plt.ylabel('v$F_v$(Hz Jy)')
                plt.xlabel('wavelength(um)')
        
        #for i in range(len(micwaves)):
        
        #plt.xticks(xticks,['g','i','j','h','k','3.6','4.5',str(10**0)])
            
            #print(micwaves[i],ydata[i])
            
                plt.legend()
        
        

        
                plt.show()
                
                if i%11==0:
                    plt.figure()
                    

        

def main():
    
    print('Press i for interactive, m for m stars or c for c stars')
    
    choice=input()
    
    c = crossall('NGC147_pand.csv')
    c.make_cmd('ngc147_sed.csv')
    c.pandaviz_extinction()
    
    if choice=='i':
        
        print('Here is a CMD of all points crossmatched with Pandas, Spitzer and UKIRT')
    
            
        c = crossall('NGC147_pand.csv')
        c.make_cmd('ngc147_sed.csv')
        c.pandaviz_extinction()
        
        print('Please input coordinates for a box containing the stars you want to see an SED for\n')
        
        print('xmin: ')
        xmin=input()
        
        print('xmax: ')
        xmax=input()
        
        print('ymin: ')
        ymin=input()
        
        print('ymax: ')
        ymax=input()
    
        c.choose_stars(xmin,xmax,ymin,ymax)
        
    elif choice=='c':
        
        c.choose_stars(1.33,20,0,18)
        
    elif choice=='m':
        
        c.choose_stars(1,1.33,0,18)

    c.plot_sed()
    
main()