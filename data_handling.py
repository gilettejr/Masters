#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 15:14:03 2019

@author: cameronrobertson
"""

import numpy as np
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord
from astropy import units as u

#class for reading in and plotting ascii file data

class asymphr:
    
    #class initialised, distance as optional input for plotting absolute magnitude CMDs. 
    
    def __init__(self,filename,distance):
        

        #distance and filename saved as class attributes
        
        self.filename = filename
        self.distance = distance
        
        #required file opened
        
        file = open(filename,'r')
        
        #data loaded and added as further attribute
    
        data = np.loadtxt(file)
    
        file.close()
        
        #data loaded for use with other methods
        
        self.data = data
        
    #method to assign column array titles and format ra and dec
        
    def loadascii(self):
        
        #function for converting hhmmss and ddmmss format into decimal degrees for ra and dec
        
        def coordtransform(rah,ram,ras,dd,dm,ds):
            
            #transformations carried out
            
            long = rah + ram/60 + ras/3600
            lat = dd + dm/60 + ds/3600
            
            #array containing transformed co-ordinates returned
    
            return(np.array([long,lat]))
        
        #class attribute called to assign columns
        
        data=self.data
        
        #columns assigned to variables for easer
        
        rah,ram,ras,dd,dm,ds = data[:,0],data[:,1],data[:,2],data[:,3],data[:,4],data[:,5]
        
        #function called to transform co-ordinates for plotting ease
        
        coords = coordtransform(rah,ram,ras,dd,dm,ds)
        
        #all remaining useful data columns assigned to variables
        
        xj,yj,jmag,jerr,jcis,xh,yh,hmag,herr,hcis,xk,yk,kmag,kerr,kcis = data[:,7],data[:,8],data[:,9],data[:,10],data[:,11],data[:,12],data[:,13],data[:,14],data[:,15],data[:,16],data[:,17],data[:,18],data[:,19],data[:,20],data[:,21]
        
        #all column variables assigned as class attributes 
        
        self.ra = coords[0]
        self.dec = coords[1]
        
        self.xj = xj
        self.yj = yj
        self.jmag = jmag
        self.jerr = jerr
        self.jcis = jcis
        self.xh = xh
        self.yh = yh
        self.hmag = hmag
        self.herr = herr
        self.hcis = hcis
        self.xk = xk
        self.yk = yk
        self.kmag = kmag
        self.kerr = kerr
        self.kcis = kcis
    
    #method for making cuts based on photometric quality
    
    def ciscuts(self):
        
        #loop to read cis number for each object
        
        for i in range(len(self.jcis)):
            
            #logic statement returns TRUE if cis number in any waveband is not -1 or -2
            
            if (self.jcis[i] != -1.0 and self.jcis[i]!=-2.0) or (self.hcis[i] != -1.0 and self.hcis[i]!=-2.0) or (self.kcis[i] != -1.0 and self.kcis[i]!=-2.0):
                
                #if TRUE returned, data not plotted
                
                self.ra[i] = np.nan
                self.dec[i] = np.nan
                self.xj[i] = np.nan
                self.yj[i] = np.nan
                self.jmag[i] = np.nan
                self.jerr[i] = np.nan
                self.jcis[i] = np.nan
                self.xh[i] = np.nan
                self.yh[i] = np.nan
                self.hmag[i] = np.nan
                self.herr[i] = np.nan
                self.hcis[i] = np.nan
                self.xk[i] = np.nan
                self.yk[i] = np.nan
                self.kmag[i] = np.nan
                self.kerr[i] = np.nan
                self.kcis[i] = np.nan
                
    #method for carrying out extinction corrections on data
                                        
    def extinction(self):
        
    #corrections made for each waveband, dependent on galaxy being opertated upon. Values taken from NED database    
        
        if self.filename == 'lot_n147.unique':
            
            self.jmag = self.jmag - 0.122
            self.hmag = self.hmag - 0.077
            self.kmag = self.kmag - 0.052
        
        elif self.filename =='lot_n185.unique':
            
            self.jmag = self.jmag - 0.130
            self.hmag = self.hmag - 0.083
            self.kmag = self.kmag - 0.056
            
            
        elif self.filename =='M32.asc':
            
            self.jmag = self.jmag - 0.044
            self.hmag = self.hmag - 0.028
            self.kmag = self.kmag - 0.019
            
            
        elif self.filename =='n205_two.asc':
            
            self.jmag = self.jmag - 0.044
            self.hmag = self.hmag - 0.028
            self.kmag = self.kmag - 0.019
            
            
        
        
           
        #print(np.shape((data[:,7:])))
        #print(np.shape(data))
        #xj = data[:,7]
        #print(xj)
        
#fileformat: na,na,na,na,na,na,na,x,y,jmag,jerr,jcis,x,y,hmag,herr,hcis,x,y,kmag,kerr,kcis

#class for reading gaia data tables in format described below

#gaiaformat:ra,raunc,dec,decunc,pmra,pmraunc,pmdec,pmdecunc
        
#class inherits from asymphr, initialises by loading datafile into numpy array
class gaiadata(asymphr):
    
    #method for assigning data into columns similarly to asymphr.loadascii
    
    def loadfile(self):
        data=self.data
        
        #columns assigned to variables
        
        self.ra,self.raunc,self.dec,self.decunc,self.pmra,self.pmraunc,self.pmdec,self.pmdecunc = data[:,1],data[:,2],data[:,3],data[:,4],data[:,5],data[:,6],data[:,7],data[:,8]
    
    #takes 1d array in the format[pmra,pmraunc,pmdec,pmdecunc]    
    
    def crossmatch(galaxypms)
        
#function plots k-j colour magnitude diagram from asymphr class execution

def kj_cmd(target):

#figure created
    
    plt.figure()
    
    #class methods carried out to load data and perform cuts and extinction corrections
    
    target.loadascii()
    target.ciscuts()
    target.extinction()
    
    #plot formatted and labelled
    
    plt.rc('axes',labelsize = 20)
    plt.scatter(target.jmag-target.kmag,target.kmag,s=3,marker='o',alpha = 0.5)
    plt.gca().invert_yaxis()
    plt.ylabel('K')
    plt.xlabel('J-K')
    plt.title(target.filename + ' CMD')
    
    plt.show()

    #function plots spatial distribution from asymphr class execution

def spatial_plot(target):
    
    plt.figure()
    
    #cuts performed, extinction correction not necessary
    
    target.loadascii()
    target.ciscuts()
    
    #plot formatted and labelled
    
    plt.rc('axes',labelsize = 20)
    plt.scatter(target.dec,target.ra,s=3,marker = 'o',alpha = 0.5)
    plt.ylabel('RA(J200)/degrees')
    plt.xlabel('Dec(J200)/degrees')
    plt.title(target.filename + ' Spatial Plot')
    plt.show()
    
def colour_colour(target):
    
    plt.figure()
    
    #colour cuts made, extinction corrected
    
    target.loadascii()
    target.ciscuts()
    target.extinction()
    
    #plot formatted and labelled
    
    plt.rc('axes',labelsize = 20)
    plt.scatter(target.jmag - target.hmag,target.hmag-target.kmag,s=1,marker=',',alpha = 0.5)
    plt.ylabel('H-K')
    plt.xlabel('J-K')
    plt.title(target.filename + 'Colour-colour diagram')
    
    plt.show()

#statements for running graphing functions. Galaxy chosen by initialising class before execution of functions

n147 = asymphr('lot_n147.unique',0)
#kj_cmd(n147)
#colour_colour(n147)
spatial_plot(n147)

#n185 = asymphr('lot_n185.unique',0)
#kj_cmd(n185)
#spatial_plot(n185)

#n205 = asymphr('N205_two.asc',0)
#kj_cmd(n205)
#spatial_plot(n205)

#m32 = asymphr('M32.asc',0)
#kj_cmd(m32)
#spatial_plot(m32)




    
            
        
        
