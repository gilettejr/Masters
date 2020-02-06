#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 21:11:28 2019

@author: cameronrobertson
"""
import numpy as np
from astropy.coordinates import SkyCoord
from dustmaps.sfd import SFDQuery


class asymphr:
    
    #class initialised, distance as optional input for plotting absolute magnitude CMDs. 
    
    def __init__(self,filename,distance):
        

        #distance and filename saved as class attributes
        
        self.filename = filename
        self.distance = distance
        
        #ngc147data
        
        self.pmra = 0
        self.pmraunc = 0
        self.pmdec = 0
        self.pmdecunc = 0
        self.raunc = 0.0005
        self.decunc = 0.0005
        
        self.xi=0
        self.eta=0
        
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
            
            long = (rah + np.divide(ram,60) + np.divide(ras,3600))*15
            lat = dd + np.divide(dm,60) + np.divide(ds,3600)
            
            #array containing transformed co-ordinates returned
    
            return(np.array([long,lat]))
        
        #class attribute called to assign columns
        
        data=self.data
        
        #columns assigned to variables for ease
        
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
        
        #index array created for crossmmatching functionality
        
        index = []
        
        
        for i in range(len(self.ra)):
            index.append(i)
        self.index = np.array(index)
            
    
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
    
    #redundant function with less rigourous cuts than ciscuts
    
    def lowciscuts(self):
        
        #loop to read cis number for each object
        
        for i in range(len(self.jcis)):
            
            #logic statement returns TRUE if cis number in all wavebands is not -1 or -2
            
            if (self.jcis[i] != -1.0 and self.jcis[i]!=-2.0) and (self.hcis[i] != -1.0 and self.hcis[i]!=-2.0) and (self.kcis[i] != -1.0 and self.kcis[i]!=-2.0):
                
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
            
    #more thourough extinction correction method, using sfd dustmap for each star individually        
            
    def sbsextinction(self):
        
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
        
        jext = sfdred * 0.709
        hext = sfdred * 0.449
        kext = sfdred * 0.302
        
        #extinction corrections carried out on class attributes
        
        self.jmag=self.jmag - jext
        self.hmag=self.hmag - hext
        self.kmag=self.kmag - kext
    
    #method to carry out star by star extinction on object using SFD dustmap. Takes object target as input
    
        print('UKIRT Extinction corrections done')
    
    def sbsextinction_on_frame(self,target):
        
        
        #attributes set as variables for ease
        
        ra = target.ra
        dec = target.dec
        
        #astropy skycoord called, units and frame set
        
        coords=SkyCoord(ra,dec,unit='deg',frame='icrs')
        
        #sfd map chosen
        
        sfd=SFDQuery()
        
        #E(B-V) for each star loaded into array from map
        
        sfdred=sfd(coords)
           
        #corrections from table 6 of Schlafly and Finkbeiner(2011) made
        
        jext = sfdred * 0.709
        hext = sfdred * 0.449
        kext = sfdred * 0.302
        
        #extinction corrections carried out on class attributes
        
        target.jmag=target.jmag - jext
        target.hmag=target.hmag - hext
        target.kmag=target.kmag - kext
        
    #method for adding xi and eta tanget co-ordinates as attributes
    
    def create_tangent_coords(self,tangentra,tangentdec):
        
        #ra and dec attributes converted to variables in radians
        
        ra = np.radians(self.ra)
        dec = np.radians(self.dec)
        
        #tangent co-ordinates also converted to radians
        
        tanra = np.radians(tangentra)
        tandec = np.radians(tangentdec)
        
        #conversion for xi carried out
        
        xi = (np.cos(dec)*np.sin(ra-tanra))/(np.sin(dec)*np.sin(tandec) + np.cos(dec)*np.cos(tandec)*np.cos(ra-tanra))
        
        #conversion for eta carried out
        
        eta = (np.sin(dec)*np.cos(tandec)-np.cos(dec)*np.sin(tandec)*np.cos(ra-tanra))/(np.sin(dec)*np.cos(tandec)+np.cos(dec)*np.sin(tandec)*np.cos(ra-tanra))
        
        #co-ordinates converted to degrees and set as attributes
        
        self.xi = xi * (180/np.pi)
        self.eta = eta * (180/np.pi)