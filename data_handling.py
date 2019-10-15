#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 15:14:03 2019

@author: cameronrobertson
"""

import numpy as np
import matplotlib.pyplot as plt


#class for reading in and plotting ascii file data

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
            
            long = (rah + ram/60 + ras/3600)*15
            lat = dd + dm/60 + ds/3600
            
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

#gaiaformat:sourceid,ra,raunc,dec,decunc,pmra,pmraunc,pmdec,pmdecunc
        
#class initialises by loading datafile into numpy array
class gaiadata:
    
    def __init__(self,filename):
        self.filename = filename
        
        file  = open(filename,'r')
        
        #delimiter specified, since csv file being read in
        
        data = np.genfromtxt(file,delimiter = ',',missing_values='',filling_values=np.nan)
        
        file.close()
        #array assigned as data object
        self.data = data
    
    #method for assigning data into columns similarly to asymphr.loadascii, and carrying out initial data cuts
    
    def loadfile(self):
        data=self.data
        
        k=[]
        #columns assigned to variables
        #rows with high excess noise and badly measured parallaxes cut
        for i in range(len(data[:,1])):
            if(data[i,11])>2 or (data[i,6]<1):
                k.append('yeet')
                for j in range(len(data[i,:])):
                    data[i,j]=np.nan
            #rows with proper motions within uncertainty of 0 cut, distant sources
            elif(np.abs(data[i,7]/data[i,8])<1) or (np.abs(data[i,9]/data[i,10])<1):
                k.append('yeet')
                for j in range(len(data[i,:])):
                    data[i,j]=np.nan
                
                    

        
        
        self.ra,self.raunc,self.dec,self.decunc,self.para,self.para_over_err,self.pmra,self.pmraunc,self.pmdec,self.pmdecunc,self.excess_sig = data[:,1],data[:,2],data[:,3],data[:,4],data[:,5],data[:,6],data[:,7],data[:,8],data[:,9],data[:,10],data[:,11]
        
        for i in range(len(self.raunc)):
            self.raunc[i] = self.raunc[i]/3600000
            self.decunc[i] = self.decunc[i]/3600000
    #produces an array of indices with crossmatched stars, but with contradictory pm motions
    
    def crossmatch(self,odat):
        
        cross_result = []
        
        for i in range(len(odat.ra)):
            
            print(i)
            
            for j in range(len(self.ra)):
                
                if ((self.ra[j]-(self.raunc[j] + odat.raunc)<odat.ra[i]<self.ra[j]+(self.raunc[j]+odat.raunc)) and (self.dec[j]-(self.decunc[j] + odat.decunc)<odat.dec[i]<self.dec[j]+(self.decunc[j]+odat.decunc))):
                    print('Crossmatch at index' + str(i))
                    print(self.ra[j])
                    print(odat.ra[i])
                    print(self.raunc[j])
                                        
                    if self.pmra[j] == np.nan:
                        
                        cross_result.append(1)
                        
                        continue
                    
                    elif ((self.pmra[j]-(self.pmraunc[j] + odat.pmraunc)<odat.pmra<self.pmra[j]+(self.pmraunc[j]+odat.pmraunc)) and (self.pmdec[j]-(self.pmdecunc[j] + odat.pmdecunc)<odat.pmdec<self.pmdec[j]+(self.pmdecunc[j]+odat.pmdecunc))):
                        
                        cross_result.append(2)
                    
                    else:
                        
                        cross_result.append(3)
                        
                    
                else:
                    cross_result.append(0)
        
        return np.array(cross_result)
                    
    #0 - not crossmatched
    #1 - crossmatched, no pm available
    #2 - crossmatched, gaia pm within uncertainty of galactic pm
    #3 - crossmatched, pm not within uncertainty of galactic pm
    
    #retain 0, 1, 2, discard 3 as it shows a good measurement of pm, most likely foreground star
    
    #takes object being crossmatched with gaia data as input   
    
    #def crossmatch(galaxypms)
        
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

#n147 = asymphr('lot_n147.unique',0)
#kj_cmd(n147)
#colour_colour(n147)
#spatial_plot(n147)

#n185 = asymphr('lot_n185.unique',0)
#kj_cmd(n185)
#spatial_plot(n185)

#n205 = asymphr('N205_two.asc',0)
#kj_cmd(n205)
#spatial_plot(n205)

#m32 = asymphr('M32.asc',0)
#kj_cmd(m32)
#spatial_plot(m32)
    
n147 = asymphr('lot_n147.unique',0)
n147.loadascii()
n147.ciscuts()
gaian147 = gaiadata('ngc147gaiapm.csv')
gaian147.loadfile()
gaian147.crossmatch(n147)


#print(gaian147.data)

#standard coordinates
    

        
        
