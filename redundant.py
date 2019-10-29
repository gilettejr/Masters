#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 16:15:26 2019

@author: cameronrobertson
"""

import numpy as np


#fileformat: na,na,na,na,na,na,na,x,y,jmag,jerr,jcis,x,y,hmag,herr,hcis,x,y,kmag,kerr,kcis

#class for reading gaia data tables in format described below

#gaiaformat:sourceid,ra,raunc,dec,decunc,pmra,pmraunc,pmdec,pmdecunc
        
#class initialises by loading datafile into numpy array
#basically redundant due to topcat crossmatching, will still include as used to verify topcat crossmatch quality
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
    #produces an array of indices with crossmatched stars
    #
    #method to crossmatch gaia data with wfcam data based on ra and dec to within ~1arcsec
    def crossmatch(self,odat):
        
        cross_result = []
        count = []
        
        for i in range(len(odat.ra)):
            
            print(i)
            
            
            for j in range(len(self.ra)):
                
                if ((self.ra[j]-(self.raunc[j] + odat.raunc)<odat.ra[i]<self.ra[j]+(self.raunc[j]+odat.raunc)) and (self.dec[j]-(self.decunc[j] + odat.decunc)<odat.dec[i]<self.dec[j]+(self.decunc[j]+odat.decunc))):
                    print('Crossmatch at index' + str(i))
                    count.append(0)
                    cross_result.append(i)

                        


        print('Total number of crossmatches:')
        print(len(count))
        cross_result = np.array(cross_result)
        np.save('gaia_crossmatch_'+self.filename,cross_result)
    #0 - not crossmatched
    #1 - crossmatched, no pm available
    #2 - crossmatched, gaia pm within uncertainty of galactic pm
    #3 - crossmatched, pm not within uncertainty of galactic pm
    
    #retain 0, 1, 2, discard 3 as it shows a good measurement of pm, most likely foreground star
    
    #takes object being crossmatched with gaia data as input   
    
    #def crossmatch(galaxypms)
        
#function plots k-j colour magnitude diagram from asymphr class execution