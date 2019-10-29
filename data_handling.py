#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 15:14:03 2019

@author: cameronrobertson
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy import units as u
from astropy.io import fits
from astropy.table import Table,Column,QTable
from astropy.coordinates import SkyCoord
from dustmaps.sfd import SFDQuery
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
            
            
    def sbsextinction(self):
        ra = self.ra
        dec = self.dec
            
        coords=SkyCoord(ra,dec,unit='deg',frame='icrs')
        sfd=SFDQuery()
        sfdred=sfd(coords)
           
        jext = sfdred * 0.709
        hext = sfdred * 0.449
        kext = sfdred * 0.302
        
        self.jmag=self.jmag - jext
        self.hmag=self.hmag - hext
        self.kmag=self.kmag - kext
            
            
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

class topcatcross(asymphr):
    
    #creates array containing ra and dec in decimal coordinates for crossmatching. Saves to a csv file
    #must be run after loadascii
    
    
    
    def load_ascii_to_cross(self):
        

        index=self.index
        ra=self.ra
        dec=self.dec
        

            
        
        crosscat = np.array([ra,dec])
        
        crosscat = crosscat.transpose()
        
        #np.savetxt('n147wfcam.csv',crosscat,delimiter=',')
        
        t=QTable()
        t['index_no'] = index
        t['ra'] = ra*u.deg
        t['dec'] = dec*u.deg
        
        t.write(self.filename+'crosscat.fits',overwrite=True)
        
    def read_crossed_csv(self,crossfile):
        
        
        self.crossfile = crossfile
        
        crossed = pd.read_csv(crossfile)
        
        print(crossed.index_no)
        
        return crossed
        

        
        
        

def kj_cmd(target):

#figure created
    
    plt.figure()
    
    #class methods carried out to load data and perform cuts and extinction corrections
    
    target.loadascii()
    target.ciscuts()
    target.sbsextinction()
    
    #plot formatted and labelled
    
    plt.rc('axes',labelsize = 20)
    plt.scatter(target.jmag-target.kmag,target.kmag,s=3,marker='o',alpha = 0.5)
    plt.gca().invert_yaxis()
    plt.ylabel('$K_0')
    plt.xlabel('J-$K_0$')
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
    
def spatial_plot_standard(target,tangentra,tangentdec):
    
    target.loadascii()
    target.ciscuts()
    
    ra = np.radians(target.ra)
    dec = np.radians(target.dec)
    tanra = np.radians(tangentra)
    tandec = np.radians(tangentdec)
    
    xi = (np.cos(dec)*np.sin(ra-tanra))/(np.sin(dec)*np.sin(tandec) + np.cos(dec)*np.cos(tandec)*np.cos(ra-tanra))
    
    eta = (np.sin(dec)*np.cos(tandec)-np.cos(dec)*np.sin(tandec)*np.cos(ra-tanra))/(np.sin(dec)*np.cos(tandec)+np.cos(dec)*np.sin(tandec)*np.cos(ra-tanra))
    
    plt.rc('axes',labelsize = 20)
    plt.scatter(xi,eta,s=3,marker = 'o',alpha=0.5)
    plt.gca().set_ylabel(r'$\eta$')
    plt.gca().set_xlabel(r'$\xi$')

def colour_colour(target):
    
    plt.figure()
    
    #colour cuts made, extinction corrected
    
    target.loadascii()
    target.ciscuts()
    target.extinction()
    
    #plot formatted and labelled
    
    plt.rc('axes',labelsize = 20)
    plt.scatter(target.jmag - target.hmag,target.hmag-target.kmag,s=1,marker=',',alpha = 0.5)
    plt.ylabel('H-$K_0$')
    plt.xlabel('J-$K_0$')
    plt.title(target.filename + 'Colour-colour diagram')
    
    plt.show()
    
def plot_crossmatch_cmd(target,cross_indices):
    
    target.loadascii()
    target.ciscuts()
    target.extinction()
    
    crossjmag = []
    crosshmag = []
    crosskmag = []
    for j in cross_indices:
        crossjmag.append(target.jmag[j])
        crosskmag.append(target.kmag[j])
        crosshmag.append(target.hmag[j])
        
    crossjmag=np.array(crossjmag)
    crosshmag = np.array(crosshmag)
    crosskmag = np.array(crosskmag)    
                
    plt.rc('axes',labelsize=20)
    
    plt.scatter(target.jmag-target.kmag,target.kmag,s=3,marker='o',alpha=0.5,label='UKIRT Data')
    plt.scatter(crossjmag-crosskmag,crosskmag,s=3,marker='o',alpha=1,label='UKIRT Data crossmatched with Gaia DR2 catalogue')
    
    plt.gca().invert_yaxis()
    plt.ylabel('$K_0$')
    plt.xlabel('$(J-K)_0$')
    plt.legend()
    
    plt.show()

def plot_topmatch_cmd(gaiacross,incrossfile):
    crossed_table = gaiacross.read_crossed_csv(incrossfile)
    
    #noise cut
    
    for i in range(len(crossed_table.pmra)):
        if crossed_table.astrometric_excess_noise_sig[i] > 2:
            crossed_table.index_no[i]=np.nan
            
    #parallax cut
            
    #for i in range(len(crossed_table.pmra)):
        #if crossed_table.parallax_over_error[i] < 1:
            ##crossed_table.index_no[i]=np.nan
            
    #proper motion cut
            
    for i in range(len(crossed_table.pmra)):
        if (np.abs(crossed_table.pmra[i]/crossed_table.pmra_error[i])) < 1 or (np.abs(crossed_table.pmdec[i]/crossed_table.pmdec_error[i])) < 1:
            crossed_table.index_no[i]=np.nan
    print(crossed_table.pmra)
    
    gaiacross.loadascii()
    gaiacross.ciscuts()
    gaiacross.extinction()
    
    crossjmag = []
    crosshmag = []
    crosskmag = []
    
    print(crossed_table.index_no)
    
    for j in crossed_table.index_no:
        
        if np.isnan(j) == True:
            
            continue
        
        else:
            crossjmag.append(gaiacross.jmag[int(j)])
            crosskmag.append(gaiacross.kmag[int(j)])
            crosshmag.append(gaiacross.hmag[int(j)])
        
    crossjmag = np.array(crossjmag)
    crosshmag = np.array(crosshmag)
    crosskmag = np.array(crosskmag)    
                
    plt.rc('axes',labelsize=20)
    
    plt.scatter(gaiacross.jmag-gaiacross.kmag,gaiacross.kmag,s=3,marker='o',alpha=0.5,label='UKIRT Data')
    plt.scatter(crossjmag-crosskmag,crosskmag,s=3,marker='o',alpha=1,label='UKIRT Data crossmatched with Gaia DR2 catalogue')
    
    plt.gca().invert_yaxis()
    plt.ylabel('$K_0$')
    plt.xlabel('$(J-K)_0$')
    plt.legend()
    
    plt.show()
    
    
    


#statements for running graphing functions. Galaxy chosen by initialising class before execution of functions

n147 = asymphr('lot_n147.unique',0)
kj_cmd(n147)
#colour_colour(n147)
#spatial_plot(n147)
##spatial_plot_standard(n147,8.300500,48.508750)

#n185 = asymphr('lot_n185.unique',0)
#kj_cmd(n185)
#spatial_plot(n185)
#spatial_plot_standard(n185,9.741542,48.337389)

#n205 = asymphr('N205_two.asc',0)
#kj_cmd(n205)
#spatial_plot(n205)
#spatial_plot_standard(n205,10.092000,41.685306)

#m32 = asymphr('M32.asc',0)
#kj_cmd(m32)
#spatial_plot(m32)
#spatial_plot_standard(m32,10.674300,40.865287)
    
#n147 = asymphr('lot_n147.unique',0)
#n147.loadascii()
#n147.ciscuts()


#gaian147 = gaiadata('ngc147gaiapm.csv')
#gaian147.loadfile()
#gaian147.crossmatch(n147)


##n147cross_cat = np.load('gaia_crossmatch_all_cull.npy')
#plot_crossmatch_cmd(n147,n147cross_cat)

#n185cross_cat = np.load('gaia_crossmatch_ngc185gaiapm.csv.npy')
#plot_crossmatch_cmd(n185,n185cross_cat)

#m32cross_cat = np.load('gaia_crossmatch_m32gaiapm.csv.npy')
#plot_crossmatch_cmd(m32,m32cross_cat)

#n205cross_cat = np.load('gaia_crossmatch_ngc205gaiapm.csv.npy')
#plot_crossmatch_cmd(n205,n205cross_cat)

#print(gaian147.data)

#standard coordinates
    
#n147cross=topcatcross('lot_n147.unique',0)
#n147cross.loadascii()
#n147cross.ciscuts()
#n147cross.load_ascii_to_cross()



#n185.loadascii()
#n185.ciscuts()
#gaian185 = gaiadata('ngc185gaiapm.csv')
#gaian185.loadfile()
#gaian185.crossmatch(n185)

#m32.loadascii()
#m32.ciscuts()
#gaiam32 = gaiadata('m32gaiapm.csv')
#gaiam32.loadfile()
#gaiam32.crossmatch(m32)

#n205.loadascii()
##n205.ciscuts()
#gaian205 = gaiadata('ngc205gaiapm.csv')
#gaian205.loadfile()
#gaian205.crossmatch(n205)

#n147cross = topcatcross('lot_n147.unique',0)
#plot_topmatch_cmd(n147cross,'crossedn147.csv')


