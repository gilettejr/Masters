#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 15:14:03 2019

@author: cameronrobertson
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
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
            
            

        



#class for dealing with crossmatching operations and further data analysis based on crossmatched data
class topcatcross(asymphr):
    
    #creates array containing ra and dec in decimal coordinates for crossmatching. Saves to a csv file
    #must be run after self.loadascii
    
    
    
    def load_ascii_to_cross(self):
        

        #attributes set to variables for ease
        
        index=self.index
        ra=self.ra
        dec=self.dec


        #quantity table initialised
        
        t=QTable()
        
        #columns added for index, ra and dec in degrees
        
        t['index_no'] = index
        t['ra'] = ra*u.deg
        t['dec'] = dec*u.deg
        
        #qtable written and saved to a fits file
        
        t.write(self.filename+'crosscat.fits',overwrite=True)
        
    #method to simply read in and return csv file of crossed data from topcat
    #used in the crossed_data.topmatch() method
        
    def read_crossed_csv(self,crossfile):
        
        #class attribute added so that read file can be checked
        
        self.crossfile = crossfile
        
        #DataFrame object
        
        crossed = pd.read_csv(crossfile)
        
        #print(crossed.index_no)
        
        return crossed
    
    #function for deleting all crossmatched points
    
    def delete_crossed_points(self):
        
        #locates each crossed point using index identifiers from crossmatched gaia file
        
        for i in self.index_nos:
        
            #ignores NaN values
            
            if np.isnan(i)==True:
                
                continue
            
            
            else:
                
                #index numbers are floats, need to be converted to integers for identifying array elements
                
                i=int(i)
                
                #entire row of data wiped if crossmatched with gaia data
                
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
                self.eta[i] = np.nan
                self.xi[i]=np.nan
            
    #method for converting rows of data into pandas dataframe for ease
                
    def make_dataframe(self):
        
        #dataframe created
        
        wfdat=pd.DataFrame({'ra':self.ra,'dec':self.dec,'xj':self.xj,'yj':self.yj,'jmag':self.jmag,'jerr':self.jerr,'jcls':self.jcis,'xh':self.xh,'yh':self.yh,'hmag':self.hmag,'herr':self.herr,'hcls':self.hcis,'xk':self.xk,'yk':self.yk,'kmag':self.kmag,'kerr':self.kerr,'kcls':self.kcis,'xi':self.xi,'eta':self.eta})
        
        #attribute dataframe set
        
        self.wfdat=wfdat
    
    
    #returns dataframe only containing C star candidates based on hard colour cuts
    
    def select_C_stars(self,j_kcut,h_kcut,minkmag):
        
        #copy made to protect original dataframe
        
        d=self.wfdat.copy()
        
        #hard colour cuts recursively made
        
        for i in range(len(d.ra)):
            
            if d.jmag[i]-d.kmag[i] < j_kcut or d.hmag[i]-d.kmag[i] < h_kcut or d.kmag[i] > minkmag:
                
                #entire row wiped if conditions not met
                
                d.loc[i]=np.nan
                
        #altered dataframe copy returned        
                
        return d
    
    #returns dataframe only containing M star candidates based on hard colour cuts
    
    def select_M_stars(self,j_kmin,j_kmax,minkmag):
        
        #copy made to protect original dataframe
        
        d=self.wfdat.copy()
        
        #hard colour cuts recursively made
        
        for i in range(len(d.ra)):
            
            if j_kmin > d.jmag[i]-d.kmag[i] or j_kmax < d.jmag[i]-d.kmag[i] or d.kmag[i] > minkmag:
                
                #entire row wiped if conditions not met

                d.loc[i]=np.nan

        #altered dataframe copy returned
        
        return d
            
                        
        
#class containing methods for plotting various graphs from analysis output

class graphs:
    
        
    #method produces K-J CMD, taking in asymphr or asympher inherited objects as input    

    def kj_cmd(self,target):
    
    #figure created
        
        plt.figure()
                
        
        #plot formatted and labelled
        
        plt.rc('axes',labelsize = 20)
        plt.scatter(target.jmag-target.kmag,target.kmag,s=3,marker='o',alpha = 0.5)
        plt.gca().invert_yaxis()
        plt.ylabel('$K_0$')
        plt.xlabel('J-$K_0$')
        plt.title(target.filename + ' CMD')
        
        plt.show()
        
        #ngc147 cut
        
        #plt.axvline(x=1.34,linestyle=':',color='black')
        #plt.axhline(y=18,linestyle=':',color='black')
    

    #method plots spatial distribution from target class or dataframe

    def spatial_plot(self,target):
        
        #cuts performed, extinction correction not necessary
        
        
        #plot formatted and labelled
        
        plt.rc('axes',labelsize = 20)
        plt.scatter(target.dec,target.ra,s=3,marker = 'o',alpha = 0.5)
        plt.ylabel('RA(J200)/degrees')
        plt.xlabel('Dec(J200)/degrees')
        plt.title(target.filename + ' Spatial Plot')
        plt.show()
        
    #method plots raw camera xy pixel positions from target class or dataframe
        
    def xy_spatial_plot(self,target):

        
        plt.rc('axes',labelsize=20)
        plt.scatter(target.xj,target.yj,s=3,marker='o',alpha=0.5)
        plt.ylabel('Y')
        plt.xlabel('X')
        plt.title(target.filename + ' J-band pixel co-ordinates')
    
    #method plots spatial distribution of target class/dataframe using tangent co-ordinates
    
    def spatial_plot_standard(self,target,tangentra,tangentdec):
        
        #target.loadascii()
        #target.ciscuts()
        
        #co=ordinates converted fro degrees to radians
        
        ra = np.radians(target.ra)
        dec = np.radians(target.dec)
        tanra = np.radians(tangentra)
        tandec = np.radians(tangentdec)
        
        #tangent co-ordinates xi and eta constructed
        
        xi = (np.cos(dec)*np.sin(ra-tanra))/(np.sin(dec)*np.sin(tandec) + np.cos(dec)*np.cos(tandec)*np.cos(ra-tanra))
        
        eta = (np.sin(dec)*np.cos(tandec)-np.cos(dec)*np.sin(tandec)*np.cos(ra-tanra))/(np.sin(dec)*np.cos(tandec)+np.cos(dec)*np.sin(tandec)*np.cos(ra-tanra))
        
        #tangent co-ordinates converted back from radians into degrees
        
        xi = xi * (180/np.pi)
        eta = eta *(180/np.pi)
        
        plt.rc('axes',labelsize = 20)
        plt.scatter(xi,eta,s=3,marker = 'o',alpha=0.5)
        plt.gca().invert_xaxis()
        plt.gca().set_ylabel(r'$\eta$')
        plt.gca().set_xlabel(r'$\xi$')
        
    #method plots tangent point co-ordinate spatial distribution of two subsets of stars
    #process is simply duplicated version of method spatial_plot_standard for two subsets
        
    def agb_spatial_plot_standard(self,ctarget,mtarget,tangentra,tangentdec):
        
        #target.loadascii()
        #target.ciscuts()
        
        cra = np.radians(ctarget.ra)
        cdec = np.radians(ctarget.dec)
        tanra = np.radians(tangentra)
        tandec = np.radians(tangentdec)
        
        mra = np.radians(mtarget.ra)
        mdec = np.radians(mtarget.dec)

        
        cxi = (np.cos(cdec)*np.sin(cra-tanra))/(np.sin(cdec)*np.sin(tandec) + np.cos(cdec)*np.cos(tandec)*np.cos(cra-tanra))
        
        ceta = (np.sin(cdec)*np.cos(tandec)-np.cos(cdec)*np.sin(tandec)*np.cos(cra-tanra))/(np.sin(cdec)*np.cos(tandec)+np.cos(cdec)*np.sin(tandec)*np.cos(cra-tanra))
        
        mxi = (np.cos(mdec)*np.sin(mra-tanra))/(np.sin(mdec)*np.sin(tandec) + np.cos(mdec)*np.cos(tandec)*np.cos(mra-tanra))
        
        meta = (np.sin(mdec)*np.cos(tandec)-np.cos(mdec)*np.sin(tandec)*np.cos(mra-tanra))/(np.sin(mdec)*np.cos(tandec)+np.cos(mdec)*np.sin(tandec)*np.cos(mra-tanra))
        
        cxi = cxi * (180/np.pi)
        ceta = ceta *(180/np.pi)
        
        mxi = mxi * (180/np.pi)
        meta = meta *(180/np.pi)
        plt.style.use('seaborn-white')
        plt.rc('axes',labelsize = 20)
        plt.scatter(cxi,ceta,s=3,marker = 'o',alpha=0.5,label='C-Star Candidates')
        plt.scatter(mxi,meta,s=3,marker = 'o',alpha=0.5,label='M-Star Candidates')
        plt.gca().set_ylabel(r'$\eta$')
        plt.gca().set_xlabel(r'$\xi$')
        plt.gca().invert_xaxis()
        plt.legend()
        plt.show()
        
    #method for plotting single subset spatial distribution in tangent co-ordinates
    #can only be used once target.xi and target.eta attributes have already been set using the 
    #asymphr().make_tangent_coords method

    def single_agb_spatial_plot_standard(self,target):
        
        plt.style.use('seaborn-white')
        plt.rc('axes',labelsize=20)
        plt.scatter(target.xi,target.eta,marker='o',alpha=0.5,size=3)
        plt.gca().set_ylabel(r'$\eta$')
        plt.gca().set_xlabel(r'$\xi$')
        plt.gca().invert_xaxis()
        plt.show()
        
    #method for plotting colour-colour diagram from target object/subset

    def colour_colour(self,target):
        
        plt.figure()
        
        #colour cuts made, extinction corrected
        
        target.loadascii()
        target.ciscuts()
        target.sbsextinction()
        
        #plot formatted and labelled
        
        plt.rc('axes',labelsize = 20)
        #graph is in the format: (j-k)/(h-k)
        plt.scatter(target.hmag - target.kmag,target.jmag-target.kmag,s=1,marker=',',alpha = 0.5)
        plt.ylabel('(J-$K)_0$')
        plt.xlabel('(H-$K)_0$')
        plt.title(target.filename + 'Colour-colour diagram')
        
        #line for ngc147 hard colour C cut, Y-J sohn et al.
        
        plt.axvline(x=0.44,linestyle=':',color='black')
        plt.axhline(y=1.34,linestyle=':',color='black')
    
    
    
        plt.show()
    
    #redundant, replaced by plot_topmatch_cmd
    
    def plot_crossmatch_cmd(self,target,cross_indices):
        
        target.loadascii()
        target.ciscuts()
        target.sbsextinction()
        
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
    
    #plots crossmatched data from gaia DR2 with uncrossmatched data in CMD
    #seperates the two subsets and labels before graphig
    
    def plot_topmatch_cmd(self,cross):
                    
        plt.rc('axes',labelsize=20)
        
        plt.scatter(cross.jmag-cross.kmag,cross.kmag,s=3,marker='o',alpha=0.5,label='UKIRT Data')
        plt.scatter(cross.crossjmag-cross.crosskmag,cross.crosskmag,s=3,marker='o',alpha=1,label='UKIRT Data crossmatched with Gaia DR2 catalogue')
        
        
        
        plt.gca().invert_yaxis()
        plt.ylabel('$K_0$')
        plt.xlabel('$(J-K)_0$')
        plt.legend()
        
        plt.show()
    
    #plots cmd with overlaid isochrones from iso asc file as input to method
    #specifically for looking at TRGB region on CMD
    
    def isoplot(self,target,iso):
        
        #truncates isochrones to only show RGB region
        
        for i in range(len(iso.jmag)):
            
            #label=3 is the identifier for RGB phase
            
            if iso.label[i]!=3:

                iso.jmag[i]=np.nan
                iso.hmag[i]=np.nan
                iso.kmag[i]=np.nan
        
        plt.rc('axes',labelsize=20)
        
        #this section is written so that isochrones of populations with different
        #metallicities and/or ages or seperated and not treated as one continuous line
        
        #list created to keep track of age or metallicity changes in the isochrones
        indices=[]
        
        #loop must begin at i=1 since array element i-1 is used
        
        for i in range(1,len(iso.age)):
            
            #code searches for age or z change from one element to the next
            #appends the element number to indices to mark when this occurs
            
            if iso.age[i]!=iso.age[i-1] or iso.z[i]!=iso.z[i-1]:
                indices.append(i)
 
        #standard cmd of target object/subset is plotted
        
        plt.scatter(target.jmag-target.kmag,target.kmag,s=3,marker='o',alpha=0.5,label='UKIRT Data',color='black')

        #first isochrone must be plotted independently, as it is not included in subsequent plotting loop
        #plots isochrone data from start of file until the lowest value in indices, which is the 0th element
        
        plt.plot(iso.jmag[:indices[0]]-iso.kmag[:indices[0]],iso.kmag[:indices[0]],label='log(t) = ' + str(iso.age[0]) + ', Z = ' +str(iso.z[0]))
        
        #rest of the isochrones are plotted recursively, using the indices list to
        #define when a new isochrone is begun
        
        for i in range(len(indices)):
            
            #conditional added to treat the final isochrone and terminating loop
            #not including this results in list index out of bounds error
            
            if i==(len(indices)-1):
                
                #final isochrone plotted from highest indices element to the end of isochrone set file
                
                plt.plot(iso.jmag[indices[i]:]-iso.kmag[indices[i]:],iso.kmag[indices[i]:],label='log(t) = ' + str(iso.age[indices[i]]) + ', Z = ' +str(iso.z[indices[i]]))
                break
            else:
                
                #plots all isochrone lines between the first and last elements defined by indices
                
                plt.plot(iso.jmag[indices[i]:indices[i+1]]-iso.kmag[indices[i]:indices[i+1]],iso.kmag[indices[i]:indices[i+1]],label='log(t) = ' + str(iso.age[indices[i]]) + ', Z = ' +str(iso.z[indices[i]]))
        
        plt.gca().invert_yaxis()
        plt.ylabel('$K_0$')
        plt.xlabel('$(J-K)_0$')
        
        #ylim set since data is irrelevant above y=20
        
        plt.ylim(20,12)
        plt.legend()
        
        plt.show()
        
    #plots simple histogram based on j-k colour of two subset dataframes, most likely to be C and M stars    
        
    def colour_hist(self,cframe,mframe):
        
        #dropna() is crucial here as seaborn doesn't seem to like NaN values
        
        sns.distplot(cframe.jmag.dropna()-cframe.kmag.dropna(),label='C-Stars')
        sns.distplot(mframe.jmag.dropna()-mframe.kmag.dropna(),label='M-Stars')
        plt.legend()
        plt.show()
        
    #plots j-k colour histogram from one subset dataframe
    
    def single_colour_hist(self,frame):
        
        #again, dropna() required for plotting dataframe
        #is it really that hard for seaborn not to plot it???????

        #kde=True/False changes whether this plots a fitted and normalised
        #histogram or just a bog standard one
        
        sns.distplot(frame.jmag.dropna()-frame.kmag.dropna(),kde=False)
        plt.legend()
        plt.show()
    
    #plots agb section of isochrones
    #method is just for overlaying so can be used with another plot in this class
    
    def overlay_agb_tip(self,iso):
        
        #truncates isochrones to only include AGB
        
        for i in range(len(iso.jmag)):
            
            #label=8 is the AGB phase in padova isochrones
            
            if iso.label[i]!=8:

                iso.jmag[i]=np.nan
                iso.hmag[i]=np.nan
                iso.kmag[i]=np.nan
        
        #same method for seperating out the different isochrones in the set
        #as in self.isoplot
        
        indices=[]
        for i in range(1,len(iso.age)):
            if iso.age[i]!=iso.age[i-1] or iso.z[i]!=iso.z[i-1]:
                indices.append(i)
 
                
        
        

        
        plt.plot(iso.jmag[:indices[0]]-iso.kmag[:indices[0]],iso.kmag[:indices[0]],label='log(t) = ' + str(iso.age[0]) + ', Z = ' +str(iso.z[0]))
        
        for i in range(len(indices)):
            if i==(len(indices)-1):
                plt.plot(iso.jmag[indices[i]:]-iso.kmag[indices[i]:],iso.kmag[indices[i]:],label='log(t) = ' + str(iso.age[indices[i]]) + ', Z = ' +str(iso.z[indices[i]]))
                break
            else:
                plt.plot(iso.jmag[indices[i]:indices[i+1]]-iso.kmag[indices[i]:indices[i+1]],iso.kmag[indices[i]:indices[i+1]],label='log(t) = ' + str(iso.age[indices[i]]) + ', Z = ' +str(iso.z[indices[i]]))
        
        plt.gca().invert_yaxis()
        plt.ylabel('$K_0$')
        plt.xlabel('$(J-K)_0$')
        plt.legend()
        
        plt.show()
        
    #plots spatial distribution of c and m agb stars as a contour density diagram
        
    def surface_density_plot(self,cframe,mframe):
        
        
        

        
        sns.set_style("white")
        
        #m and c star dataframe tangent co-ordinates merged into combined 1d dataframes
        #also deleted NaN values in this process
        xi=np.concatenate((cframe.xi.dropna(),mframe.xi.dropna()),axis=0)
        eta=np.concatenate((cframe.eta.dropna(),mframe.eta.dropna()),axis=0)
        
        #contour density plotted
        
        sns.kdeplot(xi,eta,n_levels=30,gridsize=500)
    
    #plots spatial distribution of stars from frame dataframe
        
    def single_surface_density_plot(self,frame):
        
        
        sns.set_style("white")
        
        #contour density plot created, again deleting NaN values
        
        sns.kdeplot(frame.xi.dropna(), frame.eta.dropna(),n_levels=300,gridsize=500)

        plt.show()
    
    #plots K band luminosity function of two subsets, most likely C and M stars
    
    def k_luminosity_function(self,cframe,mframe):
        sns.set_style('white')
        
        #histograms of K band magnitudes plotted, again deleting NaN values
        
        sns.distplot(cframe.kmag.dropna(),label='C-Stars')
        sns.distplot(mframe.kmag.dropna(),label='M-Stars')
        #sns.set(xlabel='$K_0$',ylabel='Number of AGB stars')
        plt.legend()
        plt.show()
    
    #plots K band luminosity function of one subset
    
    def single_k_luminosity_function(self,frame):
        sns.set_style('dark')
        #plots K band magnitude histogram, deleting NaN values
        sns.distplot(frame.kmag.dropna())
        plt.show()

#function cuts gaia data and sets class attributes for crossmatched data

#class for dealing with identifying data crossmatched by topcat
        
class crossed_data:
    
    #method sets crossmatched data from incrossfile to gaiacross class attributes
    #gaiacross must be topcatcross object
    
    def topmatch(self,gaiacross,incrossfile):
        
        #method called to read in crossmatched file as dataframe
        
        crossed_table = gaiacross.read_crossed_csv(incrossfile)
        
        #noise cut based on Gaia DR2 values and recomendations from the Gaia project website
        
        for i in range(len(crossed_table.pmra)):
            if crossed_table.astrometric_excess_noise_sig[i] > 2:
                crossed_table.index_no[i]=np.nan
                
        #parallax cut to remove any badly  measured parallaxes
                
        for i in range(len(crossed_table.pmra)):
            if crossed_table.parallax_over_error[i] < 1:
                crossed_table.index_no[i]=np.nan
                
        #proper motion cut to remove any proper motions pointing to distant objects (i.e. pm~0)
                
        for i in range(len(crossed_table.pmra)):
            if (np.abs(crossed_table.pmra[i]/crossed_table.pmra_error[i])) < 1 or (np.abs(crossed_table.pmdec[i]/crossed_table.pmdec_error[i])) < 1:
                crossed_table.index_no[i]=np.nan

        
        #WFCAM data is loaded and the usual cuts made
        gaiacross.loadascii()
        gaiacross.ciscuts()
        
        #extinction is not carried out at this stage since cuts are to be made before
        #extinction corrections
        
        #lists for crossmatched photometric data created
        
        crossjmag = []
        crosshmag = []
        crosskmag = []
        
        #index_no column used to identify gaia crossmatched data
        
        for j in crossed_table.index_no:
            
            #Nan values skipped
            
            if np.isnan(j) == True:
                
                continue
            
            else:
                
                #index number elements are initially floats, must be converted to integers
                #for list indexing
                
                j=(int(j))
                
                #empty crossmatching lists filled with crossmatched data
                
                crossjmag.append(gaiacross.jmag[j])
                crosskmag.append(gaiacross.kmag[j])
                crosshmag.append(gaiacross.hmag[j])
        
        #lists converted to numpy arrays
        
        crossjmag = np.array(crossjmag)
        crosshmag = np.array(crosshmag)
        crosskmag = np.array(crosskmag)
        
        #arrays of crossmatched data set as gaiacross attributes
        
        gaiacross.crossjmag=crossjmag
        gaiacross.crosshmag=crosshmag
        gaiacross.crosskmag=crosskmag
        
        #index numbers also set as attribute for future troubleshooting
        
        gaiacross.index_nos=crossed_table.index_no


        

        

    






    
    


#class for running graphing functions
#mainly involves calling the relevant graphs() method after class initialisation
class basic_graphs:
    
    #galaxy defined upon initialisation
    
    def __init__(self,galaxy):
        
        #galaxy set as class attribute
        
        self.galaxy=galaxy
        
        #conditional statements used to select correct datafile for the chosen galaxy
        
        if galaxy=='ngc147':
        
            n147 = asymphr('lot_n147.unique',0)

            
        elif galaxy== 'ngc185':
            n147 = asymphr('lot_n185.unique',0)

            
        elif galaxy=='ngc205':
            n147 = asymphr('N205_two.asc',0)

            
        elif galaxy=='m32':
            n147=topcatcross('M32.asc',0)
            
        elif galaxy=='andromeda':
            n147=topcatcross('lot_m31.unique',0)
            
        #error returned if galaxy not recognised 
        
        else:
            print('That is not a galaxy, enjoy the errors')
        
        #data loaded and cls cuts and extinction corrections carried out 
        
        n147.loadascii()
        n147.ciscuts()
        n147.sbsextinction()
        
        #data set as class attribute for running class methods
        
        self.n147=n147
        self.plotter=graphs()
    #uses the graphs().kj_cmd method to plot a CMD using class attribute data
        
    def plot_cmd(self):
        self.plotter.kj_cmd(self.n147)
    
    #uses the graphs().colour_colour method to plot a cc diagram using class attribute data
    
    def plot_cc(self):

        self.plotter.colour_colour(self.n147)
        
    #uses the graphs().spatial_plot method to plot a spatial distribution using class attribute data
    
    def plot_spatial(self):

        self.plotter.spatial_plot(self.n147)
        
    #uses the graphs().xy_spatial_plot method to plot an xy pixel distribution using class attribute data    
        
    def plot_xy_spatial(self):
        self.plotter.xy_spatial_plot(self.n147)
        
#class for setting up data before topcat crossmatching

class topcatstuff:
    
    #uses topcatcross class methods to convert ra and dec in named
    #file to decimal degrees, and converts the file from .csv to .fits
    
    def fit_for_cross(self,file):
        n147cross=topcatcross(file,0)
        n147cross.loadascii()
        n147cross.ciscuts()
        n147cross.load_ascii_to_cross()
        
    
#class containing methods for separating out classes of star depending on
#phase or other user defined cuts. Automatically removes stars crossmatched from
#gaia DR2
class make_subsets:
    
    #class initialised with defined galaxy. Deletes gaia crossmatched data and
    #converts data into dataframe format, as well as adding tangent co-ordinates
    #still a work in progress, only works for ngc147 right now!
    
    def __init__(self,galaxy):
        
        self.galaxy=galaxy
        
        rmatch=crossed_data()
        
        if galaxy=='ngc147':
        
            n147cross = topcatcross('lot_n147.unique',0)
            rmatch.topmatch(n147cross,'crossedn147.csv')
            
        elif galaxy== 'ngc185':
            n147cross = topcatcross('lot_n185.unique',0)
            rmatch.topmatch(n147cross,'crossedn185.csv')
            
        elif galaxy=='ngc205':
            n147cross=topcatcross('N205_two.asc',0)
            rmatch.topmatch(n147cross,'crossedn205.csv')
            
        elif galaxy=='m32':
            n147cross=topcatcross('M32.asc',0)
            rmatch.topmatch(n147cross,'crossedm32.csv')
        
        
            
        
        n147cross.delete_crossed_points()
        n147cross.create_tangent_coords(8.300500,48.50850)
        n147cross.make_dataframe()
        
        self.n147cross=n147cross
    
    def mc(self,agb):
        
        n147cross=self.n147cross
        
        if agb=='m':
        
        
            selection=n147cross.select_M_stars(1.0,1.34,18)
    
        elif agb=='c':
            
            selection=n147cross.select_C_stars(1.34,0.44,18)
    
            
            
        else:
            print('Not a valid class')
        
        n147cross.sbsextinction_on_frame(selection)
        self.subset=selection
        
    #def rgb(self,xrange,yrange):
        

  
class run_both:

    def __init__(self,galaxy):
        runm=make_subsets(galaxy)
        runc=make_subsets(galaxy)
        runm.mc('m')
        runc.mc('c')
        self.mframe=runm.subset
        self.cframe=runc.subset
        self.galaxy=galaxy
    def plot_both_spatial(self):
        plotter=graphs()
        
        galaxy=self.galaxy
        
        if galaxy=='ngc147':
            tra=8.3005
            tdec=48.50873889
        elif galaxy=='ngc185':
            tra=9.74154167
            tdec=48.33737778
        
        elif galaxy=='ngc205':
            tra=10.09189356
            tdec=41.68541564
            
        elif galaxy=='m32':
            tra=10.67427083
            tdec=40.86516944
        else:
            print('Not a valid object')
        
        
        plotter.agb_spatial_plot_standard(self.cframe,self.mframe,tra,tdec)
    
    def plot_both_lum(self):
        plotter=graphs()
        plotter.k_luminosity_function(self.cframe,self.mframe)
    
    def plot_both_contour(self):
        plotter=graphs()
        plotter.surface_density_plot(self.cframe,self.mframe)

    def plot_both_colour_hist(self):
        plotter=graphs()
        plotter.colour_hist(self.cframe,self.mframe)
        
    def c_over_m(self):
        m_no=[]
        for i in self.mframe.kmag:
            if np.isnan(i)==False:
                m_no.append(0)
        c_no=[]   
        for j in self.cframe.kmag:
            if np.isnan(j)==False:
                c_no.append(0)
        
        C=len(c_no)
        M=len(m_no)
        
        ratio = C/M
        
        print(ratio)
        
    def c_over_m_grad(self,border,tra,tdec):
        
        in_m_no=[]
        for i in range(len(self.mframe.kmag)):
            if np.isnan(self.mframe.ra[i])==False and np.sqrt((self.mframe.ra[i]-tra)**2+(self.mframe.dec[i]-tdec)**2) < border/3600 :
                in_m_no.append(0)
        in_c_no=[]   
        for i in range(len(self.mframe.kmag)):
            if np.isnan(self.cframe.ra[i])==False and np.sqrt((self.cframe.ra[i]-tra)**2+(self.cframe.dec[i]-tdec)**2) < border/3600 :
                in_c_no.append(0)
        
        in_C=len(in_c_no)
        in_M=len(in_m_no)
        
        out_m_no=[]
        for i in range(len(self.mframe.kmag)):
            if np.isnan(self.mframe.ra[i])==False and np.sqrt((self.mframe.ra[i]-tra)**2+(self.mframe.dec[i]-tdec)**2) > border/3600 :
                out_m_no.append(0)
        out_c_no=[]  
        
        for i in range(len(self.mframe.kmag)):
            if np.isnan(self.cframe.ra[i])==False and np.sqrt((self.cframe.ra[i]-tra)**2+(self.cframe.dec[i]-tdec)**2) > border/3600 :
                out_c_no.append(0)
        
        out_C=len(out_c_no)
        out_M=len(out_m_no)
        
        self.inCM=(in_C/in_M)
        self.outCM=(out_C/out_M)
        
        print(self.inCM)
        print(self.outCM)
        
    def CM_to_FEH(self,CM):
        
        FEH=-1.39 -0.47*np.log10(CM)
        
        return FEH
        
        
        
        

        
    
class run_single:
    
    def __init__(self,galaxy,subset):
        if subset !='m' and subset !='c':
            print('run_single argument must be string(m) or string(c)')
        run=make_subsets(subset,galaxy)
        
        self.frame=run.subset
        
    def plot_single_spatial(self):
        plotter=graphs()
        plotter.single_agb_spatial_plot_standard(self.frame)
        
    def plot_single_lum(self):
        plotter=graphs()
        plotter.single_k_luminosity_function(self.frame)
        
        
    def plot_single_contour(self):
        plotter=graphs()
        plotter.single_surface_density_plot(self.frame)
        
    def plot_single_colour_hist(self):
        plotter=graphs()
        plotter.single_colour_hist(self.frame)
        

            
        
class import_isos(basic_graphs):
    
    
    def read_in(self,isofile):
        
        def apparent(M,dguess):
            m = M + 5*np.log10(dguess/10)
            return m
        
        if self.galaxy=='ngc147':
            self.distance=780000
        
        distance=self.distance
        
        i=asymphr(isofile,0)
        
        isos=pd.DataFrame({'z':i.data[:,0],'age':i.data[:,2],'label':i.data[:,9],'jmag':i.data[:,27],'hmag':i.data[:,28],'kmag':i.data[:,29]})
        
        isos.jmag=apparent(isos.jmag,distance)
        isos.hmag=apparent(isos.hmag,distance)
        isos.kmag=apparent(isos.kmag,distance)
        
        i.z=isos.z
        i.age=isos.age
        i.label=isos.label
        i.jmag=isos.jmag
        i.hmag=isos.hmag
        i.kmag=isos.kmag
        
        self.i=i
        #plotter=graphs()
        #plotter.isoplot(self.n147,i)
        
class run_isos:
    
    def __init__(self,galaxy,isofile):
        imp=import_isos(galaxy)
        imp.read_in(isofile)
        self.imp=imp
        
    def plot_iso_cmd(self):
        plotter=graphs()
        plotter.isoplot(self.imp.n147,self.imp.i)
        
    def plot_iso_overlay(self):
        plotter=graphs()
        plotter.overlay_agb_tip(self.imp.i)
        

        

        
def main():  
        
    r=basic_graphs('ngc147')
    r.plot_cc()


main()


#radial gradiant: C/M = 0.22 <70'', C/M = 0.26 >70''