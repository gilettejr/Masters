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
        
    def sbsextinction_on_frame(self,target):
        
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
        
        
    def create_tangent_coords(self,tangentra,tangentdec):
        
        ra = np.radians(self.ra)
        dec = np.radians(self.dec)
        tanra = np.radians(tangentra)
        tandec = np.radians(tangentdec)
        
        xi = (np.cos(dec)*np.sin(ra-tanra))/(np.sin(dec)*np.sin(tandec) + np.cos(dec)*np.cos(tandec)*np.cos(ra-tanra))
        
        eta = (np.sin(dec)*np.cos(tandec)-np.cos(dec)*np.sin(tandec)*np.cos(ra-tanra))/(np.sin(dec)*np.cos(tandec)+np.cos(dec)*np.sin(tandec)*np.cos(ra-tanra))
        
        self.xi = xi * (180/np.pi)
        self.eta = eta * (180/np.pi)
            
            
        #print(np.shape((data[:,7:])))
        #print(np.shape(data))
        #xj = data[:,7]
        #print(xj)
        
#fileformat: na,na,na,na,na,na,na,x,y,jmag,jerr,jcis,x,y,hmag,herr,hcis,x,y,kmag,kerr,kcis

#class for reading gaia data tables in format described below

#gaiaformat:sourceid,ra,raunc,dec,decunc,pmra,pmraunc,pmdec,pmdecunc
        
#class initialises by loading datafile into numpy array
#basically redundant due to topcat crossmatching, will still include as used to verify topcat crossmatch quality



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
        
    #method to deal with csv file from topcat crossmatching
    #pretty much redundant as topmatch function does much the same thing
        
    def read_crossed_csv(self,crossfile):
        
        #class attribute added so that read file can be checked
        
        self.crossfile = crossfile
        
        #DataFrame object
        
        crossed = pd.read_csv(crossfile)
        
        #print(crossed.index_no)
        
        return crossed
    
    #function for deleting all crossmatched points
    
    def delete_crossed_points(self):
        for i in self.index_nos:
            
            if np.isnan(i)==True:
                
                continue
            
            
            else:
                
                i=int(i)
                

                
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
                
                
    def make_dataframe(self):
        wfdat=pd.DataFrame({'ra':self.ra,'dec':self.dec,'xj':self.xj,'yj':self.yj,'jmag':self.jmag,'jerr':self.jerr,'jcls':self.jcis,'xh':self.xh,'yh':self.yh,'hmag':self.hmag,'herr':self.herr,'hcls':self.hcis,'xk':self.xk,'yk':self.yk,'kmag':self.kmag,'kerr':self.kerr,'kcls':self.kcis,'xi':self.xi,'eta':self.eta})
        self.wfdat=wfdat
    
    
    #returns dataframe only containing C star candidates based on hard colour cuts
    def select_C_stars(self,j_kcut,h_kcut,minkmag):
        d=self.wfdat.copy()
        
        
        for i in range(len(d.ra)):
            
            if d.jmag[i]-d.kmag[i] < j_kcut or d.hmag[i]-d.kmag[i] < h_kcut or d.kmag[i] > minkmag:
                d.loc[i]=np.nan
                
                
        return d
                
    def select_M_stars(self,j_kmin,j_kmax,minkmag):
        
        d=self.wfdat.copy()
        
        for i in range(len(d.ra)):
            
            if j_kmin > d.jmag[i]-d.kmag[i] or j_kmax < d.jmag[i]-d.kmag[i] or d.kmag[i] > minkmag:

                d.loc[i]=np.nan

                
        return d
            
                        
        
            

class graphs:
    
        
        

    def kj_cmd(self,target):
    
    #figure created
        
        plt.figure()
        
        #class methods carried out to load data and perform cuts and extinction corrections
        
        #target.loadascii()
        ##target.ciscuts()
        #target.sbsextinction()
        
        #plot formatted and labelled
        
        plt.rc('axes',labelsize = 20)
        plt.scatter(target.jmag-target.kmag,target.kmag,s=3,marker='o',alpha = 0.5)
        plt.gca().invert_yaxis()
        plt.ylabel('$K_0$')
        plt.xlabel('J-$K_0$')
        plt.title(target.filename + ' CMD')
        
        plt.show()
        
        #ngc147 cut
        
        plt.axvline(x=1.34,linestyle=':',color='black')
        plt.axhline(y=18,linestyle=':',color='black')
    

    #function plots spatial distribution from asymphr class execution

    def spatial_plot(self,target):
        
        plt.figure()
        
        #cuts performed, extinction correction not necessary
        
        #target.loadascii()
        #target.ciscuts()
        
        #plot formatted and labelled
        
        plt.rc('axes',labelsize = 20)
        plt.scatter(target.dec,target.ra,s=3,marker = 'o',alpha = 0.5)
        plt.ylabel('RA(J200)/degrees')
        plt.xlabel('Dec(J200)/degrees')
        plt.title(target.filename + ' Spatial Plot')
        plt.show()
    
    def spatial_plot_standard(self,target,tangentra,tangentdec):
        
        #target.loadascii()
        #target.ciscuts()
        
        ra = np.radians(target.ra)
        dec = np.radians(target.dec)
        tanra = np.radians(tangentra)
        tandec = np.radians(tangentdec)
        
        xi = (np.cos(dec)*np.sin(ra-tanra))/(np.sin(dec)*np.sin(tandec) + np.cos(dec)*np.cos(tandec)*np.cos(ra-tanra))
        
        eta = (np.sin(dec)*np.cos(tandec)-np.cos(dec)*np.sin(tandec)*np.cos(ra-tanra))/(np.sin(dec)*np.cos(tandec)+np.cos(dec)*np.sin(tandec)*np.cos(ra-tanra))
        
        xi = xi * (180/np.pi)
        eta = eta *(180/np.pi)
        
        plt.rc('axes',labelsize = 20)
        plt.scatter(xi,eta,s=3,marker = 'o',alpha=0.5)
        plt.gca().set_ylabel(r'$\eta$')
        plt.gca().set_xlabel(r'$\xi$')
        
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
        plt.legend()
        plt.show()

    def single_agb_spatial_plot_standard(self,target):
        
        plt.style.use('seaborn-white')
        plt.rc('axes',labelsize=20)
        plt.scatter(target.xi,target.eta,marker='o',alpha=0.5,size=3)
        plt.gca().set_ylabel(r'$\eta$')
        plt.gca().set_xlabel(r'$\xi$')
        plt.gca().invert_xaxis()
        plt.show()
        
        

    def colour_colour(self,target):
        
        plt.figure()
        
        #colour cuts made, extinction corrected
        
        target.loadascii()
        target.ciscuts()
        target.sbsextinction()
        
        #plot formatted and labelled
        
        plt.rc('axes',labelsize = 20)
        plt.scatter(target.hmag - target.kmag,target.jmag-target.kmag,s=1,marker=',',alpha = 0.5)
        plt.ylabel('(J-$K)_0$')
        plt.xlabel('(H-$K)_0$')
        plt.title(target.filename + 'Colour-colour diagram')
        
        #line for ngc147 hard C cut, Y-J sohn et al.
        
        plt.axvline(x=0.44,linestyle=':',color='black')
        plt.axhline(y=1.34,linestyle=':',color='black')
    
    
    
        plt.show()
    
    #redundant, replaces with plot_topmatch_cmd
    
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
        
    def plot_topmatch_cmd(self,cross):
                    
        plt.rc('axes',labelsize=20)
        
        plt.scatter(cross.jmag-cross.kmag,cross.kmag,s=3,marker='o',alpha=0.5,label='UKIRT Data')
        plt.scatter(cross.crossjmag-cross.crosskmag,cross.crosskmag,s=3,marker='o',alpha=1,label='UKIRT Data crossmatched with Gaia DR2 catalogue')
        
        
        
        plt.gca().invert_yaxis()
        plt.ylabel('$K_0$')
        plt.xlabel('$(J-K)_0$')
        plt.legend()
        
        plt.show()
        
    def isoplot(self,target,iso):
        
        for i in range(len(iso.jmag)):
            if iso.label[i]!=3:

#m32.loadascii()
#m32.ciscuts()
#gaiam32 = gaiadata('m32gaiapm.csv')
#gaiam32.loadfile()
#gaiam32.crossmatch(m32)

#n205.loadascii()
##n205.ciscuts()
#gaian205 = gaiadata('ngc205gaiapm.csv')
                iso.jmag[i]=np.nan
                iso.hmag[i]=np.nan
                iso.kmag[i]=np.nan
        
        plt.rc('axes',labelsize=20)
        
        plt.scatter(target.jmag-target.kmag,target.kmag,s=3,marker='o',alpha=0.5,label='UKIRT Data')
        plt.scatter(iso.jmag-iso.kmag,iso.kmag,s=1,marker='o',label='Padova isochrones')
        
        plt.gca().invert_yaxis()
        plt.ylabel('$K_0$')
        plt.xlabel('$(J-K)_0$')
        plt.legend()
        
        plt.show()
        
    def surface_density_plot(self,cframe,mframe):
        
        #xi=np.concatenate((cframe.xi.dropna(),mframe.xi.dropna()),axis=0)
        #eta=np.concatenate((cframe.eta.dropna(),mframe.eta.dropna()),axis=0)
        
        sns.set_style("white")
        sns.kdeplot(mframe.xi.dropna(), mframe.eta.dropna(),n_levels=300,gridsize=500)
        #sns.kdeplot(cframe.xi.dropna(),cframe.eta.dropna(),n_levels=30,gridsize=500)
        #sns.kdeplot(xi,eta,n_levels=30,gridsize=500)
        
    def single_surface_density_plot(self,frame):
        
        #xi=np.concatenate((cframe.xi.dropna(),mframe.xi.dropna()),axis=0)
        #eta=np.concatenate((cframe.eta.dropna(),mframe.eta.dropna()),axis=0)
        
        sns.set_style("white")
        sns.kdeplot(frame.xi.dropna(), frame.eta.dropna(),n_levels=300,gridsize=500)

        plt.show()
        
    def k_luminosity_function(self,cframe,mframe):
        sns.set_style('white')
        sns.distplot(cframe.kmag.dropna(),label='C-Stars')
        sns.distplot(mframe.kmag.dropna(),label='M-Stars')
        #sns.set(xlabel='$K_0$',ylabel='Number of AGB stars')
        plt.legend()
        plt.show()
        
    def single_k_luminosity_function(self,frame):
        sns.set_style('dark')
        sns.distplot(frame.kmag.dropna())
        plt.show()

#function cuts gaia data and sets class attributes for crossmatched data

class isochrone_plots(asymphr):
    
    def readisochrones(self):
        
        def apparent(M,dguess):
            m = M + 5*np.log10(dguess/10)
            return m
        
        
            
            
        

        distance=self.distance
        
        isos=pd.DataFrame({'label':self.data[:,9],'jmag':self.data[:,30],'hmag':self.data[:,31],'kmag':self.data[:,32]})
        
        isos.jmag=apparent(isos.jmag,distance)
        isos.hmag=apparent(isos.hmag,distance)
        isos.kmag=apparent(isos.kmag,distance)
        
        self.label=isos.label
        self.jmag=isos.jmag
        self.hmag=isos.hmag
        self.kmag=isos.kmag
        
        
class crossed_data:

    def topmatch(self,gaiacross,incrossfile):
        crossed_table = gaiacross.read_crossed_csv(incrossfile)
        
        #noise cut
        
        for i in range(len(crossed_table.pmra)):
            if crossed_table.astrometric_excess_noise_sig[i] > 2:
                crossed_table.index_no[i]=np.nan
                
        #parallax cut
                
        for i in range(len(crossed_table.pmra)):
            if crossed_table.parallax_over_error[i] < 1:
                crossed_table.index_no[i]=np.nan
                
        #proper motion cut
                
        for i in range(len(crossed_table.pmra)):
            if (np.abs(crossed_table.pmra[i]/crossed_table.pmra_error[i])) < 1 or (np.abs(crossed_table.pmdec[i]/crossed_table.pmdec_error[i])) < 1:
                crossed_table.index_no[i]=np.nan
        #print(crossed_table.pmra)
        
        gaiacross.loadascii()
        gaiacross.ciscuts()
        #gaiacross.sbsextinction()
        
        crossjmag = []
        crosshmag = []
        crosskmag = []
        
        #print(crossed_table.index_no)
        
        for j in crossed_table.index_no:
            
            if np.isnan(j) == True:
                
                continue
            
            else:
                j=(int(j))
                crossjmag.append(gaiacross.jmag[j])
                crosskmag.append(gaiacross.kmag[j])
                crosshmag.append(gaiacross.hmag[j])
            
        crossjmag = np.array(crossjmag)
        crosshmag = np.array(crosshmag)
        crosskmag = np.array(crosskmag)
        
        gaiacross.crossjmag=crossjmag
        gaiacross.crosshmag=crosshmag
        gaiacross.crosskmag=crosskmag
        
    
        
        gaiacross.index_nos=crossed_table.index_no


        

        

    






    
    


#statements for running graphing functions. Galaxy chosen by initialising class before execution of functions
class basic_graphs:
    def __init__(self,galaxy):
        if galaxy=='ngc147':
        
            n147 = asymphr('lot_n147.unique',0)

            
        elif galaxy== 'ngc185':
            n147 = asymphr('lot_n185.unique',0)

            
        elif galaxy=='ngc205':
            n147 = asymphr('N205_two.asc',0)

            
        elif galaxy=='m32':
            n147=topcatcross('M32.asc',0)
        n147.loadascii()
        n147.ciscuts()
        n147.sbsextinction()
        self.n147=n147
        
    def plot_cmd(self):
        plotter=graphs()
        plotter.kj_cmd(self.n147)
        
    def plot_cc(self):
        plotter=graphs()
        plotter.colour_colour(self.n147)
    
    def plot_spatial(self):
        plotter=graphs()
        plotter.spatial_plot(self.n147)
        

class topcatstuff:
    
    def fit_for_cross(self):
        n147cross=topcatcross('lot_n147.unique',0)
        n147cross.loadascii()
        n147cross.ciscuts()
        n147cross.load_ascii_to_cross()
        
    
#run these to remove crossmatched stars and separate out m and c stars
class make_subsets:
    def __init__(self,agb,galaxy):
        

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
            n147cross==topcatcross('M32.asc',0)
            rmatch.topmatch(n147cross,'crossedm32.csv')
        
        
            
        
        n147cross.delete_crossed_points()
        n147cross.create_tangent_coords(8.300500,48.50850)
        n147cross.make_dataframe()
        
        if agb=='m':
        
        
            selection=n147cross.select_M_stars(1.0,1.34,18)
    
        elif agb=='c':
            
            selection=n147cross.select_C_stars(1.34,0.44,18)
    
            
            
        else:
            print('Not a valid class')
        
        n147cross.sbsextinction_on_frame(selection)
        self.subset=selection
        

  
class run_both:

    def __init__(self,galaxy):
        runm=make_subsets('m',galaxy)
        runc=make_subsets('c',galaxy)
        self.mframe=runm.subset
        self.cframe=runc.subset
      
    def plot_both_spatial(self):
        plotter=graphs()
        plotter.agb_spatial_plot_standard(self.cframe,self.mframe,8.300500,48.50850)
    
    def plot_both_lum(self):
        plotter=graphs()
        plotter.k_luminosity_function(self.cframe,self.mframe)
    
    def plot_both_contour(self):
        plotter=graphs()
        plotter.surface_density_plot(self.cframe,self.mframe)
    
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
        

r=basic_graphs('m32')
r.plot_cc()
