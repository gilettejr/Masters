#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 21:07:01 2019

@author: cameronrobertson
"""
import pandas as pd
import numpy as np
from astropy.table import QTable
from astropy import units as u
from asymphr import asymphr

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
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
        
        full_indices=np.arange(len(self.ra))
        
        unindex_nos=np.setdiff1d(full_indices,self.top_index_nos)
        
        self.unindex_nos=unindex_nos
        
        
        
        for i in self.top_index_nos:
        
            #ignores NaN values
            
            if np.isnan(i)==True:
                
                continue
            
            
            else:
                
                #index numbers are floats, need to be converted to integers for identifying array elements
                
                i=int(i)
                
                #entire row of data wiped if crossmatched with gaia data
                #self.index[i]=np.nan
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
                #self.eta[i] = np.nan
                #self.xi[i]=np.nan
            
            
        #function for deleting all uncrossmatched points
    
    def delete_uncrossed_points(self):
        
        #locates each crossed point using index identifiers from crossmatched gaia file
        
        for i in self.viz_index_nos:
            i=int(i)
        print(self.viz_index_nos)
        full_indices=np.arange(len(self.ra))
        print(full_indices)
        unindex_nos=np.setdiff1d(full_indices,self.viz_index_nos)
        self.unindex_nos=unindex_nos
        print(unindex_nos)
        
        for i in unindex_nos:
        
            #ignores NaN values
            
            if np.isnan(i)==True:
                
                continue
            
            
            else:
                
                #index numbers are floats, need to be converted to integers for identifying array elements
                
                
                #entire row of data wiped if crossmatched with gaia data
                #self.index[i]=np.nan
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
                
    def delete_uncrossed_points_in_defined_area(self,vizier_fov_file):
        
        vizcat=pd.read_csv(vizier_fov_file)
        
        ra=vizcat._RAJ2000
        dec=vizcat._DEJ2000
        
        vertex1=(ra[np.argmin(ra)],dec[np.argmin(ra)])
        vertex2=(ra[np.argmax(dec)],dec[np.argmax(dec)])
        vertex3=(ra[np.argmax(ra)],dec[np.argmax(ra)])
        vertex4=(ra[np.argmin(dec)],dec[np.argmin(dec)])
        

        
        fov=Polygon([vertex1,vertex2,vertex3,vertex4])
        
        #locates each crossed point using index identifiers from crossmatched gaia file
        
        for i in self.viz_index_nos:
            i=int(i)
        #print(self.viz_index_nos)
        full_indices=np.arange(len(self.ra))
        #print(full_indices)
        unindex_nos=np.setdiff1d(full_indices,self.viz_index_nos)
        self.unindex_nos=unindex_nos
        #print(unindex_nos)
        

        
        for i in unindex_nos:
        
            #ignores NaN values
            
            if np.isnan(i)==True:
                
                continue
            
            #formula to stop culling outside of spitzer fov
            elif fov.contains(Point(self.ra[i],self.dec[i])) == False:
                continue
            
            else:
                
                #index numbers are floats, need to be converted to integers for identifying array elements
                
                
                #entire row of data wiped if crossmatched with gaia data
                #self.index[i]=np.nan
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
            

    def gaia_viz_cull_points(self,topcatcross,vizcross):
        
        cull_index_nos=np.intersect1d(topcatcross.unindex_nos,vizcross.index_nos)
        
        full_indices=np.arange(len(self.ra))
        cull_unindex_nos=np.setdiff1d(full_indices,cull_index_nos)
        
        for i in cull_unindex_nos:
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
        
                        #self.eta[i] = np.nan
                        #self.xi[i]=np.nan
                
    #method for converting rows of data into pandas dataframe for ease
                
    def make_dataframe(self):
        
        #dataframe created
        
        wfdat=pd.DataFrame({'ra':self.ra,'dec':self.dec,'xj':self.xj,'yj':self.yj,'jmag':self.jmag,'jerr':self.jerr,'jcls':self.jcis,'xh':self.xh,'yh':self.yh,'hmag':self.hmag,'herr':self.herr,'hcls':self.hcis,'xk':self.xk,'yk':self.yk,'kmag':self.kmag,'kerr':self.kerr,'kcls':self.kcis,'xi':self.xi,'eta':self.eta})
        
        #attribute dataframe set
        
        self.wfdat=wfdat
    
    
    #returns dataframe only containing C star candidates based on hard colour cuts
    
    def select_C_stars(self,j_kmin,h_kmin,j_hmin,minkmag):
        
        #copy made to protect original dataframe
        
        d=self.wfdat.copy()
        
        #hard colour cuts recursively made
        
        for i in range(len(d.ra)):
            
            if np.isnan(d.ra[i]) == True:
                
                continue
            
            
            elif d.kmag[i] > minkmag or d.jmag[i]-d.kmag[i] < j_kmin or d.hmag[i]-d.kmag[i] < h_kmin  or d.jmag[i]-d.hmag[i] < j_hmin:
                
                #entire row wiped if conditions not met
                
                d.loc[i]=np.nan
                
        #altered dataframe copy returned        
                
        return d
    
    #returns dataframe only containing M star candidates based on hard colour cuts
    
    def select_M_stars(self,j_kmin,j_kmax,h_kmax,minkmag):
        
        #copy made to protect original dataframe
        
        d=self.wfdat.copy()
        
        #hard colour cuts recursively made
        
        for i in range(len(d.ra)):
            
            if np.isnan(d.ra[i]) == True:
                
                continue
            
            
            elif d.kmag[i] > minkmag or j_kmin > d.jmag[i]-d.kmag[i] or d.jmag[i]-d.kmag[i] > j_kmax or  d.hmag[i]-d.kmag[i] > h_kmax:
                
                #entire row wiped if conditions not met

                d.loc[i]=np.nan

        #altered dataframe copy returned
        
        return d
            
                        
    def select_rgb_region(self,xmin,xmax,ymin,ymax):
        
        d=self.wfdat.copy()
        colour = d.jmag-d.kmag
        kmag=d.kmag
        
        for i in range(len(kmag)):
            
            
            if xmin < colour[i] < xmax and ymin<kmag[i]<ymax:

                continue

            else:
                d.loc[i]=np.nan
        return d


#class for dealing with identifying data crossmatched by topcat
        
class crossed_data:
    
    #method sets crossmatched data from incrossfile to gaiacross class attributes
    #gaiacross must be topcatcross object
    
    def topmatch(self,gaiacross,incrossfile):
        
        #method called to read in crossmatched file as dataframe
        
        crossed_table = gaiacross.read_crossed_csv(incrossfile)
        
        #noise cut based on Gaia DR2 values and recomendations from the Gaia project website
        
        for i in range(len(crossed_table.parallax_over_error)):
            if crossed_table.astrometric_excess_noise_sig[i] > 2:
                crossed_table.index_no[i]=np.nan
                
        #parallax cut to remove any badly  measured parallaxes
                
        for i in range(len(crossed_table.parallax_over_error)):
            if crossed_table.parallax_over_error[i] < 1:
                crossed_table.index_no[i]=np.nan
                
        #proper motion cut to remove any proper motions pointing to distant objects (i.e. pm~0)
                
        #for i in range(len(crossed_table.pmra)):
            #if (np.abs(crossed_table.pmra[i]/crossed_table.pmra_error[i])) < 1 or (np.abs(crossed_table.pmdec[i]/crossed_table.pmdec_error[i])) < 1:
                #crossed_table.index_no[i]=np.nan

        
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
        
        gaiacross.topcrossjmag=crossjmag
        gaiacross.topcrosshmag=crosshmag
        gaiacross.topcrosskmag=crosskmag
        
        #index numbers also set as attribute for future troubleshooting
        
        gaiacross.top_index_nos=crossed_table.index_no

    def vizmatch(self,gaiacross,incrossfile):
        
        #method called to read in crossmatched file as dataframe
        
        crossed_table = gaiacross.read_crossed_csv(incrossfile)
        
        #WFCAM data is loaded and the usual cuts made
        
        gaiacross.loadascii()
        gaiacross.ciscuts()
        
        #lists for crossmatched photometric data created
        
        crossjmag = []
        crosshmag = []
        crosskmag = []

        #index_no column used to identify vizier crossmatched data
    
        #set index number to integer for indexing
    
        for j in crossed_table.index_no:
            j=(int(j))
            
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
        
        gaiacross.viz_index_nos=crossed_table.index_no
            
            
            
        #create magnitude arrays to visualise crossmatched points on cmd
        
        





    
    



        
#class for setting up data before topcat crossmatching

class topcatstuff:
    
    #uses topcatcross class methods to convert ra and dec in named
    #file to decimal degrees, and converts the file from .csv to .fits
    
    def fit_for_cross(self,file):
        n147cross=topcatcross(file,0)
        n147cross.loadascii()
        n147cross.ciscuts()
        n147cross.load_ascii_to_cross()
        

t=topcatstuff()

t.fit_for_cross('N205_new_trimmed.unique')