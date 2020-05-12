#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 20:46:42 2019

@author: cameronrobertson
"""
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from asymphr import asymphr
from astropy import stats
from astropy.stats.histogram import knuth_bin_width
from scipy.stats import gaussian_kde
from HESSCMD import plotHess
from crossmatching_utils import topcatcross
#class containing methods for plotting various graphs from analysis output

class graphs:
    
    def __init__(self):
        
        def quadrat(a,b):
            
            c=np.sqrt(a**2+b**2)
            
            return c
        self.mag=quadrat
            
            
        
    #method produces K-J CMD, taking in asymphr or asympher inherited objects as input    

    def kj_cmd(self,target):
    
    #figure created
        
        plt.figure()
                
        
        #plot formatted and labelled

        #plt.style.use('presentation')
        plt.rc('axes',labelsize = 15)
        plt.scatter(target.jmag-target.kmag,target.kmag,s=3,marker='o',alpha = 1,color='black')
        
        jerrval=[]
        herrval=[]
        kerrval=[]
        
        yerrpos=np.array([20,19,18,17,16,15,14,13])
        xerrpos=np.zeros(len(yerrpos))
        xerrpos=xerrpos-1
        
        for j in yerrpos:
            
            for i in range(len(target.jmag)):
                
                if np.isnan(target.jmag[i])==False:
                    
                    if (j-0.05) <target.kmag[i] < (j+0.05):
                        
                        jerrval.append(target.jerr[i])
                        herrval.append(target.herr[i])
                        kerrval.append(target.kerr[i])
                        
                        break
        
        jerrval=np.array(jerrval)
        herrval=np.array(herrval)
        kerrval=np.array(kerrval)
        
        print(kerrval)
        
        xerrval=self.mag(jerrval,kerrval)
        yerrval=kerrval
                
        #plt.errorbar(xerrpos,yerrpos,xerr=xerrval,yerr=yerrval,ecolor='red',linestyle='none')
        #plt.ylim(11,21)
        #plt.xlim(-1.5,4.5)
        plt.gca().invert_yaxis()
        plt.ylabel('$K_0$')
        plt.xlabel('$J_0$-$K_0$')
        
        #plt.title(target.filename + ' CMD')

        plt.show()
        
        #ngc147 cut
        
        #plt.axvline(x=1.34,linestyle=':',color='black')
        #plt.axhline(y=18,linestyle=':',color='black')
        
    def kj_cmd_select(self,target,xmin,xmax,ymin,ymax):
    
    #figure created
        
        plt.figure()
                
        
        #plot formatted and labelled

        #plt.style.use('presentation')
        plt.rc('axes',labelsize = 20)
        plt.scatter(target.jmag-target.kmag,target.kmag,s=3,marker='o',alpha = 1)
        plt.gca().set_autoscale_on(False)
        plt.plot([xmin,xmax],[ymin,ymin],color='black')
        plt.plot([xmax,xmax],[ymin,ymax],color='black')
        plt.plot([xmax,xmin],[ymax,ymax],color='black')
        plt.plot([xmin,xmin],[ymax,ymin],color='black')
        plt.gca().invert_yaxis()
        plt.ylabel('K')
        plt.xlabel('J-K')
        

        #plt.title(target.filename + ' CMD')
        #plt.style.use('ggplot')
        plt.show()
        
        #ngc147 cut
        
        #plt.axvline(x=1.34,linestyle=':',color='black')
        #plt.axhline(y=18,linestyle=':',color='black')
    
    def hk_cmd(self,target):
    
    #figure created
        
        plt.figure()
                
        
        #plot formatted and labelled

        #plt.style.use('presentation')
        plt.rc('axes',labelsize = 15)
        plt.scatter(target.hmag-target.kmag,target.kmag,s=3,marker='o',alpha = 1,color='black')
        
        jerrval=[]
        herrval=[]
        kerrval=[]
        
        yerrpos=np.array([20,19,18,17,16,15])
        xerrpos=np.zeros(len(yerrpos))
        xerrpos=xerrpos-1
        
        for j in yerrpos:
            
            for i in range(len(target.jmag)):
                
                if np.isnan(target.jmag[i])==False:
                    
                    if (j-0.05) <target.kmag[i] < (j+0.05):
                        
                        jerrval.append(target.jerr[i])
                        herrval.append(target.herr[i])
                        kerrval.append(target.kerr[i])
                        
                        break
        
        jerrval=np.array(jerrval)
        herrval=np.array(herrval)
        kerrval=np.array(kerrval)
        
        print(kerrval)
        
        xerrval=self.mag(herrval,kerrval)
        yerrval=kerrval
                
        plt.errorbar(xerrpos,yerrpos,xerr=xerrval,yerr=yerrval,ecolor='red',linestyle='none')
        plt.ylim(11,21)
        plt.xlim(-1.5,4.5)
        plt.gca().invert_yaxis()
        plt.ylabel('$K_0$')
        plt.xlabel('$H_0$-$K_0$')
        
        #plt.title(target.filename + ' CMD')

        plt.show()
        
        #ngc147 cut
        
        #plt.axvline(x=1.34,linestyle=':',color='black')
        #plt.axhline(y=18,linestyle=':',color='black')
    #method plots spatial distribution from target class or dataframe

    def spatial_plot(self,target):
        
        #cuts performed, extinction correction not necessary
        
        self.tra=10.09189356
        self.tdec=41.68541564
        
        #plot formatted and labelled
        plt.figure()
        
        plt.rc('axes',labelsize = 20)
        plt.scatter(target.ra,target.dec,s=3,marker = 'o',color='blue',alpha = 1)
        #plt.plot([self.tra,10.68470833],[self.tdec,41.26875])
        plt.xlabel('RA(J200)/deg')
        plt.ylabel('Dec(J200)/deg')
        #plt.title(target.filename + ' Spatial Plot')
        plt.style.use('ggplot')
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
        
        plt.rc('axes',labelsize = 15)
        plt.plot(xi,eta,linestyle='none',marker = 'o',markersize=1,color='black')
        plt.xlim(-0.3,0.3)
        plt.ylim(-0.3,0.3)
        plt.gca().invert_xaxis()
        plt.gca().set_ylabel(r'$\eta$/degrees')
        plt.gca().set_xlabel(r'$\xi$/degrees')
        
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
        #plt.style.use('seaborn-white')
        plt.rc('axes',labelsize = 15)
        plt.scatter(cxi,ceta,s=3,marker = 'o',alpha=1,color='black',label='C-Star Candidates')
        plt.scatter(mxi,meta,s=3,marker = 'o',alpha=1,color='black',label='M-Star Candidates')
        plt.xlim(-0.3,0.3)
        plt.ylim(-0.3,0.3)
        plt.gca().set_ylabel(r'$\eta$\degrees')
        plt.gca().set_xlabel(r'$\xi$\degrees')
        plt.gca().invert_xaxis()
        #plt.legend()
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
        
       
        
        #plot formatted and labelled
        
        plt.rc('axes',labelsize = 15)
        #graph is in the format: (j-k)/(h-k)
        plt.scatter(target.hmag - target.kmag,target.jmag-target.hmag,s=3,marker='o',color='black')
        plt.ylabel('$J_0$-$H_0$')
        plt.xlabel('$H_0$-$K_0$')
        #plt.title(target.filename + 'Colour-colour diagram')
        
        #plt.style.use('ggplot')
        
        #line for ngc147 hard colour C cut, Y-J sohn et al.
        
        #plt.axvline(x=0.44,linestyle=':',color='black')
        #plt.axhline(y=1.34,linestyle=':',color='black')
    
    
    
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
                    
        plt.rc('axes',labelsize=15)
        
        plt.scatter(cross.jmag-cross.kmag,cross.kmag,s=3,marker='o',alpha=1,color='black',label='Non-contaminants')
        plt.scatter(cross.topcrossjmag-cross.topcrosskmag,cross.topcrosskmag,s=3,marker='o',color='red',alpha=1,label='Foreground contaminants')
        
        
        
        plt.gca().invert_yaxis()
        plt.ylabel('$K_0$')
        plt.xlabel('$J_0$-$K_0$')
        #plt.style.use('ggplot')
        plt.legend()
        
        plt.show()
        
    #plots crossmatched data from vizier catalogue with uncrossmatched data in cmd
    #separates two subsets and labels before graphing
    
    def plot_vizmatch_cmd(self,cross):
                    
        plt.rc('axes',labelsize=15)
        
        plt.scatter(cross.jmag-cross.kmag,cross.kmag,s=3,marker='o',alpha=1,color='red',label='Foreground contaminants')
        plt.scatter(cross.crossjmag-cross.crosskmag,cross.crosskmag,s=3,marker='o',alpha=1,color='black',label='Non-contaminants')
        
        
        
        plt.gca().invert_yaxis()
        plt.ylabel('$K_0$')
        plt.xlabel('$J_0$-$K_0$')
        #plt.style.use('ggplot')
        plt.legend()
        
        plt.show()
        
    #graphing method to plot cmd from WFCAM data with Gaia crossmatched data subtracted, and Vizier uncrossmatched data subtra
    def plot_viztop_cmd(self,target,topcross):
        plt.rc('axes',labelsize=15)
        
        


        plt.scatter(target.jmag-target.kmag,target.kmag,s=3,marker='o',color='red',alpha=1,label='Foreground contaminants')
        plt.scatter(topcross.jmag-topcross.kmag,topcross.kmag,s=3,marker='o',alpha=1,color='black',label='Non-contaminants')

                
        plt.gca().invert_yaxis()
        plt.ylabel('$K_0$')
        plt.xlabel('$J_0$-$K_0$')
        plt.legend()
        
        plt.show()
        
    def plot_culled_viztop_cmd(self,cross):
        plt.scatter(cross.jmag-cross.kmag,cross.kmag,s=3,marker='o',alpha=1,label='UKIRT data culled from crossmatching with Gaia DR2 and Vizier IR Catalogue')
        
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
        sns.distplot(mframe.jmag.dropna()-mframe.kmag.dropna(),kde=False,label='M-Stars')
        sns.distplot(cframe.jmag.dropna()-cframe.kmag.dropna(),kde=False,bins=120,label='C-Stars')

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
        
    def cmd_kde(self,upper,lower,colour):
        plt.figure()
        
        print(colour)
        binwidth = stats.knuth_bin_width(colour)
        kde_data=colour
        #print(binwidth)
        bins=np.arange(min(kde_data ), max(kde_data ) + binwidth, binwidth)
        x_eval = np.linspace(kde_data.min() - 1.0, kde_data.max() + 1.0,500)
        kde=gaussian_kde(kde_data,bw_method=binwidth)
        
        plt.plot(x_eval,kde(x_eval),'k',lw=2,label='KDE')
        plt.hist(kde_data,bins=bins,density=True,label='Binned Data')
        plt.legend()
        plt.xlabel('$(J-K)_0$')
        plt.ylabel('Normalised Density')
        plt.show()
        
    def cmd_hist(self,upper,lower,colour, cmd_colour):
        
        
        plt.figure()
        binwidth = stats.knuth_bin_width(colour)
        kde_data=colour
        bins=np.arange(min(kde_data ), max(kde_data ) + binwidth, binwidth)

        kde_data=colour
        
        plt.hist(kde_data,bins=bins,label='Binned Data')
        plt.title(cmd_colour)

        
        
        
    
    def rgb_tip_finder(self,frame):
        data,bins=np.histogram(frame.kmag.dropna(),bins=30)
        xfunc=bins
        Xfunc=[]
        Yfunc=data
        
        for i in range(len(xfunc)-1):
            Xfunc.append(xfunc[i]+((xfunc[i+1])-xfunc[i])/2)
        
        print(len(Xfunc))
        print(len(Yfunc))
        plt.hist(frame.kmag,bins=15)
        fit=np.polyfit(Xfunc,Yfunc,3)
        print(fit)

    
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
        
    def unbound_hist(self,cframe,mframe):
        jmag=np.concatenate((cframe.jmag.dropna(),mframe.jmag.dropna()))
        kmag=np.concatenate((cframe.kmag.dropna(),mframe.kmag.dropna()))
        hmag=np.concatenate((cframe.hmag.dropna(),mframe.hmag.dropna()))
        
        plt.rc('axes',labelsize=20)
        
        plt.hist2d(jmag-kmag,hmag-kmag,bins=100)
        plt.ylabel('(H-K)$_0$')
        plt.xlabel('(J-K)$_0$')

#function cuts gaia data and sets class attributes for crossmatched data
        
#class for running graphing functions
#mainly involves calling the relevant graphs() method after class initialisation
class basic_graphs:
    
    #galaxy defined upon initialisation
    
    def __init__(self,galaxy):
        
        #galaxy set as class attribute
        
        self.galaxy=galaxy
        
        #conditional statements used to select correct datafile for the chosen galaxy
        
        if galaxy=='ngc147':
        
            n147 = topcatcross('lot_n147.unique',0)

            
        elif galaxy== 'ngc185':
            n147 = topcatcross('lot_n185.unique',0)

            
        elif galaxy=='ngc205':
            n147 = topcatcross('N205_new_trimmed.unique',0)

            
        elif galaxy=='m32':
            n147=topcatcross('M32_new.asc',0)
            
        elif galaxy=='andromeda':
            n147=topcatcross('lot_m31.unique',0)
            
        elif galaxy=='ngc205_old':
            n147=topcatcross('N205_two.asc',0)
            
        #error returned if galaxy not recognised 
        
        else:
            print('That is not a galaxy, enjoy the errors')
        
        #data loaded and cls cuts and extinction corrections carried out 
        
        n147.loadascii()
        n147.ciscuts()
        n147.sbsextinction()
        n147.make_dataframe()
        
        #data set as class attribute for running class methods
        
        self.n147=n147
        self.plotter=graphs()
    #uses the graphs().kj_cmd method to plot a CMD using class attribute data
        
    def plot_kj_cmd(self):
        self.plotter.kj_cmd(self.n147)
        
    def plot_hk_cmd(self):
        self.plotter.hk_cmd(self.n147)
    
    #uses the graphs().colour_colour method to plot a cc diagram using class attribute data
    
    def plot_cc(self):

        self.plotter.colour_colour(self.n147)
        
    #uses the graphs().spatial_plot method to plot a spatial distribution using class attribute data
    
    def plot_spatial(self):

        self.plotter.spatial_plot(self.n147)
        
    #uses the graphs().xy_spatial_plot method to plot an xy pixel distribution using class attribute data    
        
    def plot_xy_spatial(self):
        self.plotter.xy_spatial_plot(self.n147)
        
    def plot_hess(self):
        n147=self.n147.wfdat
        
        if self.galaxy=='ngc205' or self.galaxy=='ngc185':
        
            plotHess(n147.kmag.dropna(),n147.jmag.dropna()-n147.kmag.dropna(),colormap=plt.cm.plasma)
            
        elif self.galaxy=='andromeda':
            
            plotHess(n147.kmag.dropna(),n147.jmag.dropna()-n147.kmag.dropna(),colormap=plt.cm.plasma,levels=150)
            
            
    def plot_tan(self):
        
        if self.galaxy=='ngc147':
            tanra=8.3005
            tandec=48.5087389
        elif self.galaxy=='ngc185':
            tanra=9.7415417
            tandec=48.3373778
        elif self.galaxy=='ngc205'or self.galaxy=='ngc205_old':
            tanra=10.09189356
            tandec=41.68541564
        
        elif self.galaxy=='m32':
            tandec=40.8651694
            tanra=10.6742708
        
        self.plotter.spatial_plot_standard(self.n147,tanra,tandec)
            
    def plot_cchess(self):
        n147=self.n147.wfdat
        
        plotHess(n147.jmag.dropna()-n147.hmag.dropna(),n147.hmag.dropna()-n147.kmag.dropna(),colormap=plt.cm.plasma,levels=150)

        