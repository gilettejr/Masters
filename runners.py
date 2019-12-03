#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 20:56:08 2019

@author: cameronrobertson
"""
import numpy as np
from asymphr import asymphr
from graphing_class import graphs
from subset_utils import make_subsets
from iso_utils import import_isos
from crossmatching_utils import topcatcross,crossed_data


class run_cross:
    
    def __init__(self,galaxy):
        
        self.galaxy=galaxy
        
        #object for reading in and deleting crossmatches  initialised
        
        rmatch=crossed_data()
        
        #datafile chosen depending on galaxy. Crossmatched data
        #identified and read in by rmatch object
        
        if galaxy=='ngc147':
        
            n147topcross = topcatcross('lot_n147.unique',0)
            n147vizcross=topcatcross('lot_n147.unique',0)
            n147culled=topcatcross('lot_n147.unique',0)
            n147culled.loadascii()
            n147culled.ciscuts()
            rmatch.topmatch(n147topcross,'crossedn147.csv')
            rmatch.vizmatch(n147vizcross,'n147_vizcross.csv')
            
            
            
        elif galaxy== 'ngc185':
            n147cross = topcatcross('lot_n185.unique',0)
            rmatch.topmatch(n147cross,'crossedn185.csv')
            
        elif galaxy=='ngc205':
            n147cross=topcatcross('N205_two.asc',0)
            rmatch.topmatch(n147cross,'crossedn205.csv')
            
        elif galaxy=='m32':
            n147cross=topcatcross('M32.asc',0)
            rmatch.topmatch(n147cross,'crossedm32.csv')
        
        
        #crossmathed points deleted
        
        
        #tangent point topcatcross class attributes created
        
        n147topcross.create_tangent_coords(8.300500,48.50850)
        
        #dataframe set as topcatcross class attribute
        
        n147topcross.make_dataframe()
        n147vizcross.make_dataframe()
        n147culled.make_dataframe()
        
        #topcatcross object set as self class object attribute
        
        self.n147topcross=n147topcross
        self.n147vizcross=n147vizcross
        self.n147culled=n147culled
    def plot_topcross_cmd(self):
        plotter=graphs()
        plotter.plot_topmatch_cmd(self.n147topcross)
        
    def plot_vizcross_cmd(self):
        plotter=graphs()
        plotter.plot_vizmatch_cmd(self.n147vizcross)
        
    def plot_combined_cmd(self):
        self.n147topcross.delete_crossed_points()
        self.n147vizcross.delete_uncrossed_points()
        plotter=graphs()
        plotter.plot_viztop_cmd(self.n147topcross,self.n147vizcross)
        
    def plot_culled_cmd(self):
        self.n147topcross.delete_crossed_points()
        self.n147vizcross.delete_uncrossed_points()     
        self.n147culled.gaia_viz_cull_points(self.n147topcross,self.n147vizcross)
        self.n147culled.sbsextinction()
        
        plotter=graphs()
        plotter.plot_culled_viztop_cmd(self.n147culled)
        



class run_both:
    
    #initialising class creates two dataframes of data for m and c stars respectively
    #fully extinction corrected
    
    def __init__(self,galaxy):
        
        #two make_subsets instances created, one for each star type
        
        runm=make_subsets(galaxy)
        runc=make_subsets(galaxy)
        
        #m and c subsets created
        
        runm.mc('m')
        runc.mc('c')
        
        #subsets and galaxy set to class attributes
        
        self.mframe=runm.subset
        self.cframe=runc.subset
        self.galaxy=galaxy
        
        #plotter object set for plotting graphs
        
        self.plotter=graphs()
        

    #plots c and m star spatial distribution in tangent point coordinates
        
    def plot_both_spatial(self):

        #galaxy variable read in from class attribute
        
        galaxy=self.galaxy
        
        #tangent point coordinates defined depending on galaxy
        
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
            
            #error printed if incorrect galaxy used
            
            print('Not a valid object')
        
        #spatial distribution plotted
        self.plotter.agb_spatial_plot_standard(self.cframe,self.mframe,tra,tdec)
    
    #plots agb star luminosity functions
    
    def plot_both_lum(self):
        self.plotter.k_luminosity_function(self.cframe,self.mframe)
    
    def plot_both_contour(self):

        self.plotter.surface_density_plot(self.cframe,self.mframe)

    def plot_both_colour_hist(self):
        self.plotter.colour_hist(self.cframe,self.mframe)
        
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
        
        print(FEH)
        return(FEH)
        
    def find_CM_border(self):
        plotter=graphs()
        plotter.unbound_hist(self.cframe,self.mframe)
        
        
        
        

        
    
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
        

class run_rgb:
    
    def __init__(self,galaxy):
        run=make_subsets(galaxy)
        run.rgb()
        self.frame=run.subset
    
    def find_tip(self):
        plotter=graphs()
        plotter.rgb_tip_finder(self.frame)

