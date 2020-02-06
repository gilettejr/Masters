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
            n147cross=topcatcross('N205_new_trimmed.unique',0)
            rmatch.topmatch(n147cross,'crossedn205.csv')
            
        elif galaxy=='m32':
            n147cross=topcatcross('M32_new.asc',0)
            rmatch.topmatch(n147cross,'crossedm32.csv')
            
        elif galaxy=='andromeda':
            n147topcross = topcatcross('lot_m31.unique',0)
            n147vizcross=topcatcross('lot_m31.unique',0)
            n147culled=topcatcross('lot_m31.unique',0)
            n147culled.loadascii()
            n147culled.ciscuts()
            rmatch.topmatch(n147topcross,'crossedm31.csv')
            rmatch.vizmatch(n147vizcross,'m31_vizcross.csv')
        
        
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
        self.n147vizcross.delete_uncrossed_points_in_defined_area()
        plotter=graphs()
        plotter.plot_viztop_cmd(self.n147topcross,self.n147vizcross)
        
    def plot_culled_cmd(self):
        self.n147topcross.delete_crossed_points()
        self.n147vizcross.delete_uncrossed_points_in_defined_area('vizdat.csv')   
        self.n147culled.gaia_viz_cull_points(self.n147topcross,self.n147vizcross)
        self.n147culled.sbsextinction()
        
        plotter=graphs()
        plotter.plot_culled_viztop_cmd(self.n147culled)
        

    def kde(self):
        self.n147topcross.delete_crossed_points()
        self.n147vizcross.delete_uncrossed_points()     
        self.n147culled.gaia_viz_cull_points(self.n147topcross,self.n147vizcross)
        self.n147culled.sbsextinction()
        
        
        
class kde_separator:
    
    def __init__(self,galaxy):
        
        run=make_subsets(galaxy)
        
        run.crossed()
        
        self.frame=run.subset
        self.galaxy=galaxy
        
        self.plotter=graphs()
        
    def kde_graph_test(self):
        
        frame=self.frame
        
        h=frame.hmag
        k=frame.kmag
        j=frame.jmag
        
        
        
        upper=18.4
        lower=18
        
        while lower > 15:
            
            jk= j-k
            hk= h-k
            binjk=[]
            for i in range(len(j)): 
                
                if lower<j[i]<upper:
                    binjk.append(jk[i])
            
            binjk=np.array(binjk)
            print(len(binjk))
            self.plotter.cmd_kde(upper,lower,binjk)
            
            upper=upper-0.4
            lower=lower-0.4
            
    def kde_big(self):
        
        frame=self.frame
        
        h=frame.hmag
        k=frame.kmag
        j=frame.jmag
        
        
        
        upper=18
        lower=17
        
        left=1
        right=3000

        jk= j-k
        hk= h-k
        binjk=[]
        
        for i in range(len(j)):
            
            if lower<k[i]<upper and left<jk[i]<right:
                binjk.append(jk[i])

                
        binjk=np.array(binjk)
        print(len(binjk))

        self.plotter.cmd_kde(upper,lower,binjk)
                
        self.plotter.kj_cmd(frame)
        
    def hist(self):
        
        frame=self.frame
        
        h=frame.hmag
        k=frame.kmag
        j=frame.jmag
        
        
        
        upper=18
        lower=17
        


        jk= j-k
        hk= h-k
        binjk=[]
        binhk=[]
        
        for i in range(len(j)):
            
            if lower<k[i]<upper:
                binjk.append(jk[i])

                
        binjk=np.array(binjk)
        
        for i in range(len(h)):
            
            if lower<k[i]<upper:
                binhk.append(hk[i])
        
        print(len(binjk))

        self.plotter.cmd_hist(upper,lower,binjk,'K-J CMD')
        
        self.plotter.cmd_hist(upper,lower,binhk,'K-H CMD')
                

        
        
        

        
        
        
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
        
        self.mframe_top=runm.top_subset
        self.cframe_top=runc.top_subset
        
        self.mframe_viz=runm.viz_subset
        self.cframe_viz=runc.viz_subset
        
        
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
        
        def CM_to_FEH(CM):
        
            FEH=-1.39 -0.47*np.log10(CM)
        
            return(FEH)
        
        m_no_top=[]
        for i in self.mframe_top.kmag:
            if np.isnan(i)==False:
                m_no_top.append(0)
        c_no_top=[]   
        for j in self.cframe_top.kmag:
            if np.isnan(j)==False:
                c_no_top.append(0)
        
        C_top=len(c_no_top)
        M_top=len(m_no_top)
        
        ratio_top = C_top/M_top
        FEH_top=CM_to_FEH(ratio_top)
        
        m_no_viz=[]
        for i in self.mframe_viz.kmag:
            if np.isnan(i)==False:
                m_no_viz.append(0)
        c_no_viz=[]   
        for j in self.cframe_viz.kmag:
            if np.isnan(j)==False:
                c_no_viz.append(0)
        
        C_viz=len(c_no_viz)
        M_viz=len(m_no_viz)
        
        ratio_viz = C_viz/M_viz
        FEH_viz=CM_to_FEH(ratio_viz)
        
        print('Gaia crossmatch: C/M = ' + str(ratio_top) + ' and [Fe/H] = ' + str(FEH_top))
        print('Spitzer crossmatch: C/M = ' + str(ratio_viz) + ' and [Fe/H] = ' + str(FEH_viz))
        
        

        
    def c_over_m_grad(self,border,tra,tdec):
        
        def CM_to_FEH(CM):
        
            FEH=-1.39 -0.47*np.log10(CM)
        
            return(FEH)
        
        in_m_no_top=[]
        for i in range(len(self.mframe_top.kmag)):
            if np.isnan(self.mframe_top.ra[i])==False and np.sqrt((self.mframe_top.ra[i]-tra)**2+(self.mframe_top.dec[i]-tdec)**2) < border/3600 :
                in_m_no_top.append(0)
        in_c_no_top=[]   
        for i in range(len(self.mframe_top.kmag)):
            if np.isnan(self.cframe_top.ra[i])==False and np.sqrt((self.cframe_top.ra[i]-tra)**2+(self.cframe_top.dec[i]-tdec)**2) < border/3600 :
                in_c_no_top.append(0)
        
        in_C=len(in_c_no_top)
        in_M=len(in_m_no_top)
        
        out_m_no_top=[]
        for i in range(len(self.mframe_top.kmag)):
            if np.isnan(self.mframe_top.ra[i])==False and np.sqrt((self.mframe_top.ra[i]-tra)**2+(self.mframe_top.dec[i]-tdec)**2) > border/3600 :
                out_m_no_top.append(0)
        out_c_no_top=[]  
        
        for i in range(len(self.mframe_top.kmag)):
            if np.isnan(self.cframe_top.ra[i])==False and np.sqrt((self.cframe_top.ra[i]-tra)**2+(self.cframe_top.dec[i]-tdec)**2) > border/3600 :
                out_c_no_top.append(0)
        
        out_C=len(out_c_no_top)
        out_M=len(out_m_no_top)
        
        self.inCM=(in_C/in_M)
        self.outCM=(out_C/out_M)
        
        print(self.inCM)
        print(self.outCM)
        

        
        in_m_no_viz=[]
        for i in range(len(self.mframe_viz.kmag)):
            if np.isnan(self.mframe_viz.ra[i])==False and np.sqrt((self.mframe_viz.ra[i]-tra)**2+(self.mframe_viz.dec[i]-tdec)**2) < border/3600 :
                in_m_no_viz.append(0)
        in_c_no_viz=[]   
        for i in range(len(self.mframe_viz.kmag)):
            if np.isnan(self.cframe_viz.ra[i])==False and np.sqrt((self.cframe_viz.ra[i]-tra)**2+(self.cframe_viz.dec[i]-tdec)**2) < border/3600 :
                in_c_no_viz.append(0)
        
        in_C_viz=len(in_c_no_viz)
        in_M_viz=len(in_m_no_viz)
        
        out_m_no_viz=[]
        for i in range(len(self.mframe_viz.kmag)):
            if np.isnan(self.mframe_viz.ra[i])==False and np.sqrt((self.mframe_viz.ra[i]-tra)**2+(self.mframe_viz.dec[i]-tdec)**2) > border/3600 :
                out_m_no_viz.append(0)
        out_c_no_viz=[]  
        
        for i in range(len(self.mframe_viz.kmag)):
            if np.isnan(self.cframe_viz.ra[i])==False and np.sqrt((self.cframe_viz.ra[i]-tra)**2+(self.cframe_viz.dec[i]-tdec)**2) > border/3600 :
                out_c_no_viz.append(0)
        
        out_C_viz=len(out_c_no_viz)
        out_M_viz=len(out_m_no_viz)
        
        self.inCM_viz=(in_C_viz/in_M_viz)
        self.outCM_viz=(out_C_viz/out_M_viz)
        
        print(self.inCM_viz)
        print(self.outCM_viz)
        
        print('Using Gaia Crossmatching, inner galaxy (<70"), C/M=' + str(self.inCM) + ', [Fe/H]='+str(CM_to_FEH(self.inCM)))
        print('Using Gaia Crossmatching, outer galaxy (>70"), C/M=' + str(self.outCM) + ', [Fe/H]='+str(CM_to_FEH(self.outCM)))
        
        print('Using Spitzer Crossmatching, inner galaxy (<70"), C/M=' + str(self.inCM_viz) + ', [Fe/H]='+str(CM_to_FEH(self.inCM_viz)))
        print('Using Spitzer Crossmatching, outer galaxy (>70"), C/M=' + str(self.outCM_viz) + ', [Fe/H]='+str(CM_to_FEH(self.outCM_viz)))
        
        
        

        
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

