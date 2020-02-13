#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 20:56:08 2019

@author: cameronrobertson
"""
import numpy as np
import matplotlib.pyplot as plt
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
        
        if galaxy=='ngc147':
        
            self.mframe_viz=runm.viz_subset
            self.cframe_viz=runc.viz_subset
        
        
            self.tra=8.3005
            self.tdec=48.5087389
        
        elif galaxy=='ngc185':
            self.tra=9.7415417
            self.tdec=48.3373778
            
        
        #plotter object set for plotting graphs
        
        self.plotter=graphs()
        self.galaxy=galaxy
        

    #plots c and m star spatial distribution in tangent point coordinates
        
    def plot_both_spatial(self):

        #galaxy variable read in from class attribute
        
        galaxy=self.galaxy
        
        #tangent point coordinates defined depending on galaxy
        
        tra=self.tra
        tdec=self.tdec
        
        #spatial distribution plotted
        self.plotter.agb_spatial_plot_standard(self.cframe,self.mframe,tra,tdec)
    
    #plots agb star luminosity functions
    
    def plot_both_lum(self):
        self.plotter.k_luminosity_function(self.cframe,self.mframe)
    
    def plot_both_contour(self):

        self.plotter.surface_density_plot(self.cframe,self.mframe)

    def plot_both_colour_hist(self):
        self.plotter.colour_hist(self.cframe,self.mframe)
        
    def c_over_m_topviz(self):
        
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
        
    def c_over_m_top(self):
        
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
        
        

        
    def c_over_m_grad_topviz(self,border):
        
        def CM_to_FEH(CM):
        
            FEH=-1.39 -0.47*np.log10(CM)
        
            return(FEH)
            
        tra=self.tra
        tdec=self.tdec
        
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
        
        
        print('Using Gaia Crossmatching, inner galaxy (<70"), C/M=' + str(self.inCM) + ', [Fe/H]='+str(CM_to_FEH(self.inCM)))
        print('Using Gaia Crossmatching, outer galaxy (>70"), C/M=' + str(self.outCM) + ', [Fe/H]='+str(CM_to_FEH(self.outCM)))
        
        print('Using Gaia and Spitzer Crossmatching, inner galaxy (<70"), C/M=' + str(self.inCM_viz) + ', [Fe/H]='+str(CM_to_FEH(self.inCM_viz)))
        print('Using Gaia Spitzer Crossmatching, outer galaxy (>70"), C/M=' + str(self.outCM_viz) + ', [Fe/H]='+str(CM_to_FEH(self.outCM_viz)))
        
    def c_over_m_grad_top(self,border):
        
        def CM_to_FEH(CM):
        
            FEH=-1.39 -0.47*np.log10(CM)
        
            return(FEH)
            
        tra=self.tra
        tdec=self.tdec
        
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
        
        print('Using Gaia Crossmatching, inner galaxy (<70"), C/M=' + str(self.inCM) + ', [Fe/H]='+str(CM_to_FEH(self.inCM)))
        print('Using Gaia Crossmatching, outer galaxy (>70"), C/M=' + str(self.outCM) + ', [Fe/H]='+str(CM_to_FEH(self.outCM)))
        
    def mc_grad_colours(self,border):
        
        tra=self.tra
        tdec=self.tdec
        
        in_mframe_top_jmag=[]
        in_mframe_top_hmag=[]
        in_mframe_top_kmag=[]
        out_mframe_top_jmag=[]
        out_mframe_top_hmag=[]
        out_mframe_top_kmag=[]
        for i in range(len(self.mframe_top.kmag)):
            
            if np.isnan(self.mframe_top.ra[i])==True:
                
                continue
            
            elif np.sqrt((self.mframe_top.ra[i]-tra)**2+(self.mframe_top.dec[i]-tdec)**2) < border/3600 :
                
                in_mframe_top_jmag.append(self.mframe_top.jmag[i])
                in_mframe_top_hmag.append(self.mframe_top.hmag[i])
                in_mframe_top_kmag.append(self.mframe_top.kmag[i])
            
            else:
                out_mframe_top_jmag.append(self.mframe_top.jmag[i])
                out_mframe_top_hmag.append(self.mframe_top.hmag[i])
                out_mframe_top_kmag.append(self.mframe_top.kmag[i])
        
        
        in_m_top_j=np.array(in_mframe_top_jmag)
        in_m_top_h=np.array(in_mframe_top_hmag)
        in_m_top_k=np.array(in_mframe_top_kmag)
        in_m_colour=in_m_top_j-in_m_top_k
        
        out_m_top_j=np.array(out_mframe_top_jmag)
        out_m_top_h=np.array(out_mframe_top_hmag)
        out_m_top_k=np.array(out_mframe_top_kmag)
        out_m_colour=out_m_top_j-out_m_top_k
        
        in_cframe_top_jmag=[]
        in_cframe_top_hmag=[]
        in_cframe_top_kmag=[]
        out_cframe_top_jmag=[]
        out_cframe_top_hmag=[]
        out_cframe_top_kmag=[]
        for i in range(len(self.cframe_top.kmag)):
            
            if np.isnan(self.cframe_top.ra[i])==True:
                
                continue
            
            elif np.sqrt((self.cframe_top.ra[i]-tra)**2+(self.cframe_top.dec[i]-tdec)**2) < border/3600 :
                
                in_cframe_top_jmag.append(self.cframe_top.jmag[i])
                in_cframe_top_hmag.append(self.cframe_top.hmag[i])
                in_cframe_top_kmag.append(self.cframe_top.kmag[i])
            
            else:
                out_cframe_top_jmag.append(self.cframe_top.jmag[i])
                out_cframe_top_hmag.append(self.cframe_top.hmag[i])
                out_cframe_top_kmag.append(self.cframe_top.kmag[i])
        
        
        in_c_top_j=np.array(in_cframe_top_jmag)
        in_c_top_h=np.array(in_cframe_top_hmag)
        in_c_top_k=np.array(in_cframe_top_kmag)
        in_c_colour=in_c_top_j-in_c_top_k
        
        out_c_top_j=np.array(out_cframe_top_jmag)
        out_c_top_h=np.array(out_cframe_top_hmag)
        out_c_top_k=np.array(out_cframe_top_kmag)
        out_c_colour=out_c_top_j-out_c_top_k
        
        fig,axs=plt.subplots(3,1)
        
        
        axs[0].scatter(in_m_colour,in_m_top_k,label='Inner Region M stars',marker='.',alpha=0.8)
        axs[0].scatter(in_c_colour,in_c_top_k,label='Inner Region C stars',marker='.',alpha=0.8)
        axs[1].scatter(out_m_colour,out_m_top_k,label='Outer Region M stars',marker='.',alpha=0.8)
        axs[1].scatter(out_c_colour,out_c_top_k,label='Outer Region C stars',marker='.',alpha=0.8)
        axs[2].scatter(self.mframe_top.ra,self.mframe_top.dec,label='M stars',marker='.')
        axs[2].scatter(self.cframe_top.ra,self.cframe_top.dec,label='C stars',marker='.')
        circle1=plt.Circle((tra,tdec),border/3600,alpha=0.5,color='red')
        axs[2].add_artist(circle1)
        
        axs[0].set_ylabel('$K_0$')
        axs[0].set_xlabel('$J_0$-$K_0$')
        
        if self.galaxy=='ngc147':
            xmin=0.9
            xmax=4.27
            ymin=15.7
            ymax=18
            
        elif self.galaxy=='ngc185':
            xmin=0.95
            xmax=3.54
            ymin=15.1
            ymax=17.8
            
        
        axs[0].set_xlim(xmin,xmax)
        axs[0].set_ylim(ymin,ymax)
        axs[1].set_xlim(xmin,xmax)
        axs[1].set_ylim(ymin,ymax)
            
        axs[0].invert_yaxis()
        
        axs[1].set_ylabel('$K_0$')
        axs[1].set_xlabel('$J_0$-$K_0$')
        
        axs[1].invert_yaxis()
        
        axs[2].set_ylabel('Dec')
        axs[2].set_xlabel('RA')
        
        axs[0].legend()
        axs[1].legend()
        axs[2].legend()
        plt.show()
            
            
        

        
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

