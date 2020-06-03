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

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


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
            n147vizcross.make_dataframe()
            n147culled.make_dataframe()
            self.n147vizcross=n147vizcross
            self.n147culled=n147culled
            
            
            
        elif galaxy== 'ngc185':
            n147topcross = topcatcross('lot_n185.unique',0)
            rmatch.topmatch(n147topcross,'crossedn185.csv')
            
        elif galaxy=='ngc205':
            n147topcross=topcatcross('N205_new_trimmed.unique',0)
            rmatch.topmatch(n147topcross,'crossedn205.csv')
            
        elif galaxy=='m32':
            n147topcross=topcatcross('M32_new.asc',0)
            rmatch.topmatch(n147topcross,'crossedm32.csv')
            
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
        n147topcross.sbsextinction()
        n147topcross.create_tangent_coords(8.300500,48.50850)
        
        #dataframe set as topcatcross class attribute
        
        n147topcross.make_dataframe()

        
        #topcatcross object set as self class object attribute
        
        self.n147topcross=n147topcross

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
        
    def plot_gaia_culled_cmd(self):
        self.n147topcross.delete_crossed_points()
        plotter=graphs()
        a=asymphr('N205_new_trimmed.unique',0)
        a.loadascii()
        a.ciscuts()
        a.sbsextinction()
        
        plotter.plot_viztop_cmd(a,self.n147topcross)
        
    def plot_culled_cmd(self):
        self.n147topcross.delete_crossed_points()
        self.n147vizcross.delete_uncrossed_points_in_defined_area('vizdat.csv')   
        self.n147culled.gaia_viz_cull_points(self.n147topcross,self.n147vizcross)
        self.n147culled.sbsextinction()
        
        plotter=graphs()
        plotter.plot_culled_viztop_cmd(self.n147culled)
    
    def plot_gaiacrossed_cmd(self):
        self.n147topcross.delete_crossed_points()
        colour=self.n147topcross.jmag-self.n147topcross.kmag
        #for i in range(len(self.n147topcross.kmag)):
            #if self.n147topcross.kmag[i] > 18 or self.n147topcross.kmag[i] < 17 or colour[i] > 1:
                #self.n147topcross.kmag[i] = np.nan
                #self.n147topcross.jmag[i] = np.nan
        #points=[]
        
        #for i in self.n147topcross.kmag:
            #if np.isnan(i)==False:
                #points.append(0)
        #print(str(len(points)) + ' datapoints')
        plotter=graphs()
        plotter.kj_cmd(self.n147topcross)
        
    def plot_gaiacrossed_cc(self):
        self.n147topcross.delete_crossed_points()

        
        plotter=graphs()
        plotter.colour_colour(self.n147topcross)
    
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
        
        
        
        upper=18.2
        lower=18
        
        while lower > 15.6:
            
            jk= j-k
            hk= h-k
            binjk=[]
            for i in range(len(j)): 
                
                if lower<j[i]<upper:
                    binjk.append(jk[i])
            
            binjk=np.array(binjk)
            print(len(binjk))
            self.plotter.cmd_kde(upper,lower,binjk)
            
            upper=upper-0.2
            lower=lower-0.2
            
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
            
        elif galaxy=='ngc205':
            
            self.tra=10.09189356
            self.tdec=41.68541564
            
        
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
        self.plotter.agb_spatial_plot_standard(self.cframe_top,self.mframe_top,tra,tdec)
    
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
    
    def c_over_m_grad_plot(self,slice_size,outer_limit):
        
        def CM_to_FEH(CM):
        
            FEH=-1.39 -0.47*np.log10(CM)
        
            return(FEH)
            
        tra=self.tra
        tdec=self.tdec
        
        grad=slice_size
        m_no_grad=[0]
        c_no_grad=[0]
        while grad < outer_limit:
            
            in_m_no_top=[]
            for i in range(len(self.mframe_top.kmag)):
                if np.isnan(self.mframe_top.ra[i])==False and np.sqrt((self.mframe_top.ra[i]-tra)**2+(self.mframe_top.dec[i]-tdec)**2) < grad/3600 :
                    in_m_no_top.append(0)
                    
            in_c_no_top=[]
            for i in range(len(self.mframe_top.kmag)):
                if np.isnan(self.cframe_top.ra[i])==False and np.sqrt((self.cframe_top.ra[i]-tra)**2+(self.cframe_top.dec[i]-tdec)**2) < grad/3600 :
                    in_c_no_top.append(0)
        
                
            m_no_grad.append(len(in_m_no_top))
            c_no_grad.append(len(in_c_no_top))
        

        
        
            grad=grad+slice_size
            
 
        

           
        m_slices=[]
        for i in range(1,len(m_no_grad)):
            m_slices.append(m_no_grad[i]-m_no_grad[i-1])
        c_slices=[]
        for j in range(1,len(c_no_grad)):
            c_slices.append(c_no_grad[j]-c_no_grad[j-1])
            
        print(m_slices)
        print(c_slices)
        
        m_slices=np.array(m_slices)    
        c_slices=np.array(c_slices)
        
        CM_slices=c_slices/m_slices
        
        FEH_slices=CM_to_FEH(CM_slices)
        
        xdata=[]
        for i in range(slice_size,(len(m_slices)+1)*slice_size,slice_size):
            xdata.append(i)
            
        xdata=np.array(xdata)
        
        print(xdata)
        print(FEH_slices)
        plt.rc('axes',labelsize=15)
        plt.errorbar(xdata,FEH_slices,xerr=0.2,yerr=0.2,capsize=10,ecolor='black',linestyle='none',marker='x',color='black',markersize=10)
        m,b=np.polyfit(xdata,FEH_slices,1)
        plt.plot(xdata,m*xdata+b)
        plt.xlabel('Radial Distance/arcsec')
        plt.ylabel('[Fe/H]')
        
        plt.show()
        
        
    #same as above method, but splits so I can look at how m31 gradient affects things    
    def c_over_m_split_plot(self,slice_size,outer_limit):
        
        def CM_to_FEH(CM):
        
            FEH=-1.39 -0.47*np.log10(CM)
        
            return(FEH)
            
        tra=self.tra
        tdec=self.tdec
        
        vertex1=(9.697,41.389)
        vertex2=(10.4909,41.389)
        vertex3=(10.4909,41.989)
        
        nandvertex1=(tra,tdec)
        nandvertex2=(9.697,tdec)
        nandvertex3=(9.697,41.989)
        nandvertex4=(tra,41.989)
        
        andvertex1=(tra,tdec)
        andvertex2=(10.4909,tdec)
        andvertex3=(10.4909,41.389)
        andvertex4=(tra,41.389)
        
        
        
        
        
        
        grad_seg=Polygon([vertex1,vertex2,vertex3])
        
        nand_seg=Polygon([nandvertex1,nandvertex2,nandvertex3,nandvertex4])
        and_seg=Polygon([andvertex1,andvertex2,andvertex3,andvertex4])   
        
        mframe_and=self.mframe_top.copy()
        cframe_and=self.cframe_top.copy()
        mframe_nand=self.mframe_top.copy()
        cframe_nand=self.cframe_top.copy()
        
        for i in range(len(mframe_and.ra)):
            
            if and_seg.contains(Point(mframe_and.ra[i],mframe_and.dec[i]))== False:
                
                mframe_and.loc[i]=np.nan
                
        for i in range(len(mframe_nand.ra)):
            
            if nand_seg.contains(Point(mframe_nand.ra[i],mframe_nand.dec[i]))== False:
                
                mframe_nand.loc[i]=np.nan
                
        for i in range(len(cframe_and.ra)):
            
            if and_seg.contains(Point(cframe_and.ra[i],cframe_and.dec[i]))==False:
                
                cframe_and.loc[i]=np.nan
                
        for i in range(len(cframe_nand.ra)):
            
            if nand_seg.contains(Point(cframe_nand.ra[i],cframe_nand.dec[i]))==False:
                
                cframe_nand.loc[i]=np.nan
            
                

        
        grad=slice_size
        m_no_grad=[0]
        c_no_grad=[0]
        while grad < outer_limit:
            
            in_m_no_top=[]
            for i in range(len(mframe_and.kmag)):
                if np.isnan(mframe_and.ra[i])==False and np.sqrt((mframe_and.ra[i]-tra)**2+(mframe_and.dec[i]-tdec)**2) < grad/3600 :
                    in_m_no_top.append(0)
                    
            in_c_no_top=[]
            for i in range(len(cframe_and.kmag)):
                if np.isnan(cframe_and.ra[i])==False and np.sqrt((cframe_and.ra[i]-tra)**2+(cframe_and.dec[i]-tdec)**2) < grad/3600 :
                    in_c_no_top.append(0)
        
                
            m_no_grad.append(len(in_m_no_top))
            c_no_grad.append(len(in_c_no_top))
        

        
        
            grad=grad+slice_size
            
 
        

           
        m_slices=[]
        for i in range(1,len(m_no_grad)):
            m_slices.append(m_no_grad[i]-m_no_grad[i-1])
        c_slices=[]
        for j in range(1,len(c_no_grad)):
            c_slices.append(c_no_grad[j]-c_no_grad[j-1])
            
        print(m_slices)
        print(c_slices)
        
        m_slices=np.array(m_slices)    
        c_slices=np.array(c_slices)
        
        CM_slices=c_slices/m_slices
        
        FEH_slices=CM_to_FEH(CM_slices)
        
        xdata=[]
        for i in range(slice_size,(len(m_slices)+1)*slice_size,slice_size):
            xdata.append(i)
            
        xdata=np.array(xdata)
        
        print(xdata)
        print(FEH_slices)
        
        xdata_and=xdata
        FEH_slices_and=FEH_slices
        
        grad=slice_size
        m_no_grad=[0]
        c_no_grad=[0]
        while grad < outer_limit:
            
            in_m_no_top=[]
            for i in range(len(mframe_nand.kmag)):
                if np.isnan(mframe_nand.ra[i])==False and np.sqrt((mframe_nand.ra[i]-tra)**2+(mframe_nand.dec[i]-tdec)**2) < grad/3600 :
                    in_m_no_top.append(0)
                    
            in_c_no_top=[]
            for i in range(len(cframe_nand.kmag)):
                if np.isnan(cframe_nand.ra[i])==False and np.sqrt((cframe_nand.ra[i]-tra)**2+(cframe_nand.dec[i]-tdec)**2) < grad/3600 :
                    in_c_no_top.append(0)
        
                
            m_no_grad.append(len(in_m_no_top))
            c_no_grad.append(len(in_c_no_top))
        
        

        
        
            grad=grad+slice_size
            
 
        

           
        m_slices=[]
        for i in range(1,len(m_no_grad)):
            m_slices.append(m_no_grad[i]-m_no_grad[i-1])
        c_slices=[]
        for j in range(1,len(c_no_grad)):
            c_slices.append(c_no_grad[j]-c_no_grad[j-1])
            
        print(m_slices)
        print(c_slices)
        
        m_slices=np.array(m_slices)    
        c_slices=np.array(c_slices)
        
        CM_slices=c_slices/m_slices
        
        FEH_slices=CM_to_FEH(CM_slices)
        
        xdata=[]
        for i in range(slice_size,(len(m_slices)+1)*slice_size,slice_size):
            xdata.append(i)
            
        xdata=np.array(xdata)
        
        print(xdata)
        print(FEH_slices)
        
        xdata_nand=-xdata
        FEH_slices_nand=FEH_slices
        
        xdata=np.concatenate((xdata_nand,xdata_and))
        FEH=np.concatenate((FEH_slices_nand,FEH_slices_and))
        
        xdatleft=[]
        FEHleft=[]
        xdatright=[]
        FEHright=[]
        for i in range(len(xdata)):
            if xdata[i]<0:
                xdatleft.append(xdata[i])
                FEHleft.append(FEH[i])
            elif xdata[i]>0:
                xdatright.append(xdata[i])
                FEHright.append(FEH[i])
        
        xdatleft=np.array(xdatleft)
        FEHdatleft=np.array(FEHleft)
        xdatright=np.array(xdatright)
        FEHdatright=np.array(FEHright)
        
        plt.rc('axes',labelsize=15)
        plt.errorbar(xdata,FEH,yerr=0.2,ecolor='black',marker='x',color='black',linestyle='none',markersize=10,capsize=10)
        m1,b1=np.polyfit(xdatleft,FEHdatleft,1)
        m2,b2=np.polyfit(xdatright,FEHdatright,1)
        plt.plot(xdatleft,m1*xdatleft+b1)
        plt.plot(xdatright,m2*xdatright+b2)
        
        
                
        plt.xlabel('Radial Distance/arcsec')
        plt.ylabel('[Fe/H]')
        
        plt.show()
                
        
                
            
            
            
    
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

