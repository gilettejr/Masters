#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 21:09:10 2019

@author: cameronrobertson
"""

#class containing methods for separating out classes of star depending on
#phase or other user defined cuts. Automatically removes stars crossmatched from
#gaia DR2

from crossmatching_utils import crossed_data,topcatcross


class make_subsets:
    
    #class initialised with defined galaxy. Deletes gaia crossmatched data and
    #converts data into dataframe format, as well as adding tangent co-ordinates
    #still a work in progress, only works for ngc147 right now!
    
    def __init__(self,galaxy):
        
        self.galaxy=galaxy
        
        #object for reading in and deleting crossmatches  initialised
        
        rmatch=crossed_data()
        
        #datafile chosen depending on galaxy. Crossmatched data
        #identified and read in by rmatch object
        
        if galaxy=='ngc147':
        
            n147cross = topcatcross('lot_n147.unique',0)
            rmatch.topmatch(n147cross,'crossedn147.csv')
            n147vizcross=topcatcross('lot_n147.unique',0)
            rmatch.topmatch(n147vizcross,'crossedn147.csv')
            rmatch.vizmatch(n147vizcross,'n147_vizcross.csv')
            
            n147cross.delete_crossed_points()
            n147vizcross.delete_crossed_points()
            n147vizcross.delete_uncrossed_points_in_defined_area('vizdat.csv')
            n147cross.sbsextinction()
            n147vizcross.sbsextinction()
            
            n147cross.make_dataframe()
            n147vizcross.make_dataframe()
            
            self.n147cross=n147cross
            self.n147vizcross=n147vizcross


            
        elif galaxy== 'ngc185':
            n147cross = topcatcross('lot_n185.unique',0)
            rmatch.topmatch(n147cross,'crossedn185.csv')

            
            n147cross.delete_crossed_points()
            n147cross.sbsextinction()
            
            n147cross.make_dataframe()
            self.n147cross=n147cross

            
        elif galaxy=='ngc205':
            n147cross=topcatcross('N205_new_trimmed.unique',0)
            rmatch.topmatch(n147cross,'crossedn205.csv')
            
            n147cross.delete_crossed_points()
            n147cross.sbsextinction()
            
            n147cross.make_dataframe()
            self.n147cross=n147cross
            
        elif galaxy=='m32':
            n147cross=topcatcross('M32_new.asc',0)
            rmatch.topmatch(n147cross,'crossedm32.csv')
            
        elif galaxy=='andromeda':
            n147cross=topcatcross('lot_m31.unique',0)
            rmatch.topmatch(n147cross,'crossedm31.csv')
        
        
        #crossmathed points deleted
        


        
        #tangent point topcatcross class attributes created
        
        #n147cross.create_tangent_coords(8.300500,48.50850)
        
        #dataframe set as topcatcross class attribute
        

        
        #topcatcross object set as self class object attribute
        
        
    
    #function for selecting specifically m and c stars based on hard cuts
    
    def mc(self,agb):
        
        #attribute object set to variable for ease
        
        n147cross=self.n147cross

        
        if self.galaxy=='ngc147':
            
            n147vizcross=self.n147vizcross
        
        #appropriate topcatcross method utilised to take subset data from
        #wfcat dataframe attribute and select only c/m star candidates
        
            if agb=='m':
            
                #balanced selection
            
                selection_top=n147cross.select_M_stars(1.0,1.33,0.44,0.83,18)
                
                selection_viz=n147vizcross.select_M_stars(1.0,1.33,0.44,0.83,18)
                
                #low c selection
                
                #selection_top=n147cross.select_M_stars(0.95,1.38,0.55,18)
                
                #selection_viz=n147vizcross.select_M_stars(0.95,1.38,0.55,18)
                
            elif agb=='c':
                
                #balanced
                
                selection_top=n147cross.select_C_stars(1.33,0.44,0.83,18)
                selection_viz=n147vizcross.select_C_stars(1.33,0.44,0.83,18)
                
                #low c selection
                
                #selection_top=n147cross.select_C_stars(1.38,0.55,0.85,18)
                #selection_viz=n147vizcross.select_C_stars(1.38,0.55,0.85,18)
                

            #error printed if neither 'm' nor 'c' is chosen
                
            else:
                print('Not a valid class')
            
            self.top_subset=selection_top
            self.viz_subset=selection_viz
        
        #different selection parameters for different galaxy
        
        elif self.galaxy=='ngc185':
        
            if agb=='m':
                #parameters chosen using paper and visual inspection of cmd
                #should also verify with isochrones and other red giant branch finding methods
                selection_top = n147cross.select_M_stars(0.96,1.27,0.42,0.8,17.8)
                
            elif agb=='c':
                
                selection_top = n147cross.select_C_stars(1.27,0.42,0.8,17.8)
                
            self.top_subset=selection_top
            
        elif self.galaxy=='ngc205':
            
            if agb=='m':
                
                selection_top=n147cross.select_M_stars(0.97,1.44,0.52,0.8,18)
            
            elif agb=='c':
                
                selection_top=n147cross.select_C_stars(1.44,0.52,0.8,18)
                
            self.top_subset=selection_top
            
            
        
        #extinction carried out after subsets chosen
        
        #n147cross.sbsextinction_on_frame(selection)
        
        #subset set as class attribute
        

        
    #def rgb(self,xrange,yrange):
    
    def rgb(self):
        
        n147cross=self.n147cross
        n147cross.sbsextinction()
        
        selection=n147cross.select_rgb_region(0.975,1.12,17.75,18.3)
        self.subset=selection
        
    def crossed(self):
        n147cross=self.n147cross
        n147cross.sbsextinction()
        
        self.subset=n147cross
        

        