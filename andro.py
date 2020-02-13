#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 15:31:10 2020

@author: cameronrobertson
"""

from crossmatching_utils import topcatcross,crossed_data
import pandas as pd

#takes in m31 file as input, applies cls and extinction cuts, saves to csv file
def clsextcut():
    
    f=topcatcross('lot_m31.unique',0)
        
    f.loadascii()
    f.ciscuts()
    f.sbsextinction()
    f.create_tangent_coords(10.6847083,41.26875)
    f.make_dataframe()
    
    data=f.wfdat
    
    data.to_parquet('m31_int_files/int1')

def topcrossed():
    rmatch=crossed_data()
    
    n147topcross = topcatcross('lot_m31.unique',0)
    rmatch.topmatch(n147topcross,'crossedm31.csv')
    n147topcross.delete_crossed_points()
    
    n147topcross.create_tangent_coords(10.6847083,41.26875)
    n147topcross.make_dataframe()
    
    data=n147topcross.wfdat
    
    print(data)
    
    data.to_parquet('m31_int_files/int2')


topcrossed()

    

class m31data:
    
    #class for andromeda methods
    def __init__(self,intno):
        
        data=pd.read_parquet('m31_int_files/int' + str(intno))
        self.data=data

#m31=m31data(1)

        