#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 15:14:03 2019

@author: cameronrobertson
"""

from runners import run_both,run_rgb,run_cross,kde_separator
from graphing_class import basic_graphs,graphs
from crossmatching_utils import topcatstuff
from iso_utils import import_isos
#class for reading in and plotting ascii file data

##185 centre: 9.7415417, 48.3373778
##147 centre: 8.3005, 48.5087389
##205 centre : 10.09189356, 41.68541564

#andromeda, roughly 8 deg field

#for kde method, start at 18, go upwards in bins of mag=0.4

def main():  
        
    
    
    #g=kde_separator('ngc205')
    #g.kde_graph_test()
    
    e=run_both('ngc205')
    e.c_over_m_split_plot(100,1000)
    

    

    

    
#0.25
    
    

    
main()


#radial gradiant: C/M = 0.22 <70'', C/M = 0.26 >70''

#ho Nhung paper