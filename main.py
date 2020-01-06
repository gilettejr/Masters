#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 15:14:03 2019

@author: cameronrobertson
"""

from runners import run_both,run_rgb,run_cross,kde_separator
from graphing_class import basic_graphs
from crossmatching_utils import topcatstuff
#class for reading in and plotting ascii file data

##185 centre: 9.7415417, 48.3373778
##147 centre: 8.3005, 48.5087389

#andromeda, roughly 8 deg field

#for kde method, start at 18, go upwards in bins of mag=0.4

def main():  
        

    r=kde_separator('ngc147')
    r.kde_graph_test()
    
    

    
main()


#radial gradiant: C/M = 0.22 <70'', C/M = 0.26 >70''

#ho Nhung paper