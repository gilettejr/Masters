#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 15:14:03 2019

@author: cameronrobertson
"""

from runners import run_both,run_rgb
from graphing_class import basic_graphs
#class for reading in and plotting ascii file data



        
        
def main():  
        
    #r=run_both('ngc185')
    #r.c_over_m_grad(70,9.7415417,48.3373778)
    #r.CM_to_FEH(r.inCM)
    #r.CM_to_FEH(r.outCM)
    
    f=run_both('ngc147')
    f.c_over_m_grad(70,8.3005,48.5087389)
    f.CM_to_FEH(f.inCM)
    f.CM_to_FEH(f.outCM)



main()


#radial gradiant: C/M = 0.22 <70'', C/M = 0.26 >70''