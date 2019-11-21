#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 21:10:12 2019

@author: cameronrobertson
"""
import numpy as np
import pandas as pd
from graphing_class import basic_graphs
from asymphr import asymphr

class import_isos(basic_graphs):
    
    
    def read_in(self,isofile):
        
        def apparent(M,dguess):
            m = M + 5*np.log10(dguess/10)
            return m
        
        if self.galaxy=='ngc147':
            self.distance=780000
        
        distance=self.distance
        
        i=asymphr(isofile,0)
        
        isos=pd.DataFrame({'z':i.data[:,0],'age':i.data[:,2],'label':i.data[:,9],'jmag':i.data[:,27],'hmag':i.data[:,28],'kmag':i.data[:,29]})
        
        isos.jmag=apparent(isos.jmag,distance)
        isos.hmag=apparent(isos.hmag,distance)
        isos.kmag=apparent(isos.kmag,distance)
        
        i.z=isos.z
        i.age=isos.age
        i.label=isos.label
        i.jmag=isos.jmag
        i.hmag=isos.hmag
        i.kmag=isos.kmag
        
        self.i=i
        #plotter=graphs()
        #plotter.isoplot(self.n147,i)