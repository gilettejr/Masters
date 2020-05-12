#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 13:54:07 2020

@author: cameronrobertson
"""

from andro import m31data

import numpy as np

class andro_utils(m31data):
    
    def make_c_m(self,fore,jklim,hklim,jhlim,trgb):
        
        data=self.data
        mdata=data.copy()
        cdata=data.copy()
        
        jk=data.jmag-data.kmag
        hk=data.hmag-data.kmag
        jh=data.jmag-data.hmag
        
        for i in range(len(mdata.kmag)):
            
            if np.isnan(mdata.kmag)==True:
                
                continue
            
            elif hk[i]>hklim or jk[i]>jklim or jh[i]>jhlim or mdata[i].kmag>trgb or jk<fore:
                mdata.loc[i]=np.nan
                
        for i in range(len(cdata.kmag)):
            
            if np.isnan(cdata.kmag)==True:
                
                continue
            
            elif hk[i]<hklim or jk[i]<jklim or jh[i]<jhlim or mdata[i].kmag>trgb:
                cdata.loc[i]=np.nan
                
        self.cdata=cdata
        self.mdata=mdata
        
    def find_FEH(self):
            
        def CM_to_FEH(CM):
        
            FEH=-1.39 -0.47*np.log10(CM)
        
            return(FEH)
        
        mdata=self.mdata
        cdata=self.cdata
        
        m_no_top=[]
        for i in mdata.kmag:
            if np.isnan(i)==False:
                m_no_top.append(0)
        c_no_top=[]   
        for j in cdata.kmag:
            if np.isnan(j)==False:
                c_no_top.append(0)
        
        C_top=len(c_no_top)
        M_top=len(m_no_top)
        
        ratio_top = C_top/M_top
        FEH_top=CM_to_FEH(ratio_top)
        
        
        print('Gaia crossmatch: C/M = ' + str(ratio_top) + ' and [Fe/H] = ' + str(FEH_top))
    
    
