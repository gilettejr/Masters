#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 21:10:12 2019

@author: cameronrobertson
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from graphing_class import basic_graphs
from asymphr import asymphr
from celluloid import Camera

class import_isos:
    
    
    def __init__(self,galaxy,isofile):
        
        def apparent(M,dguess):
            m = M + 5*np.log10(dguess/10)
            return m
        
        self.galaxy = galaxy
        
        if self.galaxy=='ngc147':
            
            a=asymphr('lot_n147.unique',0)
            
            a.loadascii()
            a.ciscuts()
            a.sbsextinction()
            self.asymph=a
            
            
            self.distance=780000
            
        elif self.galaxy=='ngc185':
            
            self.distance= 630000
            
        
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
        
        self.iso=i
        
    def overlay_agb_tip(self):
        
        #truncates isochrones to only include AGB
        
        iso=self.iso
        
        for i in range(len(iso.jmag)):
            
            #label=8 is the AGB phase in padova isochrones
            
            if iso.label[i]!=8:

                iso.jmag[i]=np.nan
                iso.hmag[i]=np.nan
                iso.kmag[i]=np.nan
        
        #same method for seperating out the different isochrones in the set
        #as in self.isoplot
        
        indices=[]
        for i in range(1,len(iso.age)):
            if iso.age[i]!=iso.age[i-1] or iso.z[i]!=iso.z[i-1]:
                indices.append(i)
 
                
        
        

        
        plt.plot(iso.jmag[:indices[0]]-iso.kmag[:indices[0]],iso.kmag[:indices[0]],label='log(t) = ' + str(iso.age[0]) + ', Z = ' +str(iso.z[0]))
        
        for i in range(len(indices)):
            if i==(len(indices)-1):
                plt.plot(iso.jmag[indices[i]:]-iso.kmag[indices[i]:],iso.kmag[indices[i]:],label='log(t) = ' + str(iso.age[indices[i]]) + ', Z = ' +str(iso.z[indices[i]]))
                break
            else:
                plt.plot(iso.jmag[indices[i]:indices[i+1]]-iso.kmag[indices[i]:indices[i+1]],iso.kmag[indices[i]:indices[i+1]],label='log(t) = ' + str(iso.age[indices[i]]) + ', Z = ' +str(iso.z[indices[i]]))
        
        #plt.gca().invert_yaxis()
        plt.ylabel('$K_0$')
        plt.xlabel('$(J-K)_0$')
        plt.legend()
        
        #plt.show()
        #plotter=graphs()
        #plotter.isoplot(self.n147,i)
        
    def create_animation(self):
        
        asymph=self.asymph
        
        iso=self.iso
        
        fig=plt.figure()

        
        
        camera=Camera(fig)
        
        for i in range(len(iso.jmag)):
            
            if iso.label[i]!=8:
                iso.jmag[i]=np.nan
                iso.hmag[i]=np.nan
                iso.kmag[i]=np.nan
                
        indices=[]
        for i in range(1,len(iso.age)):
            if iso.age[i]!=iso.age[i-1] or iso.z[i]!=iso.z[i-1]:
                indices.append(i)
 
                
        
        
        
        
        t = plt.plot(iso.jmag[:indices[0]]-iso.kmag[:indices[0]],iso.kmag[:indices[0]],color='blue')
        plt.legend(t,['log(t) = ' + str(iso.age[0]) + ', Z = ' +str(iso.z[0])])
        
        camera.snap()
        
        for i in range(len(indices)):
            if i==(len(indices)-1):
                t = plt.plot(iso.jmag[indices[i]:]-iso.kmag[indices[i]:],iso.kmag[indices[i]:],color='blue')
                plt.legend(t,['log(t) = ' + str(iso.age[indices[i]]) + ', Z = ' +str(iso.z[indices[i]])])
                
                camera.snap()
                break
            else:
                t = plt.plot(iso.jmag[indices[i]:indices[i+1]]-iso.kmag[indices[i]:indices[i+1]],iso.kmag[indices[i]:indices[i+1]],color='blue')
                plt.legend(t,['log(t) = ' + str(iso.age[indices[i]]) + ', Z = ' +str(iso.z[indices[i]])])
                
                camera.snap()
                
        plt.scatter(asymph.jmag-asymph.kmag,asymph.kmag,marker='.',color='black')
        
        plt.ylim(12,20)
        plt.xlim(0.5,5)
        
        plt.gca().invert_yaxis()
        
        plt.ylabel('$K_0$')
        plt.xlabel('$J_0$-$K_0$')
        
        print('Saving')
        animation=camera.animate()
        
        animation.save('initial_iso.mp4')