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
from crossmatching_utils import topcatcross
from asymphr import asymphr
from celluloid import Camera

class import_isos:
    
    
    def __init__(self,galaxy,isofile):
        
        def apparent(M,dguess):
            m = M + 5*np.log10(dguess/10)
            return m
        
        def merr(dist,derr):
            
            m=5*np.log10(dist/10)-5*np.log10((dist-derr)/10)
            
            return m
        def mag(a,b):
            c=np.sqrt(a**2+b**2)
            return c
        
        self.merr=merr
        self.mag=mag
        
        self.galaxy = galaxy
        self.isofile=isofile
        
        if self.galaxy=='ngc147':
            
            a=topcatcross('lot_n147.unique',0)
            
        

            
            
            self.distance=680000
            self.disterr=30000
            
            
        elif self.galaxy=='ngc185':
            
            a=topcatcross('lot_n185.unique',0)
            
            self.distance = 630000
            self.disterr = 30000
        
        elif self.galaxy=='ngc205':
            
            a=topcatcross('N205_new_trimmed.unique',0)
            
            self.distance=820000
            self.disterr=30000
            
        elif self.galaxy=='andromeda':
            
            a=topcatcross('lot_n147.unique',0)
            
        
            self.distance=785000
            self.disterr=25000
        a.loadascii()
        a.ciscuts()
        a.sbsextinction()
        a.make_dataframe()
        self.dat=a
            
        
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
 
                
        
        plt.rc('axes',labelsize=15)

        
        plt.plot(iso.jmag[:indices[0]]-iso.kmag[:indices[0]],iso.kmag[:indices[0]],label='log(t) = ' + str(iso.age[0]) + ', Z = ' +str(iso.z[0]))
        
        for i in range(len(indices)):
            if i==(len(indices)-1):
                plt.plot(iso.jmag[indices[i]:]-iso.kmag[indices[i]:],iso.kmag[indices[i]:],label='log(t) = ' + str(iso.age[indices[i]]) + ', Z = ' +str(iso.z[indices[i]]))
                break
            else:
                plt.plot(iso.jmag[indices[i]:indices[i+1]]-iso.kmag[indices[i]:indices[i+1]],iso.kmag[indices[i]:indices[i+1]],label='log(t) = ' + str(iso.age[indices[i]]) + ', Z = ' +str(iso.z[indices[i]]))
        
        #plt.gca().invert_yaxis()
        plt.ylabel('$K_0$')
        plt.xlabel('$J_0$-$K_0$')
        plt.legend()
        
        #plt.show()
        #plotter=graphs()
        #plotter.isoplot(self.n147,i)
        
    def overlay_agb_single(self):
        iso=self.iso
        
        for i in range(len(iso.jmag)):
            
            #label=8 is the AGB phase in padova isochrones
            
            if iso.label[i]!=8:

                iso.jmag[i]=np.nan
                iso.hmag[i]=np.nan
                iso.kmag[i]=np.nan
                
        plt.rc('axes',labelsize=15)
        
        plt.plot(iso.jmag-iso.kmag,iso.kmag,color='green',label='log(t)= ' + str(iso.age[0])+ ', Z= ' +str(iso.z[0]))
        
        plt.ylabel('$K_0$')
        plt.xlabel('$J_0$-$K_0$')
        plt.legend()
        
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
                
        plt.scatter(dat.jmag-dat.kmag,dat.kmag,marker='.',color='black')
        
        plt.ylim(12,20)
        plt.xlim(0.5,5)
        
        plt.gca().invert_yaxis()
        
        plt.ylabel('$K_0$')
        plt.xlabel('$J_0$-$K_0$')
        
        print('Saving')
        animation=camera.animate()
        
        animation.save('initial_iso.mp4')
        
        
    def quant_age(self):

        
        data=self.dat
        
        xerror=np.sqrt((data.jerr)**2+(data.kerr)**2)
        
        yerror=self.mag(self.merr(self.distance,self.disterr),data.kerr)
        
        if self.galaxy=='ngc147':
        
            #cdat=data.select_C_stars(1.33,0.44,0.83,18)
            cdat=data.select_M_stars(1.0,1.33,0.44,0.83,18)
            
            for i in range(len(cdat.jmag)):
                
                if np.isnan(cdat.jerr[i])==False:
                    print(xerror[i])
                    print(yerror[i])
            
        elif self.galaxy=='ngc185':
            
            #cdat=data.select_C_stars(1.27,0.42,0.8,17.8)
            cdat=data.select_M_stars(0.96,1.27,0.42,0.8,17.8)
            
        elif self.galaxy=='ngc205':
            
            #cdat=data.select_C_stars(1.44,0.52,0.8,18)
            cdat=data.select_M_stars(0.97,1.44,0.52,0.8,18)
            
        
        
        
        cj=cdat.jmag
        ch=cdat.hmag
        ck=cdat.kmag
        ccol=cj-ck
        
        count=0
        
        for i in cj:
            if np.isnan(i)==False:
                
                count=count+1
                
        print(count)
        
        iso=self.iso
        
        fig,axes=plt.subplots(2)

        
        
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
 
                
        
        
        
        

        
        isoj=iso.jmag[:indices[0]].values
        isoh=iso.hmag[:indices[0]].values
        isok=iso.kmag[:indices[0]].values
        isocol=isoj-isok
        agehitlist=[]
        age_data=[]
        nhits=[]
        
        for j in range(len(cj)):
            
            if np.isnan(cj[j])==True:
                continue
            
            for k in range(len(isoj)):
                
                if cj[j] > isoj[k]-yerror[j] and cj[j] < isoj[k] + yerror[j] and ccol[j] > isocol[k]-xerror[j] and ccol[j] < isocol[k] + xerror[j]:
                    
                    nhits.append(0)
                    break
            
        
        agehitlist.append(len(nhits)-1)

        age_data.append(iso.age[0])
        
        print(len(agehitlist))
        
        t = axes[0].plot(iso.jmag[:indices[0]]-iso.kmag[:indices[0]],iso.kmag[:indices[0]],color='blue')
        g=axes[1].plot(age_data,agehitlist,color='red')
        

        
        
        axes[0].legend(t,['log(t) = ' + str(iso.age[0]) + ', Z = ' +str(iso.z[0])])
        axes[1].legend(g,['log(t) = ' + str(iso.age[0]) + ', Z = ' +str(iso.z[0])])
                    
        
        camera.snap()
        

        
        for i in range(len(indices)):
            
            nhits=[]
            

            
            if i==(len(indices)-1):

                
                isoj=iso.jmag[indices[i]:].values
                isoh=iso.hmag[indices[i]:].values
                isok=iso.kmag[indices[i]:].values
                isocol=isoj-isok
                
                for j in range(len(cj)):
                    
                    if np.isnan(cj[j])==True:
                        
                        continue
            
                    for k in range(len(isoj)):
                
                        if cj[j] > isoj[k]-yerror[j] and cj[j] < isoj[k] + yerror[j] and ccol[j] > isocol[k]-xerror[j] and ccol[j] < isocol[k] + xerror[j]:
                    
                            nhits.append(0)
                            break
                
                age_data.append(iso.age[indices[i]])
                agehitlist.append(len(nhits)-1)
                
                t = axes[0].plot(iso.jmag[indices[i]:]-iso.kmag[indices[i]:],iso.kmag[indices[i]:],color='blue')
                g=axes[1].plot(age_data,agehitlist,color='red')
                axes[0].legend(t,['log(t) = ' + str(iso.age[indices[i]]) + ', Z = ' +str(iso.z[indices[i]])])
                axes[1].legend(g,['log(t) = ' + str(iso.age[indices[i]]) + ', Z = ' +str(iso.z[indices[i]])])
                
                camera.snap()
                
                agehitlist.append(len(nhits)-1)
                
                print(len(agehitlist))
                
                break


            else:

                
                isoj=iso.jmag[indices[i]:indices[i+1]].values
                isoh=iso.hmag[indices[i]:indices[i+1]].values
                isok=iso.kmag[indices[i]:indices[i+1]].values
                isocol=isoj-isok
            
                
                for j in range(len(cj)):
                    
                    
                    if np.isnan(cj[j])==True:
                        
                        continue
            
                    for k in range(len(isoj)):

                        if cj[j] > isoj[k]-yerror[j] and cj[j] < isoj[k] + yerror[j] and ccol[j] > isocol[k]-xerror[j] and ccol[j] < isocol[k] + xerror[j]:
                    
                            nhits.append(0)
                            break

                

                age_data.append(iso.age[indices[i]])
                agehitlist.append(len(nhits)-1)
                
                print(len(agehitlist))
                
                
                t = axes[0].plot(iso.jmag[indices[i]:indices[i+1]]-iso.kmag[indices[i]:indices[i+1]],iso.kmag[indices[i]:indices[i+1]],color='blue')
                g = axes[1].plot(age_data,agehitlist,color='red')
        
        
                axes[0].legend(t,['log(t) = ' + str(iso.age[indices[i]]) + ', Z = ' +str(iso.z[indices[i]])])
                axes[1].legend(g,['log(t) = ' + str(iso.age[indices[i]]) + ', Z = ' +str(iso.z[indices[i]])])
    
        
                camera.snap()
            

                
        axes[0].scatter(data.jmag-data.kmag,data.kmag,marker='.',color='black')
        
        axes[0].set_ylim(12,20)
        axes[0].set_xlim(0,5)
        
        axes[0].invert_yaxis()
        
        axes[0].set_ylabel('$K_0$')
        axes[0].set_xlabel('$J_0$-$K_0$')
        
        print('Saving')
        animation=camera.animate()
        
        if self.isofile=='Isochrones/animation_lowz_205.dat':
        
            animation.save('iso_low_185.mp4')
            
        elif self.isofile=='Isochrones/animation_highz_205.dat':
            
            animation.save('iso_high_185.mp4')
            
        else:
            
            animation.save('nothiing.mp4')
            
    def plot_solar(self):
        plt.rc('axes',labelsize=15)
        plt.plot(self.iso.jmag-self.iso.hmag,self.iso.kmag)
        
        plt.xlabel('J-K')
        plt.ylabel('K')
        plt.gca().invert_yaxis()