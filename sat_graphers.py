from graphing_class import graphs
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
class sat(graphs):

    def __init__(self,galaxy):
        
        def quadrat(a,b):
            
            c=np.sqrt(a**2+b**2)
            
            return c
        self.mag=quadrat
        
        infilenames=['ngc147','ngc185','ngc205','m32']
        
        for i in infilenames:
            if galaxy==i:
                infile=('pm33_4rem/' + i)
                break
        
        data=pd.read_parquet(infile)
        
        self.galaxy=galaxy
        self.data=data
    
    def plot_tan(self,stars='all',marker='o',markersize=1,color='black'):
        
        data=self.data
        
        plt.rc('axes',labelsize = 15)
        plt.plot(data.xi,data.eta,linestyle='none',marker = marker,markersize=markersize,color=color)

        plt.gca().invert_xaxis()
        plt.gca().set_ylabel(r'$\eta$/degrees')
        plt.gca().set_xlabel(r'$\xi$/degrees')
    
    def plot_hist_slices(self,maxm=20,minm=16,slicemag=0.3):
        
        slices=np.linspace(minm,maxm,int((maxm-minm)/slicemag),endpoint=False)

        colours=[]
        kmags=[]
        tdata=self.data.copy()
        
        for i in slices:
            
            data=tdata.copy()
            

            length=[]
            for j in range(len(data.kmag)):

                
                if np.isnan(data.kmag[j])==False:
                    length.append(0)
                    
                
                if data.kmag[j] > (i + slicemag) or data.kmag[j] < i:

                    data.loc[j]=np.nan
                    #data.jmag[j]=np.nan
                    
            print(len(length))
            

            colours.append(data.jmag-data.kmag)
            kmags.append(data.kmag)

        fig,axs=plt.subplots(2,1)
            
        axs[0].hist(colours[1].dropna())
        axs[1].plot(colours[1],kmags[1],linestyle='none',marker='o',markersize='1',color='black')
        axs[1].plot(colours[1],kmags[1],linestyle='none',marker='o',markersize='1',color='red')
                    
                    
                    
                    
                
                
                
            
            
        
