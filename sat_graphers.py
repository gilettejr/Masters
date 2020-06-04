from graphing_class import graphs
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
class sat(graphs):

    def __init__(self,galaxy,pandas=False,spitzer=False):
        
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

        if pandas == True:
            
            for i in infilenames:
                if galaxy==i:
                    pinfile=('pandas/' + i + '_pand' )
                    break
                pdata=pd.read_parquet(pinfile)
                self.pdata=pdata
        
        self.galaxy=galaxy
        self.data=data
        self.infilenames=infilenames
    def plot_tan(self,stars='all',marker='o',markersize=1,color='black'):
        
        data=self.data
        
        plt.rc('axes',labelsize = 15)
        plt.plot(data.xi,data.eta,linestyle='none',marker = marker,markersize=markersize,color=color)

        plt.gca().invert_xaxis()
        plt.gca().set_ylabel(r'$\eta$/degrees')
        plt.gca().set_xlabel(r'$\xi$/degrees')
    
    def plot_panda_cmd(self,marker='o',markersize=1,color='black'):
        
        pdata=self.pdata
        
        plt.rc('axes',labelsize = 15)
        plt.plot(pdata.g-pdata.i,pdata.i,linestyle='none',marker=marker,markersize=markersize,color=color)
        
        plt.gca().invert_yaxis()
        plt.ylabel('$i_0$')
        plt.xlabel('$g_0$-$i_0$')

                    
                
            
            
        
