from sat_graphers import sat
import numpy as np
import pandas as pd
class ratio_utils(sat):
    
    def __init__(self,*args,**kwargs):
        super(ratio_utils,self).__init__(*args,**kwargs)
        
        def CM_to_FEH(CM):
            FEH=-1.39 -0.47*np.log10(CM)
        
            return(FEH)            
        jhc=np.array([0.81,0.70,0.77,10])
        hkc=np.array([0.50,0.46,0.57,8])
        jkc=np.array([1.33,1.30,1.42,1.60])
        jk4m=np.array([0.99,0.988,0.96,0.86])
        trgb=([18,17.8,17.8,18])
        
        cuts=np.array([jhc,hkc,jkc,jk4m,trgb])
        cut_frame=[]
        for i in range(len(cuts)):

            for j in range(len(cuts[i])):
                if self.infilenames[j]==self.galaxy:
                    cut_frame.append(cuts[i][j])
                    break
        cut_frame=np.array(cut_frame)
        
        self.jhc=cut_frame[0]
        self.hkc=cut_frame[1]
        self.jkc=cut_frame[2]
        self.jk4m=cut_frame[3]
        self.trgb=cut_frame[4]
        self.CM_to_FEH=CM_to_FEH

        data=self.data
        hk=data.hmag-data.kmag
        jh=data.jmag-data.hmag
        jk=data.jmag-data.kmag
        mdata=data.copy()
        cdata=data.copy()


        for i in data.index:
            if data.kmag[i] > self.trgb or jk[i] < self.jk4m:
                mdata.loc[i]=np.nan
                cdata.loc[i]=np.nan
            elif (hk[i] > self.hkc and jh[i] > self.jhc)  or jk[i] > self.jkc:
                mdata.loc[i]=np.nan
                
            else:
                cdata.loc[i]=np.nan
        
        self.cdata=cdata.dropna()
        self.mdata=mdata.dropna()
    
    def define_marginals(self,jkmult=1,hkmult=1,jhmult=1):
        
        for j in range(len(self.infilenames)):
            
            if self.galaxy==self.infilenames[j]:
                k=j
                break
            
        
        jkdata=self.data.copy()
        hkdata=self.data.copy()

        jherrs=[0,0,0,0]
        hkerrs=[0,0,0,0]
        jkerrs=[0,0,0,0]
        
        hk=self.data.hmag-self.data.kmag
        jh=self.data.jmag-self.data.hmag
        jk=self.data.jmag-self.data.kmag
        
        
        
        
        jkmmdata=self.data.copy()
        jkMdata=self.data.copy()
        jkcmdata=self.data.copy()
        jkCdata=self.data.copy()
        
        hkmmdata=self.data.copy()
        hkMdata=self.data.copy()
        hkcmdata=self.data.copy()
        hkCdata=self.data.copy()
        

        
        for i in jkdata.index:
            if jkdata.kmag[i] > self.trgb or jk[i] < self.jk4m:
                jkMdata.loc[i]=np.nan
                jkmmdata.loc[i]=np.nan
                jkcmdata.loc[i]=np.nan
                jkCdata.loc[i]=np.nan
                jkdata.loc[i]=np.nan
                
                hkMdata.loc[i]=np.nan
                hkmmdata.loc[i]=np.nan
                hkcmdata.loc[i]=np.nan
                hkCdata.loc[i]=np.nan
                hkdata.loc[i]=np.nan
                
                
            elif jk[i]>self.jkc-jkerrs[k]:
                jkMdata.loc[i]=np.nan
        
        jkMdata=jkMdata.dropna()
        jkmmdata=jkmmdata.dropna()
        jkcmdata=jkcmdata.dropna()
        jkCdata=jkCdata.dropna()
        jkdata=jkdata.dropna()
        
        hkMdata=jkMdata.dropna()
        hkmmdata=jkmmdata.dropna()
        hkcmdata=jkcmdata.dropna()
        hkCdata=jkCdata.dropna()
        hkdata=jkdata.dropna()
        
        for i in jkdata.index:
            if jk[i]<self.jkc - jkerrs[k] or jk[i] > self.jkc:
                jkmmdata.loc[i]=np.nan
        jkmmdata=jkmmdata.dropna()
        
        for i in jkdata.index:
            if jk[i] < self.jkc or jk[i] > self.jkc + jkerrs[k]:
                jkcmdata.loc[i]=np.nan
        jkcmdata=jkcmdata.dropna()
        
        for i in jkdata.index:
            if jk[i] > self.jkc + jkerrs[k]:
                jkCdata.loc[i].dropna()
        jkCdata=jkCdata.dropna()
        
        
        for i in hkdata.index:
            
            if hk[i] > self.hkc-hkerrs[k] and jh[i] > self.jhc-jherrs[k]:
                hkMdata.loc[i]=np.nan
        hkMdata=hkMdata.dropna()
        
        for i in hkdata.index:
            
            if (hk[i] > self.hkc and jh[i] > self.jhc) or (hk[i] < self.hkc-hkerrs[k] or jh[i] < self.jhc-jherrs[k]):
                hkmmdata.loc[i]=np.nan
                
        hkmmdata=hkmmdata.dropna()
        
        for i in hkdata.index:
            
            if (hk[i] > self.hkc+ hkerrs[k] and jh[i] > self.jhc + jherrs[k]) or (hk[i]<self.hkc or jh[i] < self.jhc):
                hkcmdata.loc[i]=np.nan
        hkcmdata=hkcmdata.dropna()
        
        for i in hkdata.index:
            
            if hk[i] < self.hkc+hkerrs[k] or jh[i] > self.jhc+jherrs[k]:
                hkCdata.loc[i]=np.nan
        hkCdata=hkCdata.dropna()
                



            
        
        
    
    def total_CM_FEH(self):
        
        cm=np.divide(len(self.cdata),len(self.mdata))
        FEH= self.CM_to_FEH(cm)
        
        print('Average C/M ratio = ' + str(cm) +', average [Fe/H] = ' + str(FEH))
        
        