from sat_graphers import sat
import numpy as np
import pandas as pd
class ratio_utils(sat):
    
    def __init__(self,*args,**kwargs):
        super(ratio_utils,self).__init__(*args,**kwargs)
        
    
    
        jhc=np.array([0,0,0,10])
        hkc=np.array([0,0,0,8])
        jkc=np.array([0,0,0,4])
        jk4c=np.array([0,0,0,2])
        trgb=([0,0,0,18])
        
        cuts=np.array([jhc,hkc,jkc,jk4c,trgb])
        cut_frame=[]
        for i in range(len(cuts)):
            print(i)
            for j in range(len(cuts[i])):
                if self.infilenames[j]==self.galaxy:
                    cut_frame.append(cuts[i][j])
                    break
        cut_frame=np.array(cut_frame)
        
        jhc=cut_frame[0]
        hkc=cut_frame[1]
        jkc=cut_frame[2]
        jk4c=cut_frame[3]
        trgb=cut_frame[4]