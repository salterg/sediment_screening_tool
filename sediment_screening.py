# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 14:01:03 2023

@author: gsalter
"""

import pandas as pd
import numpy as np
from itertools import compress

all_data=pd.read_csv('sample_time_series_single_policy.csv')
paria_data=pd.read_csv('par_sed.csv',header=None).values
novpar=.9*np.squeeze(np.sort(paria_data[0:-2,0],axis=0))
aprpar=.9*np.squeeze(np.sort(paria_data[0:-1,1],axis=0))

trace=np.unique(all_data.trace_number)
hydrology=np.unique(all_data.hydrology)
demand=np.unique(all_data.demand)
initial_condition=np.unique(all_data.initial_condition)

run_matrix=np.zeros((len(trace)*len(hydrology)*len(demand)*len(initial_condition),4),dtype=np.int8)
counter=0
for j1 in range(0,len(trace)):
    for j2 in range(0,len(hydrology)):
        for j3 in range(0,len(demand)):
            for j4 in range(0,len(initial_condition)):
                run_matrix[counter,0]=j1
                run_matrix[counter,1]=j2
                run_matrix[counter,2]=j3
                run_matrix[counter,3]=j4
                counter=counter+1

for j in range(0,np.size(run_matrix)):
    run_true= (all_data.trace_number==trace[run_matrix[j,0]]) & (all_data.hydrology==hydrology[run_matrix[j,1]]) & (all_data.demand==demand[run_matrix[j,2]]) & (all_data.initial_condition==initial_condition[run_matrix[j,3]])
    idx=list(compress(range(len(run_true)),run_true))
    variable=all_data.time_series_slot_name.iloc[idx]
    idx2=list(compress(range(len(variable)), variable=='Powell Outflow'))
    idx3=list(compress(range(len(variable)), variable=='Powell Monthly Pool Elevation'))
    
    t=all_data.timestep.iloc[idx].iloc[idx2]
    Q_month=all_data.trace_value.iloc[idx].iloc[idx2].to_numpy()
    Powell_elev=all_data.trace_value.iloc[idx].iloc[idx3].to_numpy()
    
    if j==0:
        export_Tmt_allruns=np.zeros((len(run_matrix),len(t)))


    t=pd.to_datetime(t)
    Qmax_cfs=np.zeros(np.size(Q_month))
    Qmin_cfs=np.zeros(np.size(Q_month))
    MD=t.dt.daysinmonth.to_numpy()
    M=t.dt.month.to_numpy()
    rampM={1:.009,2:.009,3:.009,4:.009,5:.009, 6:.01,7:.01,8:.01, 9:.009,10:.009,11:.009,12:.009}
    for i in range(len(Q_month)):
        dailyrampcon=np.min([8000,np.multiply(Q_month[i],rampM[M[i]]),np.max([0,2*(25000-Q_month[i]/MD[i]/24/.082644)]),np.max([0,2*(-5000+Q_month[i]/MD[i]/24/.082644)])])
        A=np.array([[12, 12] , [1, -1]])
        b=np.array([[Q_month[i]/MD[i]/.082644],[dailyrampcon]])
        X=np.linalg.lstsq(A,b,rcond=None)[0]
        Qmax_cfs[i]=X[0]
        Qmin_cfs[i]=X[1]
        
    #%%
    D=(1/1000)*2**np.arange(-3.875,1,.25)   
    F_bed=np.array([0.0093293 , 0.017285,  0.030448  , 0.047813 ,  0.064944 ,  0.078534  , 0.089295   , 0.096769, \
           0.099388,  0.096762 , 0.089068 , 0.077232,  0.063017 , 0.048413 , 0.035042  , 0.0239, \
           0.015355 , 0.0092881 , 0.005287  , 0.0028307])
        
    Qmax_cms=np.vstack(Qmax_cfs)*0.028316847
    Qmin_cms=np.vstack(Qmin_cfs)*0.028316847
    Ck_vol_max= (4.0000e-27*(Qmax_cms<=707.9212)+((707.9212**4)/(707.9212**1.7))*4.0000e-27*(Qmax_cms>707.9212))*F_bed*((Qmax_cms<=707.9212)*Qmax_cms**4+(Qmax_cms>707.9212)*Qmax_cms**1.7)*(D**-3)
    Ck_vol_min= (4.0000e-27*(Qmin_cms<=707.9212)+((707.9212**4)/(707.9212**1.7))*4.0000e-27*(Qmin_cms>707.9212))*F_bed*((Qmin_cms<=707.9212)*Qmin_cms**4+(Qmin_cms>707.9212)*Qmin_cms**1.7)*(D**-3)
    Qsk_cms=.5*Ck_vol_max*Qmax_cms+.5*Ck_vol_min*Qmin_cms
    Qs_cms=np.sum(Qsk_cms,1)
    export_Tmt=np.cumsum(Qs_cms*MD*24*60*60*2650/1000/1000)
    export_Tmt_allruns[j,:]=export_Tmt
    
    #%%
    octs=list(compress(range(len(M)),M==10))
    novs=list(compress(range(len(M)),M==11))
    juls=list(compress(range(len(M)),M==7))
    mars=list(compress(range(len(M)),M==3))
    aprs=list(compress(range(len(M)),M==4))
    decs=list(compress(range(len(M)),M==12))
    
    HFE_fall_sed=np.zeros(len(novs))
    HFE_fall_implement=np.zeros(len(novs))
    HFE_spring_sed=np.zeros(len(aprs))
    HFE_spring_implement=np.zeros(len(aprs))
    
    
    if j==0:
        t_nov=t.iloc[novs]
        HFE_sed_allruns=np.zeros((len(run_matrix),len(t_nov)))
        HFE_implement_allruns=np.zeros((len(run_matrix),len(t_nov)))

        
    for i in range(0,len(novs)):
        HFE_fall_sed[i]=np.max((0,1-np.where(np.concatenate((novpar,[10**20]),axis=0)>(export_Tmt[novs[i]]-export_Tmt[juls[i]]+294))[0][0]/len(novpar)))
        if Powell_elev[octs[i]]<3550:
            HFE_fall_implement[i]=0
        else:
            HFE_fall_implement[i]=HFE_fall_sed[i]
            
    for i in range(1,len(aprs)):
        HFE_spring_sed[i]=np.max((0,1-np.where(np.concatenate((aprpar,[10**20]),axis=0)>(export_Tmt[aprs[i]]-export_Tmt[decs[i-1]]+294))[0][0]/len(aprpar)))
        if Powell_elev[mars[i]]<3525:
            HFE_spring_implement[i]=0
        else:
            HFE_spring_implement[i]=HFE_spring_sed[i]
            
    HFE_sed_allruns[j,:]=HFE_fall_sed+(1-HFE_fall_sed)*HFE_spring_sed
    HFE_implement_allruns[j,:]=HFE_fall_implement+(1-HFE_fall_implement)*HFE_spring_implement

    print(j/np.size(run_matrix))