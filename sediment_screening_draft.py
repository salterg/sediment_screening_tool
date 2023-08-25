# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 14:01:03 2023

@author: gsalter
"""

import pandas as pd
import numpy as np
from itertools import compress

all_data=pd.read_csv('sample_time_series_single_policy.csv') #crss output policy file
paria_data=pd.read_csv('par_sed.csv',header=None).values #paria input data

#parse paria input distributions for use in HFE probability calculation, distributions include extra 10% for UMC lesser tributaries
novpar=.9*np.squeeze(np.sort(paria_data[0:-2,0],axis=0)) #assume 90% of load (lower-bound HFE planning purposes)
aprpar=.9*np.squeeze(np.sort(paria_data[0:-1,1],axis=0))
aprseispar=.9*np.squeeze(np.sort(paria_data[0:-2,2],axis=0)) 
mayseispar=.9*np.squeeze(np.sort(paria_data[0:-2,3],axis=0)) 
junseispar=.9*np.squeeze(np.sort(paria_data[0:-2,4],axis=0)) 

#product of four variables below gives full list of unique scenarios
trace=np.unique(all_data.trace_number)
hydrology=np.unique(all_data.hydrology)
demand=np.unique(all_data.demand)
initial_condition=np.unique(all_data.initial_condition)

#create a matrix storing values of above for all unique scenarios
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

#%% load up unique scenarios
#loop through all unique scenarios, listed in run_matrix
for j in range(0,len(run_matrix)):
    
    #run_true selects a unique scenario
    run_true= (all_data.trace_number==trace[run_matrix[j,0]]) & (all_data.hydrology==hydrology[run_matrix[j,1]]) & (all_data.demand==demand[run_matrix[j,2]]) & (all_data.initial_condition==initial_condition[run_matrix[j,3]])
    idx=list(compress(range(len(run_true)),run_true))
    variable=all_data.time_series_slot_name.iloc[idx]
    #find variables of interest for the unique scenario
    idx2=list(compress(range(len(variable)), variable=='Powell Outflow'))
    idx3=list(compress(range(len(variable)), variable=='Powell Monthly Pool Elevation'))
    
    t=all_data.timestep.iloc[idx].iloc[idx2]
    Q_month=all_data.trace_value.iloc[idx].iloc[idx2].to_numpy()
    Powell_elev=all_data.trace_value.iloc[idx].iloc[idx3].to_numpy()
    t=pd.to_datetime(t)
    

    #%% sediment mass balance metric
    
    #find max and min daily flows 
    Qmax_cfs=np.zeros(np.size(Q_month)) #initialize
    Qmin_cfs=np.zeros(np.size(Q_month))
    MD=t.dt.daysinmonth.to_numpy() #days per month
    M=t.dt.month.to_numpy() #month number (1-12)
    rampM={1:.009,2:.009,3:.009,4:.009,5:.009, 6:.01,7:.01,8:.01, 9:.009,10:.009,11:.009,12:.009} #LTEMP daily fluctuation constraint
    for i in range(len(Q_month)):
        #fluctuation based on ltemp constraints (fluctuation in cfs = 0.9 or 1.0 % of monthly release in af, max 8000, instantaneous release not to exceed 25,000 cfs or less than 5000 (if they must, flow is steady)
        dailyrampcon=np.min([8000,np.multiply(Q_month[i],rampM[M[i]]),np.max([0,2*(25000-Q_month[i]/MD[i]/24/.082644)]),np.max([0,2*(-5000+Q_month[i]/MD[i]/24/.082644)])])
        A=np.array([[12, 12] , [1, -1]]) #for simplicity assume 12 hrs high 12 hrs low
        b=np.array([[Q_month[i]/MD[i]/.082644],[dailyrampcon]]) 
        X=np.linalg.lstsq(A,b,rcond=None)[0] #find max/min flows given fluctuation constraint and monthly release volume
        Qmax_cfs[i]=X[0]
        Qmin_cfs[i]=X[1]
        
    if j==0: #initialize
        MB_annual_allruns=np.zeros((len(run_matrix),len(t[M==1]))) 
        
    D=(1/1000)*2**np.arange(-3.875,1,.25) #vector of discrete grain sizes for computing xport   
    #bed grain size distribution based on mean Wright et al. 2010 values for 2002-2023
    F_bed=np.array([0.0093293 , 0.017285,  0.030448  , 0.047813 ,  0.064944 ,  0.078534  , 0.089295   , 0.096769, \
           0.099388,  0.096762 , 0.089068 , 0.077232,  0.063017 , 0.048413 , 0.035042  , 0.0239, \
           0.015355 , 0.0092881 , 0.005287  , 0.0028307]) 
        
    Qmax_cms=np.vstack(Qmax_cfs)*0.028316847 #convert cfs to cms
    Qmin_cms=np.vstack(Qmin_cfs)*0.028316847
    #calculate sed conc based on bed grain size, max/min releases calculated above, and Wright et al. 2010 model coefficients for LMC
    Ck_vol_max= (4.0000e-27*(Qmax_cms<=707.9212)+((707.9212**4)/(707.9212**1.7))*4.0000e-27*(Qmax_cms>707.9212))*F_bed*((Qmax_cms<=707.9212)*Qmax_cms**4+(Qmax_cms>707.9212)*Qmax_cms**1.7)*(D**-3)
    Ck_vol_min= (4.0000e-27*(Qmin_cms<=707.9212)+((707.9212**4)/(707.9212**1.7))*4.0000e-27*(Qmin_cms>707.9212))*F_bed*((Qmin_cms<=707.9212)*Qmin_cms**4+(Qmin_cms>707.9212)*Qmin_cms**1.7)*(D**-3)
    Qsk_cms=.5*Ck_vol_max*Qmax_cms+.5*Ck_vol_min*Qmin_cms
    Qs_cms=np.sum(Qsk_cms,1)
    paria_avg_month=np.array([40.217,12.123,15.841,3.101,2.5485,16.846,106.3,254.36,321.26,141.15,14.348,15.192]) #avg monthly paria inputs plus 10% for UMC lesser tributaries
    export_Tmt=np.cumsum(Qs_cms*MD*24*60*60*2650/1000/1000) #convert to cumulative export in thousand metric tons, this gets used for HFE probability
    MB_Tmt=-Qs_cms*MD*24*60*60*2650/1000/1000+paria_avg_month[M-1] #MC mass balance (monthly) = monthly paria input minus monthly LMC export
    
    counter=0
    MB_annual=np.zeros(len(t[M==1]))
    t_annual=t.values[(M==1)]
    for i in range(0, len(MB_Tmt)):
        MB_annual[counter]=MB_annual[counter]+MB_Tmt[i] #find annual mass balance, i.e. end of year mass balance
        if M[i]==12: #move to next year
            counter=counter+1
            
    MB_annual_allruns[j,:]=MB_annual #save this run into matrix storing all unqiue scenarios
    
    #Potential performance metrics (probably pick one)
    #Average annual mass balance 
    #threshold: percent years with mass balance less than +294 Tmt (294 Tmt = export for 60-hr hfe)
    

    #%% HFE likelihood metric
    
    mars=list(compress(range(len(M)),M==3))
    aprs=list(compress(range(len(M)),M==4))
    mays=list(compress(range(len(M)),M==5))
    juns=list(compress(range(len(M)),M==6))
    octs=list(compress(range(len(M)),M==10))
    novs=list(compress(range(len(M)),M==11))
    
    HFE_fall_sed=np.zeros(len(novs))
    HFE_fall_implement=np.zeros(len(novs))
    HFE_spring_sed=np.zeros(len(aprs))
    HFE_spring_implement=np.zeros(len(aprs))
    
    
    if j==0:
        t_nov=t.iloc[novs]
        HFE_implement_allruns=np.zeros((len(run_matrix),len(t_nov)))

        
    for i in range(0,len(novs)):
        #obtain % of years where paria input as of nov 1 is greater than export since jul 1 plus 294 Tmt (export for a 60-hr HFE)
        HFE_fall_sed[i]=np.max((0,1-np.where(np.concatenate((novpar,[10**20]),axis=0)>(export_Tmt[novs[i]]-export_Tmt[juns[i]]+294))[0][0]/len(novpar)))
        if Powell_elev[octs[i]]<3550: #no hfe's if oct EOM elevation of powell below 3550'
            HFE_fall_implement[i]=0
        else:
            HFE_fall_implement[i]=HFE_fall_sed[i]
            
    for i in range(1,len(aprs)): #note, skip spring of first year (do not calculate hfe probability for first (partial) sediment year)
        #obtain % of years where paria input as of apr 1 is greater than export since jul 1 plus 294 Tmt (export for a 60-hr HFE)
        HFE_spring_sed[i]=np.max((0,1-np.where(np.concatenate((aprpar,[10**20]),axis=0)>(export_Tmt[juns[i]]-export_Tmt[novs[i-1]]+294))[0][0]/len(aprpar)))
        if Powell_elev[mars[i]]<3525: #no hfe's if mar EOM elevation of powell below 3525'
            HFE_spring_implement[i]=0
        else:
            HFE_spring_implement[i]=HFE_spring_sed[i]
            
    #save this run into matrix storing all unqiue scenarios   
    HFE_implement_allruns[j,:]=HFE_fall_implement+(1-HFE_fall_implement)*HFE_spring_implement
    #note, probability not calculated for first (partial) year but is calculated for final (partial) year because it includes the fall implementation window where HFE's are much more likeley
    
    #Performance metric: average HFE implementation probability across all years and traces
    #Brad Butterfield's veg model will use HFE implementation probability (probably as a binary yes/no HFE based on probability >50% for implementation, please check with him)
   
    #%% SEIS HFE likelihood metric - skip this for now unless something changes

    HFE_seis4_sed=np.zeros(len(aprs))
    HFE_seis4_implement=np.zeros(len(aprs))
    HFE_seis5_sed=np.zeros(len(mays))
    HFE_seis5_implement=np.zeros(len(mays))
    HFE_seis6_sed=np.zeros(len(juns))
    HFE_seis6_implement=np.zeros(len(juns))
    
    if j==0:
        HFE_implement_seis_allruns=np.zeros((len(run_matrix),len(t_nov)))
        
    for i in range(1,len(aprs)): #probability of april hfe implementation w/ 1-yr accounting window
        HFE_seis4_sed[i]=np.max((0,1-np.where(np.concatenate((aprseispar,[10**20]),axis=0)>(export_Tmt[aprs[i]]-export_Tmt[juns[i-1]]+294))[0][0]/len(aprseispar)))
        if Powell_elev[mars[i]]<3525:
            HFE_seis4_implement[i]=0
        else:
            HFE_seis4_implement[i]=HFE_seis4_sed[i]
            
    for i in range(1,len(mays)): #probability of may hfe implementation w/ 1-yr accounting window
        HFE_seis5_sed[i]=np.max((0,1-np.where(np.concatenate((mayseispar,[10**20]),axis=0)>(export_Tmt[mays[i]]-export_Tmt[juns[i-1]]+294))[0][0]/len(mayseispar)))
        if Powell_elev[aprs[i]]<3525:
            HFE_seis5_implement[i]=0
        else:
            HFE_seis5_implement[i]=HFE_seis5_sed[i]
            
    for i in range(1,len(juns)): #probability of jun hfe implementation w/ 1-yr accounting window
        HFE_seis6_sed[i]=np.max((0,1-np.where(np.concatenate((junseispar,[10**20]),axis=0)>(export_Tmt[juns[i]]-export_Tmt[juns[i-1]]+294))[0][0]/len(junseispar)))
        if Powell_elev[mays[i]]<3525:
            HFE_seis6_implement[i]=0
        else:
            HFE_seis6_implement[i]=HFE_seis6_sed[i]
            
    #save this run into matrix storing all unqiue scenarios   
    #hfe annual probability is max of probability in any implementation month w/ 1-yr window (probabilities are not independent)
    #nov implementation is same as under LTEMP, but spring implementation probability is different due to 1-yr window     
    HFE_implement_seis_allruns[j,:]=np.max(((HFE_fall_implement),(HFE_seis4_implement),(HFE_seis5_implement),(HFE_seis6_implement)),axis=0) #max annual probability of HFE


    print(j/len(run_matrix))