#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import lotus
import numpy as np
import pandas as pd

def get_power(i,j,k,Theta,r,alpha,test=False):
    #
    # set up and run
    lotus.replace('template.f90',
        {'THETAVAL':str(Theta),'RVAL':str(r),'ALPHAVAL':str(alpha)})
    folder = 'case{}{}{}'.format(i,j,k)
    lotus.run(16,folder,test=test)
    #
    # read data and drop unwanted rows and columns
    if not test:
        df = pd.read_csv(folder+'/fort.9',delim_whitespace = True,
            names=["time","CFL","power","phi"])
        df.drop(df.index[:3], inplace=True)
        #
        # -- get average and amplitude of power
        bar_p = np.trapz(df.power,df.time)/(df.time.iloc[-1]-df.time.iloc[0])
        amp_p = 1.5748*df.power.mad()
    else:
        bar_p,amp_p = 0,0
    #
    # return dict
    return {'folder':folder,'Theta':Theta,'r':r,'alpha':alpha,
                'bar_p':bar_p,'amp_p':amp_p}

Theta_list = np.logspace(-1.3,-0.3,10)
r_list = np.linspace(0.1,1,10)
alpha_list = np.pi*np.linspace(0.05,0.5,10)

df = pd.DataFrame([
    get_power(i,j,k,Theta,r,alpha)
    for i,Theta in enumerate(Theta_list)
    for j,r in enumerate(r_list)
    for k,alpha in enumerate(alpha_list)])
df.to_csv('Script_out.csv')
print(df)
