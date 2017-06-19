#!/usr/bin/env python3
# ----------------------------------------- #
# stat.py
# ----------------------------------------- #
#
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def read_it(keys,values):
    folder = '_'.join([i+'_'+j for i,j in zip(keys,values)])
    try:
        df = pd.read_csv(folder+'/fort.9',delim_whitespace = True,
            names=["time","CFL","drag","lift","moment","x","y","phi"])
    except:
        df = pd.read_csv(folder+'/fort.9',delim_whitespace = True,
            names=["time","CFL","drag","lift","moment","y","phi"])
    df.drop(df.index[:3], inplace=True)
    for i,j in zip(keys,values):
        df[i] = j
    return df

keys = ['K_VAL','KINEMATIC_FLAG','ALPHA_VAL']
k_vals = [0.125,0.25,0.5,1.0]
kinematic_flags = ['true','false']
alpha_vals = [15]

# Sweep through the simulations
df = pd.DataFrame()
for k in k_vals:
    for kin in kinematic_flags:
        for alpha in alpha_vals:
            df = pd.concat([df,read_it(keys,[str(k),kin,str(alpha)])])
df.query('time>=0',inplace=True)
df = df.groupby(['K_VAL','KINEMATIC_FLAG'])
#
# -- plot PDF pages
figSize = (8,4)
with PdfPages('compare_k.pdf') as pdf:
    for name,group in df:
        print(name)
        if(name[1]=='false'):
            fig, ax = plt.subplots(figsize=figSize)
        group.plot(x="time",y="lift",ax=ax,label=name)
        plt.xlabel(r'$t/T$', fontsize=12)
        plt.ylabel(r'$C_L$', fontsize=12)
        if(name[1]=='true'):
            pdf.savefig()
#
# ---    end of lotusStat_simple.py     --- #
