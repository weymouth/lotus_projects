#!/usr/bin/env python3
# ----------------------------------------- #
# stat.py
# ----------------------------------------- #
#
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from itertools import product

def read_it(keys,values):
    folder = '_'.join([i+'_'+j for i,j in zip(keys,values)])
    try:
        df = pd.read_csv(folder+'/fort.9',delim_whitespace = True,
            names=["time","CFL","fx","fy","moment","t","x","y","phi"])
    except:
        return pd.DataFrame()

    df['drag'] = np.cos(df.phi)*df.fx+np.sin(df.phi)*df.fy
    df['lift'] = np.cos(df.phi)*df.fy-np.sin(df.phi)*df.fx
    df['u'] = (df.x.shift()-df.x.shift(-1))/(df.t.shift()-df.t.shift(-1))
    df['v'] = (df.y.shift()-df.y.shift(-1))/(df.t.shift()-df.t.shift(-1))
    df['mag'] = np.sqrt(df.u**2+df.v**2)
    df['alpha'] = np.arctan2(df.v,df.u)+df.phi

    df.drop(df.index[:3], inplace=True)

    for i,j in zip(keys,values):
        df[i] = j
    return df

# Sweep through the simulations
keys = ['RUN_FLAG','K_VAL','ALPHA_VAL']
run_flags = ['heave','center','gauss','edges','gust']
k_vals = [0.25,1.0]
alpha_vals = [15,30]

#
# plotting function
def plot_it(var,axis_label):
    for k, alpha in zip(k_vals,alpha_vals):
            fig, ax = plt.subplots(figsize=(8,4))
            for run in run_flags:
                df = read_it(keys,[run,str(k),str(alpha)])
                if(not df.empty):
                    df.query('time>=0',inplace=True)
                    if(run=='gust'):
                        df.plot(x="time",y=var,ax=ax,label=run,style='k')
                    else:
                        df.plot(x="time",y=var,ax=ax,label=run,style='--')
            plt.xlabel(r'$t/T$', fontsize=12)
            plt.ylabel(axis_label, fontsize=12)
            plt.title(r'$k=$'+str(k)+r', $\alpha$='+str(alpha),fontsize=12)
            plt.legend(title='run type')
            pdf.savefig()
            plt.close(fig)

with PdfPages('compare_runs.pdf') as pdf:
    plot_it('lift',r'$C_L$')
    plot_it('moment',r'$C_M$')
    plot_it('drag',r'$C_D$')
#
# ---    end of lotusStat_simple.py     --- #
