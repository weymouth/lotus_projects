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
        return pd.DataFrame()
    df.drop(df.index[:3], inplace=True)
    for i,j in zip(keys,values):
        df[i] = j
    return df

# Sweep through the simulations
keys = ['RUN_FLAG','K_VAL','ALPHA_VAL']
run_flags = ['gust','edges','center','heave']
k_vals = [0.5,1.0]
alpha_vals = [15,30]
# Sweep through the simulations
for run,k,alpha in product(run_flags,k_vals,alpha_vals):

#
# plotting function
def plot_it(var,axis_label):
    for k, alpha in product(k_vals,alpha_vals):
            fig, ax = plt.subplots(figsize=(8,4))
            for run in run_flags:
                df = read_it(keys,[run,str(k),str(alpha)])
                if(not df.empty):
                    df.query('time>=0',inplace=True)
                    df.plot(x="time",y=var,ax=ax,label=run)
            plt.xlabel(r'$t/T$', fontsize=12)
            plt.ylabel(axis_label, fontsize=12)
            plt.title(r'$k=$'+str(k)+r', $\alpha$='+str(alpha),fontsize=12)
            plt.legend(title='run type')
            pdf.savefig()

with PdfPages('compare_k.pdf') as pdf:
    plot_it('lift',r'$C_L$')
    plot_it('drag',r'$C_D$')
    plot_it('moment',r'$C_M$')
#
# ---    end of lotusStat_simple.py     --- #
