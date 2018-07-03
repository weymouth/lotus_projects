#!/usr/bin/env python3
# ----------------------------------------- #
# stat.py
# ----------------------------------------- #
#
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def read_it(folder):
    try:
        df = pd.read_csv(folder+'/fort.9',delim_whitespace = True,
            names=["time","CFL","drag","lift","side"])
        df[2:] *= 4 # quarter domain
#         df.time = df.time/(2.*np.pi)
    except:
        return pd.DataFrame()
    return df #.query('time > 1')

def plot_it(var,axis_label,list):
    fig, ax = plt.subplots(figsize=(8,4))
    for folder in list:
        df = read_it(folder)
        df.plot(x="time",y=var,ax=ax,label=folder)
    plt.xlabel(r'$t/T$', fontsize=12)
    plt.ylabel(axis_label, fontsize=12)
    plt.legend(title='folder')
    pdf.savefig()
    plt.close()

# Sweep through the simulations
Aratio = 4
beta0 = 0.25
L_De = np.sqrt(Aratio)/(2*beta0)
print(L_De)
list = ['dm'+str(dm) for dm in ['05',10,15,20,25,30,40,50,60]]
with PdfPages('compare_runs.pdf') as pdf:
    plot_it('drag',r'$C_D$',list)
    # plot_it('lift',r'$C_L$',list)

    for dm,folder in zip([5,10,15,20,25,30,40,50,60],list):
        df = read_it(folder)
        Sratio = dm/100.*Aratio*4/3*L_De
        plt.scatter(Sratio,-df.drag.mean()*Aratio/2.,color='C0')
        plt.scatter(Sratio,-df.drag.min()*Aratio/2.,color='C1')
    plt.xlabel(r'$S/D_e$')
    plt.ylabel(r'$\frac{T}{\dot m U_j}$',rotation=0,size=14,labelpad=14)
    plt.ylim(0,);plt.xlim(0,)
    pdf.savefig()
