#!/usr/bin/env python3
# ----------------------------------------- #
# stat.py
# ----------------------------------------- #
#
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

Aratio = 4
beta0 = 0.25
L_De = np.sqrt(Aratio)/(2*beta0)
coeff = 0.71
def model(t,dm):
    return coeff*(-np.sin(t)*abs(np.sin(t))
                  -np.cos(t)/(2./3.*dm/100*Aratio))
def peak(dm):
    return -min(model(np.linspace(0,2*np.pi),dm))

def str_rnd(num,d=4): return str(round(num,d))
def read_it(folder):
    try:
        df = pd.read_csv(folder+'/fort.9',delim_whitespace = True,
            names=["time","CFL","drag","lift","side"])
        df[['drag','lift','side']] *= 4 # quarter domain
    except:
        return pd.DataFrame()
    return df.query('time > 1')

def plot_it(var,axis_label,list):
    plt.figure(figsize=(6,3))
    for i,dm in enumerate(list):
        folder = 'dm'+str(dm)
        df = read_it(folder)
        Sratio = dm/100.*Aratio*4/3*L_De
        plt.plot(df.time/(2*np.pi),df[var]*Aratio/2.,label=str_rnd(Sratio,0))
        plt.plot(df.time/(2*np.pi),model(df.time,dm),'--',label='',color='C'+str(i))
    plt.xlabel(r'$t/T$', fontsize=12)
    plt.ylabel(axis_label, rotation=0, fontsize=12)
    plt.legend(title=r'$S/D_e$'); plt.tight_layout()
    pdf.savefig()
    plt.close()

# Sweep through the simulations
with PdfPages('compare_runs.pdf') as pdf:
    plot_it('drag',r'$C_D$',[10,50])

    plt.figure(figsize=(5,3))
    for dm in [5,10,15,20,25,30,40,50,60]:
        folder = 'dm'+'{:02d}'.format(dm)
        df = read_it(folder)
        Sratio = dm/100.*Aratio*4/3*L_De
        plt.scatter(Sratio,-df.drag.min()*Aratio/2.,color='C0')
    dm = np.linspace(5,60)
    plt.plot(dm/100.*Aratio*4/3*L_De,[peak(d) for d in dm],'--',color='C0')

    plt.xlabel(r'$S/D_e$')
    plt.ylabel(r'$\max(C_T)$')
    plt.ylim(0,);plt.xlim(0,); plt.tight_layout()
    plt.legend(['Bernoulli','measured'])
    pdf.savefig()
    plt.close()
