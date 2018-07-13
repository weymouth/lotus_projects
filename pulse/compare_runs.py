#!/usr/bin/env python3
# ----------------------------------------- #
# stat.py
# ----------------------------------------- #
#
import numpy as np
import pandas as pd
from scipy.integrate import trapz
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

Aratio = 4
beta0 = 0.25
L_De = np.sqrt(Aratio)/(2*beta0)
coeff = 1
def model(t,dm):
    return coeff*(np.sin(t)*abs(np.sin(t))
                  +np.cos(t)/(2./3.*dm/100*Aratio))
def peak(dm):
    return -min(model(np.linspace(0,np.pi),dm))

def mean(df_in,var,T):
    df = df_in.query('time>{}'.format(T))
    return trapz(df[var],df.time)/(df.time.iloc[-1]-df.time.iloc[0])

def str_rnd(num,d=4): return str(round(num,d))
def read_it(folder):
    try:
        df = pd.read_csv(folder+'/fort.9',delim_whitespace = True,
            names=["time","CFL","drag","lift","side","power"])[3:]
        df[['drag','lift','side','power']] *= 4 # quarter domain
    except:
        return pd.DataFrame()
    return df

def plot_hist(var,axis_label,list):
    plt.figure(figsize=(6,3))
    for i,dm in enumerate(list):
        folder = 'dm'+str(dm)
        df = read_it(folder)
        Sratio = dm/100.*Aratio*4/3*L_De
        plt.plot(df.time/(2*np.pi),df[var]*Aratio/2.,label=str_rnd(Sratio,1))
        plt.plot(df.time/(2*np.pi),model(df.time,dm),'--',label='',color='C'+str(i))
    plt.xlabel(r'$t/T$', fontsize=12)
    plt.ylabel(axis_label, rotation=0, fontsize=12)
    plt.ylim(-5,5)
    plt.legend(title=r'$S/D_e$'); plt.tight_layout()
    pdf.savefig()
    plt.close()

def plot_means(fnc,axis_label,list):
    plt.figure(figsize=(5,3))
    for dm in list:
        folder = 'dm'+str(dm)
        df = read_it(folder)
        Sratio = dm/100.*Aratio*4/3*L_De
        plt.scatter(Sratio,fnc(df,2*np.pi),color='C0')
        plt.scatter(Sratio,fnc(df,3*np.pi),color='C1')
    plt.xlabel(r'$S/D_e$'); plt.ylim(0,)
    plt.ylabel(r'$\overline{'+axis_label+r'}$');
    plt.tight_layout();plt.legend(['full cycle','jet only'])
    pdf.savefig()
    plt.close()

# Sweep through the simulations
with PdfPages('compare_runs.pdf') as pdf:
    plot_hist('drag',r'$C_D$',[10,50])

    list = [5,10,15,20,25,30,40,50,60]
    plot_means(lambda df,T:-mean(df,'drag',T)*Aratio/2.,'C_T',list)
    plot_means(lambda df,T: mean(df,'power',T)*Aratio/2.,'C_P',list)
    plot_means(lambda df,T:-mean(df,'drag',T)/mean(df,'power',T),'\eta',list)
