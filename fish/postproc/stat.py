#!/usr/bin/env python3
# ----------------------------------------- #
# stat.py
# ----------------------------------------- #
#
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
#
# -- read data
try:
    df = pd.read_csv('fort.9',delim_whitespace = True,
        names=["time","CFL","fx","fy","fz","px","py","pz","fp","pp"])
except FileNotFoundError:
    exit('stat: fort.9 not found')
try:
    mg = pd.read_csv('fort.8',delim_whitespace = True,
        names=["itr","res0","res","inf"])[2:]
except FileNotFoundError:
    exit('stat: fort.8 not found')
try:
    ln = pd.read_csv('fort.10',delim_whitespace = True,
        names=["time","x","py"])
except FileNotFoundError:
    exit('stat: fort.10 not found')
#
# -- post process signals
df['x'] = df['fx']+df['px']
df['y'] = df['fy']+df['py']
df['p'] = df['fp']+df['pp']

ln['cycle'] = np.floor(ln.time)
ln['phase'] = np.rint((ln.time-ln.cycle)*8)
#
# -- remove early (transient) and make late
df.query('time > 0.1', inplace=True)
late = df.query('time > 5')
ln.query('time > 5', inplace=True)
#
# -- group the lines
grouped = (ln.groupby(['phase','x'], as_index=False).mean() # take mean for each phase and x
            .groupby('phase')) # regroup by phase
#
# -- plot PDF pages
def plot_hist(pdf,name,label):
    ax = df.plot(x='time',y=name,figsize=(8,4))
    plt.xlabel(r'$t/T$', fontsize=12)
    plt.ylabel(label, fontsize=12)
    txt = 'mean='+str(round(late[name].mean(),4)) \
        +', mad='+str(round(1.5748*late[name].mad(),4))
    plt.text(0.5,0.01,txt,transform=ax.transAxes)
    pdf.savefig()

with PdfPages('history.pdf') as pdf:
    plot_hist(pdf,name='x',label=r'$C_D$')
    # plot_hist(pdf,name='y',label=r'$C_L$')
    plot_hist(pdf,name='p',label=r'$C_P$')
    # plot_hist(pdf,name='fx',label=r'$C_{Df}$')
    # plot_hist(pdf,name='fy',label=r'$C_{Lf}$')
    # plot_hist(pdf,name='fp',label=r'$C_{Pf}$')
    # plot_hist(pdf,name='px',label=r'$C_{Dp}$')
    # plot_hist(pdf,name='py',label=r'$C_{Lp}$')
    # plot_hist(pdf,name='pp',label=r'$C_{Pp}$')
    # plot_hist(pdf,name='CFL',label='CFL')

    fig, ax = plt.subplots(figsize=(8,4))
    for name, group in grouped:
        arg = 2.*np.pi*name/8.
        cols = (1/2-np.cos(arg)/2,0.125,1/2-np.sin(arg)/2)
        group.plot(x='x',y='py',ax=ax,label=r'$\phi=$'+'{}/8'.format(int(name)),c=cols)
    plt.xlabel(r'$x/c$', fontsize=12)
    plt.ylabel(r'$\widetilde C_{Y}$', fontsize=12)
    ax.legend(ncol=2, bbox_to_anchor=(0.4, 1.05))
    pdf.savefig()

    mg.plot(y=['res0','res','inf'],figsize=(8,4))
    plt.yscale('log')
    pdf.savefig()

    mg.plot(y='itr',figsize=(8,4))
    pdf.savefig()
