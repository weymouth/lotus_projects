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
        names=["time","x","q","y"])
except FileNotFoundError:
    exit('stat: fort.10 not found')
#
# -- post process signals
df['x'] = df['fx']+df['px']
df['y'] = df['fy']+df['py']
df['p'] = df['fp']+df['pp']
df.query('time > 0.1', inplace=True)
late = df.query('time > 5')

ln['cycle'] = np.floor(ln.time)
ln['phase'] = ((np.rint((ln.time-ln.cycle)*8)-1) % 8) +1
maxcycle = ln.cycle.max()

first = ln.query('cycle == 0').groupby('phase')
last = ln.query('cycle == @maxcycle-1').groupby('phase')
averaged = (ln.query('cycle > @maxcycle/2')
            .groupby(['phase','x'], as_index=False).mean() # take mean for each phase and x
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

def plot_phase(pdf,dat,name,label):
    fig, ax = plt.subplots(figsize=(8,4))
    for phase, group in dat:
        arg = 2.*np.pi*phase/8.
        cols = (1/2-np.cos(arg)/2,0.125,1/2-np.sin(arg)/2)
        group.plot(x='x',y=name,ax=ax,label=r'$\phi=$'+'{}/8'.format(int(phase)),c=cols)
    plt.xlabel(r'$x/c$', fontsize=12)
    plt.ylabel(label, fontsize=12)
    ax.legend(ncol=2)
    pdf.savefig()

with PdfPages('history.pdf') as pdf:
    plot_hist(pdf,name='fx',label=r'$C_D$')
    plot_hist(pdf,name='fy',label=r'$C_L$')
    plot_hist(pdf,name='p',label=r'$C_P$')

    plot_phase(pdf,first,name='q',label=r'first $ q$')
    plot_phase(pdf,first,name='y',label=r'first $y/c$')

    plot_phase(pdf,last,name='q',label=r'last $ q$')
    plot_phase(pdf,last,name='y',label=r'last $y/c$')

    plot_phase(pdf,averaged,name='q',label=r'$\widetilde q$')
    plot_phase(pdf,averaged,name='y',label=r'$\widetilde y$')

    mg.plot(y=['res0','res','inf'],figsize=(8,4))
    plt.yscale('log')
    pdf.savefig()

    mg.plot(y='itr',figsize=(8,4))
    pdf.savefig()
