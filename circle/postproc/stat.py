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
        names=["time","CFL","fx","fy","fz","px","py","pz"])
except FileNotFoundError:
    exit('stat: fort.9 not found')
try:
    mg = pd.read_csv('fort.8',delim_whitespace = True,
        names=["itr","res0","res","inf"])[2:]
except FileNotFoundError:
    exit('stat: fort.8 not found')
#
# -- post process signals
df['x'] = 1.*(df['fx']+df['px'])
df['y'] = 1.*(df['fy']+df['py'])
df['z'] = 0.
#
# -- remove early (transient) and make late
df.query('time > 2', inplace=True)
meanTime = df.time.mean()
late = df.query('time > @meanTime')
#
# -- plot PDF pages
def plot_hist(pdf,name,label):
    ax = df.plot(x='time',y=name,figsize=(8,4))
    plt.xlabel(r'$tU/L$', fontsize=12)
    plt.ylabel(label, fontsize=12)
    txt = 'mean='+str(round(late[name].mean(),4)) \
        +', mad='+str(round(1.5748*late[name].mad(),4))
    plt.text(0.5,0.01,txt,transform=ax.transAxes)
    pdf.savefig()

with PdfPages('history.pdf') as pdf:
    plot_hist(pdf,name='x',label=r'$C_D$')
    plot_hist(pdf,name='y',label=r'$C_L$')
    plot_hist(pdf,name='fx',label=r'$C_{Df}$')
    plot_hist(pdf,name='fy',label=r'$C_{Lf}$')
    plot_hist(pdf,name='px',label=r'$C_{Dp}$')
    plot_hist(pdf,name='py',label=r'$C_{Lp}$')
    plot_hist(pdf,name='CFL',label='CFL')

    mg.plot(y=['res0','res','inf'],figsize=(8,4))
    plt.yscale('log')
    pdf.savefig()

    mg.plot(y='itr',figsize=(8,4))
    pdf.savefig()
