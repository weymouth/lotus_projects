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
# read data and drop unwanted rows and columns
try:
    df = pd.read_csv('fort.9',delim_whitespace = True,
        names=["time","drag","lift","side"])
except FileNotFoundError:
    exit('stat: fort.9 not found')
#
# -- plot PDF pages
def str_rnd(num,d=4): return str(round(num,d))

def plot_hist(pdf,name,label):
    ax = df.plot(x='time',y=name,figsize=(8,4))
    plt.xlabel(r'$t/T$', fontsize=12)
    plt.ylabel(label, fontsize=12)
    mean,mad = df[name].mean(), 1.5748*df[name].mad()
    x1,x2,y1,y2 = plt.axis()
    mx,mn = min(y2,mean+3*mad),max(y1,mean-3*mad)
    plt.ylim([mn,mx])
    txt = 'mean='+str_rnd(mean)+', mad='+str_rnd(mad)
    plt.text(0.5,0.01,txt,transform=ax.transAxes)
    pdf.savefig()
    plt.close()

with PdfPages('history.pdf') as pdf:
    plot_hist(pdf,name='drag',label=r'$C_{Xp}$')
    plot_hist(pdf,name='lift',label=r'$C_{Yp}$')
    plot_hist(pdf,name='side',label=r'$C_{Zp}$')
