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
df = pd.read_csv('fort.9',delim_whitespace = True,
    names=["time","CFL","power","phi"])
df.drop(df.index[:3], inplace=True)

def average(df,var):
    return np.trapz(df[var],df.time)/(df.time.iloc[-1]-df.time.iloc[0])
#
# -- plot PDF pages
with PdfPages('history.pdf') as pdf:
    ax = df.plot(x='time',y='power',figsize=(10,4))
    plt.xlabel(r'$t/T$', fontsize=12)
    plt.ylabel(r'$C_P$', fontsize=12)

    mean,mad = average(df,'power'), 1.5748*df.power.mad()
    x1,x2,y1,y2 = plt.axis()
    mx,mn = min(y2,mean+3*mad),max(y1,mean-3*mad)
    plt.ylim([mn,mx])
    txt = 'mean={:.3g}, mad={:.3g}'.format(mean,mad)
    plt.text(0.5,0.01,txt,transform=ax.transAxes)

    pdf.savefig()
    plt.close()

#
# ---    end of lotusStat_simple.py     --- #
