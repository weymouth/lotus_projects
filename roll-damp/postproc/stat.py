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
    names=["time","CFL","moment","phi"])
df.drop(df.index[:3], inplace=True)
#
# -- define cycle, phase
df['cycle'] = np.floor(df.time)
df['phase'] = df.time-df.cycle
#
# -- remove early (transient) and last (incomplete) cycles
df.query('cycle > .0', inplace=True)
maxCycle = max(df.cycle)
df.query('cycle < @maxCycle', inplace=True)
maxCycle = max(df.cycle)
#
# -- group by cycle
grouped = df.groupby('cycle')
#
# -- plot PDF pages
figSize = (10,4)
with PdfPages('history.pdf') as pdf:
    df.plot(x='time',y='moment',figsize=figSize)
    plt.xlabel(r'$t/T$', fontsize=12)
    plt.ylabel(r'$C_M$', fontsize=12)
    pdf.savefig()
    plt.close()
#
    fig, ax = plt.subplots(figsize=figSize)
    for cycle, group in grouped:
        if cycle>maxCycle-10:
            group.plot(x='phase',y='moment',ax=ax,label=cycle)
    plt.xlabel(r'$t/T$', fontsize=12)
    plt.ylabel(r'$C_M$', fontsize=12)
    pdf.savefig()
    plt.close()

#
# ---    end of lotusStat_simple.py     --- #
