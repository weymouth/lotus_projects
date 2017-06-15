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
        names=["time","CFL","drag","lift","moment","y","phi"])
except FileNotFoundError:
    exit('stat: fort.9 not found')
df.drop(df.index[:3], inplace=True)

try:
    mg = pd.read_csv('fort.8',delim_whitespace = True,
        names=["itr","res0","res","inf"])[2:]
except FileNotFoundError:
    exit('stat: fort.8 not found')
#
# -- plot PDF pages
figSize = (10,4)
with PdfPages('history.pdf') as pdf:
    df.plot(x='time',figsize=figSize)
    plt.xlabel(r'$t/T$', fontsize=12)
    pdf.savefig()

    mg.plot(y=['res0','res','inf'],figsize=figSize)
    plt.yscale('log')
    pdf.savefig()

    mg.plot(y='itr',figsize=figSize)
    pdf.savefig()
#
# ---    end of lotusStat_simple.py     --- #
