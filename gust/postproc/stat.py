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
    names=["time","CFL","drag","lift","moment","x","y","phi"])
df.drop(df.index[:3], inplace=True)
#
# -- plot PDF pages
figSize = (10,4)
with PdfPages('history.pdf') as pdf:
    df.plot(x='time',y=["drag","lift","moment"],figsize=figSize)
    plt.xlabel(r'$t/T$', fontsize=12)
    pdf.savefig()
    df.plot(x='time',y=["x","y","phi"],figsize=figSize)
    plt.xlabel(r'$t/T$', fontsize=12)
    pdf.savefig()
#
# ---    end of lotusStat_simple.py     --- #
