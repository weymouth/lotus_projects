#!/usr/bin/env python3
# ----------------------------------------- #
# stat.py
# ----------------------------------------- #
#
import numpy as np
import pandas as pd
from scipy.integrate import trapz
import matplotlib.pyplot as plt

def average(df,var):
    return trapz(df[var],df.time)/(df.time.iloc[-1]-df.time.iloc[0])

def print_means(df,var):
    l = [average(group,var) for cycle,group in df.groupby('cycle')]
    print(var,', '.join('{:.3g}'.format(i) for i in l))

def read_it(folder,trim=True):
    try:
        df = pd.read_csv(folder+'/fort.9',delim_whitespace = True,
            names=["time","CFL","drag","lift","side","power","U"])
        df['cycle'] = np.floor(df.time/(2*np.pi))
        if trim:
            df.query('cycle<{}'.format(max(df.cycle)), inplace=True)
    except:
        exit('stat: '+folder+'/fort.9 not found')
    return df

for folder in ['dm10_U100','dm10_free','dm40_U100','dm40_free']:
    print(folder)
    df = read_it(folder)
    [print_means(df,var) for var in ['drag','power','U']]
