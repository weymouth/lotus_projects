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

def cycle_mean(df,var):
    return [average(group,var) for cycle,group in df.groupby('cycle')]
def print_means(df,var):
    print(var,', '.join('{:.3g}'.format(i) for i in cycle_mean(df,var)))

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

# df = read_it('AoAe4_dm20/free',trim=False)
# df[["drag","lift","side","power"]] *= 4/1.25 # correct the area-scale
# df[["time","CFL","drag","lift","side","power","U"]].to_csv(r'dm20_free/fort.9', header=None, index=None, sep=' ', mode='a')

R100 = cycle_mean(read_it('towU100'),'drag')[-1]
R112 = cycle_mean(read_it('towU112'),'drag')[-1]/(1.12)**2
def RU(U):
    r = (U-1)/(0.12)
    return (R100*(1-r)+R112*r)*U**3

# for dm in [5,10,15,20,25,30,40,50,60]:
#     for folder in ['dm{:02d}_'.format(dm)+s for s in ['U100','free']]:
#         print(folder)
#         df = read_it(folder)
#         [print_means(df,var) for var in ['drag','power','U']]

dm = [5,10,15,20,25,30,40]
