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

def fmt(list,s='{:.3g}'):
    return ', '.join(s.format(i) for i in list)

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

R100 = cycle_mean(read_it('towU100'),'drag')[-1]*4/1.25
R112 = cycle_mean(read_it('towU112'),'drag')[-1]/(1.12)**2*4/1.25
def RU(U):
    r = (U-1)/(0.12)
    return (R100*(1-r)+R112*r)*U**3

df = pd.DataFrame({'dm':[5,10,15,20,25,30,40,50]})
df['folder'] = ['dm{:02d}_free'.format(dm) for dm in df.dm]
for var in ['drag','power','U']:
    df[var] = [cycle_mean(read_it(f),var)[-1] for f in df.folder]
df['eta'] = RU(df.U)/df.power
print(df)
df.to_csv('free_swimming_AoAe4.csv')
