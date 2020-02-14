#!/usr/bin/env python3
# ----------------------------------------- #
# stat.py
# ----------------------------------------- #
#
import numpy as np
import pandas as pd
from scipy.integrate import trapz
import matplotlib.pyplot as plt

def ave(df,var):
    return trapz(df[var],df.time)/(df.time.iloc[-1]-df.time.iloc[0])

def amp(df,var):
    return 1.5748*df[var].mad()

def cycle_ave(df,var):
    return [ave(group,var) for cycle,group in df.groupby('cycle')]

def fmt(list,s='{:.3g}'):
    return ', '.join(s.format(i) for i in list)

def read_it(folder,trim=True,last=False):
    try:
        df = pd.read_csv(folder+'/fort.9',delim_whitespace = True,
            names=["time","CFL","drag","lift","side","power","U"])
        df['cycle'] = np.floor(df.time/(2*np.pi))
        df['phase'] = df.time/(2*np.pi)-df.cycle
        if trim:
            df.query('cycle<{}'.format(max(df.cycle)), inplace=True)
            if(len(df)==0):
                exit('stat: no full cycles')
        if last:
            df.query('cycle=={}'.format(df.cycle.max()), inplace=True)
    except:
        exit('stat: '+folder+'/fort.9 not found')
    return df

beta0,dm = 0.25,0.2

df = pd.DataFrame({'AoAe':[20,25,30,40,50,60]})
df['AR'] = df.AoAe/10.*(170*0.25+5)**2/(170*0.25)**2
df['folder'] = ['AoAe{}_U000'.format(A) for A in df.AoAe]
for var in ['drag','power']:
    df[var+'0'] = [A*ave(read_it(f,last=True),var) for A,f in zip(df.AR,df.folder)]
    df[var+'1'] = [A*amp(read_it(f,last=True),var) for A,f in zip(df.AR,df.folder)]
# df['Vf'] = 2./3.+np.sqrt(1-10./df.AoAe)-(1-10./df.AoAe)**(1.5)/3.
# df['Per'] = df.Vf*dm/2.*df.AoAe/10.
# df['SR'] = df.Per/beta0*np.sqrt(df.AoAe/10.)
# df['eta0'] = -df.drag0/df.power0
df.to_csv('AoAe_U000.csv')


df = pd.DataFrame({'AoAe':[20,30,40,50,60]})
df['folder'] = ['AoAe{}_U100'.format(A) for A in df.AoAe]
for var in ['drag','power']:
    df[var+'0'] = [ave(read_it(f,last=True),var) for f in df.folder]
    df[var+'1'] = [amp(read_it(f,last=True),var) for f in df.folder]
print(df)
df.to_csv('AoAe_U100.csv')
