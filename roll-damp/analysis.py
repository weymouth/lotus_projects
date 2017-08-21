#!/usr/bin/env python3
# ----------------------------------------- #
# stat.py
# ----------------------------------------- #
#
import numpy as np
import pandas as pd

def savefig(ax,name='plot.pdf'):
    fig = ax.get_figure()
    fig.savefig(name)

def read_from(folder):
    df = pd.read_csv(folder+'/fort.9',delim_whitespace = True,
        names=["time","CFL","moment","phi"])
    df.drop('CFL', axis=1, inplace=True)
    df.query('time >= 3', inplace=True)
    return df

def resample_evenly(df):
    'resample a DataFrame onto an evenly spaced time series'
    n = len(df.time)
    even_df = pd.DataFrame()
    even_df['time'] = np.linspace(df.time.min(),df.time.max(),n)
    for name in df:
        if(name == 'time'): continue
        even_df[name] = np.interp(even_df.time,df.time,df[name])
    return even_df

def take_fft(df):
    'take the fft of all non-time columns in a DataFrame'
    n = len(df.phi)
    T = df.time.max()-df.time.min()
    hat_df = pd.DataFrame()
    hat_df['freq'] = np.linspace(0,n/T,n)
    for name in df:
        if(name == 'time'): continue
        hat_df[name] = np.fft.fft(df[name])*2./n
    hat_df.query('freq < 3.0', inplace=True)
    return hat_df

def find_max_row(df):
    'return the DataFrame row with the maximum phi value'
    return df.ix[abs(df.phi).idxmax()]

def get_coeffs(folder):
    df = (read_from(folder)
            .pipe(resample_evenly)
            .pipe(take_fft)
            .pipe(find_max_row))
    print(folder,':')
    print(df)

# ax = (read_from('3D_D_142')
#     .pipe(resample_evenly)
#     .pipe(take_fft)
#     .plot(x='freq'))
# savefig(ax,'temp.pdf')

for folder in ['3D','2D_D_71','smag_0p4_D_71','smag_0p8_D_71']:
    get_coeffs(folder)

# for num in ['34','48','68','96','136','192']:
#     get_coeffs('smag_D_'+num)
