#!/usr/bin/env python3
# ----------------------------------------- #
#
import numpy as np
import lotus
import pandas as pd
from scipy.optimize import brent

def y_amp(folder):
    try:
        df = pd.read_csv(folder+'/fort.9',delim_whitespace = True,
            names=["time","CFL","y","fx","fy","fz","px","py","pz","fp","pp"])
    except FileNotFoundError:
        exit('stat: '+folder+'/fort.9 not found')

    return df.query('time > 5').y.mad()*1.5748

def func(x):
    lotus.replace('resonant.f90',{'MASS':'1.0','ADDED':str(x),'MODE':'3'})
    func.count += 1
    folder = 'eval_'+str(func.count)
    lotus.run(0,folder)
    y = y_amp(folder)
    print('{}: x={}, y={}'.format(func.count,x,y))
    return -y # negative to maximize amplitude

func.count = 0
res = brent(func, brack=(0,1), tol=1E-5, full_output=1, maxiter=20)
print(res)
