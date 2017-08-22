#!/usr/bin/env python3
# ----------------------------------------- #
#
import numpy as np
import lotus
import pandas as pd
from scipy.optimize import brent
from functools import lru_cache

def y_amp(folder):
    try:
        df = pd.read_csv(folder+'/fort.9',delim_whitespace = True,
            names=["time","CFL","y","fx","fy","fz","px","py","pz","fp","pp"])
    except FileNotFoundError:
        exit('stat: '+folder+'/fort.9 not found')
    return df.query('time > 5').y.mad()*1.5748

@lru_cache(maxsize=None)
def func(x):
    lotus.replace('resonant.f90',{'MASS':'1.0','ADDED':str(x),'MODE':'1'})
    func.count += 1
    folder = 'eval_'+str(func.count)
    lotus.run(16,folder)
    y = y_amp(folder)
    print('{}: x={}, y={}'.format(func.count,x,y))
    return -y # negative to maximize amplitude

def summary(thedir):
    import os
    dic = {}
    for dirs in (d for d in os.listdir(thedir) if os.path.isdir(os.path.join(thedir, d))):
        fulldir = os.path.join(thedir,dirs)
        x = read_added(fulldir)
        y = y_amp(fulldir)
        dic[x] = y

    f = open(thedir+'/summary.txt','w')
    print(thedir)
    for key in sorted(dic.keys()):
        print('{}  {}'.format(key, dic[key]))
        print('{}  {}'.format(key, dic[key]),file=f)
    f.close()

def read_added(folder):
    f = open(folder+'/lotus.f90','r')
    for text in f:
        if 'ma = ' in text:
            for word in text.split(','):
                if 'ma = ' in word:
                    return float(word.split('ma = ',1)[1])
    f.close()

func.count = 0
res = brent(func, brack=(0.12,0.32,0.72), tol=1E-3, full_output=1, maxiter=20)
print(res)
