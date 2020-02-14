#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import lotus
import os
from itertools import product
import pandas as pd

def runcase(folder,dic,resume=None,n=16):
    print(dic)
    lotus.replace('template.f90',dic)
    lotus.run(n,folder,resume,test=False)

parameters = {'DM_IN':'0.2','AR_IN':'4.0','PER_IN':'3.0',
            'F_IN':'.FALSE.','V_IN':'0.0','a_IN':'0.0','DP_IN':'4'}

df = pd.DataFrame({
    'AR':[2.0,2.0,2.0,2.0,3.0,3.0,3.0,3.0,3.0,3.0,4.0,4.0,4.0,4.0,4.0,4.0],
    'dm':[0.05,0.25,0.45,0.65,0.05,0.15,0.25,0.35,0.45,0.55,0.05,0.1,0.15,0.2,.25,.3]
})
for AR,dm in df.itertuples(index=False):
    p = parameters.copy()
    p['AR_IN'] = str(AR)
    p['DM_IN'] = str(dm)
    folder = 'AoAe{:02d}_dm{:02d}_U000'.format(int(AR*10),int(dm*100))
    if not os.path.exists(folder+'/fort.9'):
        runcase(folder,p)
