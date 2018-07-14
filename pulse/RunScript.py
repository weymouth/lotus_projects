#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import lotus
from itertools import product

for dm in [5,15,25,30,50,60]:
    folder1 = 'dm{:02d}_U100'.format(dm)
    lotus.replace('template.f90',
        {'DM_IN':str(dm/100),'F_IN':'.FALSE.','PER_IN':'2'})
    lotus.run(16,folder1,'../towU100/')

    folderf = 'dm{:02d}_free'.format(dm)
    lotus.replace('template.f90',
        {'DM_IN':str(dm/100),'F_IN':'.TRUE.','PER_IN':'8'})
    lotus.run(16,folderf,'../'+folder1+'/')
