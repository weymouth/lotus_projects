#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import lotus
from itertools import product

# run across the mass fraction
for dm in [5,10,15,20,25,30,40,50,60]:
    folder = 'dm'+str(dm)
    lotus.replace('template.f90',{'DM_IN':str(dm/100)})
    lotus.run(16,folder)
