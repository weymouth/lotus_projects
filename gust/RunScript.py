#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import subprocess
from time import sleep
from itertools import product
import os

def run_lotus(keys,values):
    "Write a lotus.f90 file and run in a folder using a set of keys and values"
    folder = '_'.join([i+'_'+j for i,j in zip(keys,values)]) # name the folder
    if(not os.path.isdir(folder)):
        print('Going to run: '+folder)
        write_lotus(dict(zip(keys,values)))                      # make the lotus file
        subprocess.check_call('runLotus 16 '+folder, shell=True)  # run the simulation

def write_lotus(dic):
    "Write the lotus file by replacing from the template"
    f1 = open('template.f90','r')
    f2 = open('lotus.f90','w')
    for line in f1:
        f2.write(replace_all(line,dic))
    f1.close()
    f2.close()
    return

def replace_all(text, dic):
    "Replace all dictionary keys in text with their value"
    for i, j in dic.items():
        text = text.replace(i, j)
    return text

# Set up the keys and values
keys = ['RUN_FLAG','K_VAL','ALPHA_VAL']
run_flags = ['gust','edges','center','heave']
k_vals = [0.25,0.5,0.75,1.0]
alpha_vals = [15,30]
# Sweep through the simulations
for k,alpha,run in product(k_vals,alpha_vals,run_flags):
    run_lotus(keys,[run,str(k),str(alpha)])
