#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import numpy


def new_f90(n: int, re: float, k_lam: int, zeta: float):
    with open("lotus.f90","r") as fileSource:
        fileLines = fileSource.readlines()

    fileLines[11] = f"    real,parameter     :: Re = {re}\n"

    if n==2:
        fileLines[19] = f"    real, parameter    :: A = 0.1*L, St_d = 0.3, k_x=0., k_z={float(k_lam)}, h_roughness=0.01\n"
        fileLines[15] = f"    integer            :: b(3) = [16,16,1]\n"
    elif n==3:
        fileLines[19] = f"    real, parameter    :: A = 0.1*L, St_d = 0.3, k_x={k_lam+0.5}, k_z={float(k_lam)}, h_roughness=0.01\n"
        fileLines[15] = f"    integer            :: b(3) = [16,16,2]\n"

    fileLines[23] = f"                          k_coeff = {1/zeta}, &\n"
    fileLines[27] = f"    integer            :: n(3), ndims={n}\n"

    # setup grid
    if k_lam % 2 ==0:
        fileLines[54] = f"      if(ndims==3) xg(3)%h = 4.\n"
    elif k_lam % 2 != 0:
        # If \k_lambda is odd we need to change the grid spacing so we don't have a half bump on the periodic boundary (L/8*delta_z = (k_lam/2+0.5)*L/k_lam)
        fileLines[54] = f"      if(ndims==3) xg(3)%h = {4*(1+1/k_lam)}\n"

    with open("lotus.f90","w") as fileOutput:
        fileOutput.writelines(fileLines)


