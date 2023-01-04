#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import numpy


def new_f90(n: int, t, c, k, lambda_x=0., lambda_z=0.):
    with open("lotus.f90","r") as fileSource:
        fileLines = fileSource.readlines()

    fileLines[14] = f"    real,parameter     :: c={float(c)}, nu=c/Re\n"
    fileLines[15] = f"    real, parameter    :: finish={t}\n"
    fileLines[20] = f"    real, parameter    :: A = 0.1*c, St_d = 0.3, k_x={lambda_x}, k_z={lambda_z}, h_roughness=0.0\n"
    fileLines[24] = f"                          k_coeff = {k}, &\n"
    fileLines[28] = f"    integer            :: n(3), ndims={n}\n"
        
    with open("lotus.f90","w") as fileOutput:
        fileOutput.writelines(fileLines)


