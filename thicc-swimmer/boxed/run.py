#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from lotus import run
import numpy as np
from changef90 import new_f90
from pathlib import Path
from scipy.optimize import brentq
import os
import sys
from lotusvis.plot_flow import *
from lotusvis.flow_field import *


def read_forces(force_file, c, interest='p'):
    names = ['t', 'dt', 'px', 'py', 'pz', 'cpx', 'cpy', 'cpz', 'vx', 'vy', 'vz', 'Ex', 'Ey', 'Ez']
    fos = np.transpose(np.genfromtxt(force_file))

    forces_dic = dict(zip(names, fos))
    t = forces_dic['t']
    u = np.array((forces_dic[interest + 'x'], forces_dic[interest + 'y']))

    ux, uy = np.squeeze(np.array(u[0])), np.squeeze(np.array(u[1]))

    return t, ux, uy


def get_thrust(folder, c):
    t, ux, _ = read_forces(os.path.join(folder, 'fort.9'), c)
    t, vx, _ = read_forces(os.path.join(folder, 'fort.9'), c, interest='v')
    return np.mean(ux[t > 4]) + np.mean(vx[t > 4])


def optimisation_func(zeta):
    new_f90(3, 10000, k_lam, zeta)
    run(512, f'{cwd}/rough')
    thrust = get_thrust(f'{cwd}/rough', 1024)
    return thrust


def optimise():
    k_0 = brentq(optimisation_func, 0.8, 2.0,
                 xtol=1e-2, maxiter=7)


if __name__ == "__main__":
    cwd = Path.cwd()
    k_lam = 16
    sol = optimise()
    print(sol)
    

