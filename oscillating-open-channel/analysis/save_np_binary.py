#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import numpy as np

from lotusvis.flow_field import ReadIn
from lotusvis.decompositions import Decompositions


def save_phase_avg(c):
    cwd = os.getcwd()
    try:
        fns = os.listdir(f"{c}")
        if 'phase_average.npy' in fns:
            print(f'Has a np binary')
        else:
            print(f'Saving phase_average.npy for 3D')
            phase_av = Decompositions(f"{c}", "fluid", length_scale=c).phase_average(3)
            np.save(f"{cwd}/phase_average.npy", phase_av)
            print(f'Saved phase_average.npy for 3D')
            del phase_av
    except FileNotFoundError:
        print(f'DONT EXIST')


def save_phase_avg_2d():
    try:
        fns = os.listdir(f"{cwd}/2D")
        if 'phase_average.npy' in fns:
            print(f'2D has a np binary')
        else:
            print(f'Saving phase_average.npy for 2D')
            phase_av = Decompositions(f"{cwd}/2D", "fluid", length_scale=1024).phase_average(3)
            np.save(f"{cwd}/2D/phase_average.npy", phase_av)
            print(f'Saved phase_average.npy for 2D')
            del phase_av
    except FileNotFoundError:
        print(f'2D DONT EXIST')


def main():
    save_phase_avg(4096)
    # save_phase_avg_2d()


if __name__ == "__main__":
    cwd = os.getcwd()
    main()
