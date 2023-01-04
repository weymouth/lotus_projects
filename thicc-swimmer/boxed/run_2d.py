#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from lotus import run
from changef90 import new_f90
from pathlib import Path

def extract_k():
    with open("lotus.f90","r") as fileSource:
        fileLines = fileSource.readlines()
    txt = fileLines[24]
    return float([s for s in txt.split(' ')][-2][:-1])


def run2d(k_lam):
    k = extract_k()
    new_f90(2, 10000, k_lam, float(1./k))
    run(256, f'{cwd}/smooth')


if __name__ == "__main__":
    cwd = Path.cwd()
    k_lam = 4
    run2d(k_lam)
    

