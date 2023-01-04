
from lotus import run
from changef90 import new_f90
import os

def res():
    Ss = [128, 256, 512, 1024, 2046]
    Ss = [64, 128]
    for S in Ss:
        new_f90(S, k_lam)
        run(8, f'res_test/{int(S)}')
        os.chdir('..')


def xz_domains():
    S = 128
    Ls = [1, 2, 4, 8]
    for l in Ls:
        new_f90(S, k_lam, L=(l, 5, l))
        run(8, f'xz_dom_test/{int(l)}')
        os.chdir('..')


def y_domain():
    S = 128
    lxz = 2  # Use the value from the previous test
    Ls = [2, 4, 8, 16, 32]
    for l in Ls:
        new_f90(S, k_lam, L=(lxz, l, lxz))
        run(8, f'y_dom_test/{int(l)}')
        os.chdir('..')

        
def time_test():
    S = 128
    L = (2, 4, 2)  # Use the value from the previous test
    new_f90(S, k_lam, L=L, t_max=2)
    run(8, f'time_test/long')


if __name__ == "__main__":
    k_lam = 2
    time_test()

