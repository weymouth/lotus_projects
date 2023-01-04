import numpy as np


if __name__ == "__main__":
    T=4
    u_m = 0.4
    omega=2*np.pi/T
    ds = np.sqrt(2*1e-6/omega)
    re_ds = u_m*ds/1e-6
    re_a = re_ds**2/2
    a = u_m/omega
    st = 2*a*1/T/u_m
    print(omega, re_a, a, st)