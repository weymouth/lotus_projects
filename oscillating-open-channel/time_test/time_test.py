import numpy as np
import os
import matplotlib.pyplot as plt

def read_forces(force_file, interest):
    names = ['t', 'dt', 'E', 'tke']
    # print(np.genfromtxt(force_file))
    # fos = np.transpose(np.genfromtxt(force_file))

    # forces_dic = dict(zip(names, fos))
    # t = forces_dic['t']
    # u = forces_dic[interest]

    return np.genfromtxt('/home/masseyj/Workspace/lotus_projects/oscillating-open-channel/time_test/long/fort.9')


def get_enstrophy(folder):
    t, u = read_forces(os.path.join(folder, 'fort.9'), interest='E')
    return np.mean(u[t > 4])


def get_tke(folder):
    t, u = read_forces(os.path.join(folder, 'fort.9'), interest='tke')
    return np.mean(u[t > 4])

def plot_ts():
    fig, ax = plt.subplots(figsize=(4,3))
    ax. set_xlabel(r"$t$")
    ax. set_xlabel(r"$E$")

    ax.plot(get_enstrophy(os.path.join(cwd, "long")))
    plt.savefig('ts.pdf', dpi=100)


if __name__ == "__main__":
    cwd = os.getcwd()
    folder = os.path.join(cwd, "long")
    read_forces(os.path.join(folder, 'fort.9'), 'E')

