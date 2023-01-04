#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
from turtle import color
import warnings
import vtk
import numpy as np
import seaborn as sns
from scipy.interpolate import interp1d
from inout import read_forces, get_power, FreqConv, get_power_rms
from matplotlib import pyplot as plt

from lotusvis.flow_field import ReadIn
# from lotusvis.plot_flow import Plots
from lotusvis.decompositions import Decompositions
from lotusvis.assign_props import AssignProps

import matplotlib.colors as colors
import numpy as np
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.style.use(["science", "grid"])


def plot_vortx_integral(directory, cs):
    fig, ax = plt.subplots(figsize=(4, 4))

    ax.set_ylabel(r'$ \int \omega_x $')
    ax.set_xlabel(r"$ \Delta y/h $")

    slice_int = np.zeros(len(cs))
    average_int = np.zeros(len(cs))
    total_int = np.zeros(len(cs))
    for idx, c in enumerate(cs):
        flow = FlowBase(f"{directory}/{c}", "fluid", length_scale=c)

        cell_volume = 4/c * 2/c * 4/c
        vortx = vorticity_x(flow)*cell_volume
        mag = np.sqrt((vortx)**2)

        slice= mag[ :, np.shape(vortx)[1]//2, :]
        avg = np.mean(mag, axis=1)

        slice_int[idx] = np.trapz(slice.flat)
        average_int[idx] = np.trapz(avg.flat)
        total_int[idx] = np.trapz(mag.flat)

    # slice_int = (slice_int-slice_int[-1])/slice_int
    # average_int = (average_int-average_int[-1])/average_int
    # total_int = (total_int-total_int[-1])/total_int

    ax.plot(200/np.array(cs), slice_int, marker="+", color="steelblue", label="Midplane slice")
    ax.plot(200/np.array(cs), average_int, marker="*", color="red", label="y average")
    # ax.plot(200/np.array(cs), total_int, marker=".", color="green", label="Total field")

    # ax.set_xscale("log")
    ax.loglog()
    ax.legend()
    plt.savefig(
        f"{directory}/figures/vortx_integral_convergence.pdf", dpi=300, transparent=False
    )
    plt.close()


def plot_vortx_body_integral(directory, cs):
    fig, ax = plt.subplots(figsize=(4, 4))

    ax.set_ylabel(r'$ \int \omega_x $')
    ax.set_xlabel(r"$ \Delta y/h $")

    slice_int = np.zeros(len(cs))
    average_int = np.zeros(len(cs))
    total_int = np.zeros(len(cs))
    for idx, c in enumerate(cs):
        flow = FlowBase(f"{directory}/{c}", "fluid", length_scale=c)

        cell_volume = 4/c * 2/c * 4/c
        vortx = vorticity_x(flow)*cell_volume

        # Trim to the body
        vortx = np.ma.masked_where(flow.X>1, vortx).filled(0)
        mag = np.sqrt((vortx)**2)

        slice= mag[ :, np.shape(vortx)[1]//2, :]
        avg = np.mean(mag, axis=1)

        slice_int[idx] = np.trapz(slice.flat)
        average_int[idx] = np.trapz(avg.flat)
        total_int[idx] = np.trapz(mag.flat)

    # slice_int = (slice_int-slice_int[-1])/slice_int
    # average_int = (average_int-average_int[-1])/average_int
    # total_int = (total_int-total_int[-1])/total_int

    # ax.plot(200/np.array(cs), slice_int, marker="+", color="steelblue", label="Midplane slice")
    ax.plot(200/np.array(cs), average_int, marker="*", color="red", label="y average")
    ax.plot(200/np.array(cs), total_int, marker=".", color="green", label="Total field (just body, no wake)")

    # ax.set_xscale("log")
    ax.loglog()
    ax.legend()
    plt.savefig(
        f"{directory}/figures/vortx_body_integral_convergence.pdf", dpi=300, transparent=False
    )
    plt.close()


def plot(flow, fn_save, c, **kwargs):
    cwd = os.getcwd()
    x = np.mean(flow.X, axis=0)
    y = np.mean(flow.Z, axis=0)

    mag = vorticity_x(flow)
    # mag = np.ma.masked_where(flow.X>1, mag).filled(0)
    mag = np.mean(mag, axis=0)
    

    plt.style.use(["science", "grid"])
    fig, ax = plt.subplots(figsize=(14, 4))
    divider = make_axes_locatable(ax)
    # Plot the window of interest
    ax.set_xlim(kwargs.get("xlim", (-0.25, 1.5)))
    ax.set_ylim(kwargs.get("ylim", (np.min(y), np.max(y))))

    lim = [-0.007, 0.007]
    # lim = kwargs.get('lims', lim)

    norm = colors.Normalize(vmin=lim[0], vmax=lim[1])
    levels = np.linspace(lim[0], lim[1], 51)

    _cmap = sns.color_palette("seismic", as_cmap=True)

    cs = ax.contourf(
        x,
        y,
        mag,
        levels=levels,
        vmin=lim[0],
        vmax=lim[1],
        norm=norm,
        cmap=_cmap,
        extend="both",
    )

    ax_cb = divider.new_horizontal(size="5%", pad=0.05)
    fig.add_axes(ax_cb)
    plt.colorbar(cs, cax=ax_cb)
    ax_cb.yaxis.tick_right()
    ax_cb.yaxis.set_tick_params(labelright=True)
    ax_cb.set_ylabel("$ \Omega_x $", rotation=0)

    # plt.setp(ax_cb.get_yticklabels()[::2], visible=False)
    ax.set_aspect(1)

    plt.savefig(fn_save, dpi=800, transparent=True)
    plt.close()


def vorticity_z(flow):
    dv_dx = np.gradient(flow.V, axis=0, edge_order=2)
    du_dy = np.gradient(flow.U, axis=1, edge_order=2)
    return dv_dx - du_dy


def vorticity_x(flow):
    dv_dz = np.gradient(flow.V, axis=2, edge_order=2)
    dw_dy = np.gradient(flow.W, axis=1, edge_order=2)
    return dv_dz - dw_dy


def vorticity_y(flow):
    du_dz = np.gradient(flow.U, axis=2, edge_order=2)
    dw_dx = np.gradient(flow.W, axis=0, edge_order=2)
    return du_dz - dw_dx


def main():
    cwd = os.getcwd()
    cs = np.array([384, 768, 896, 1024])
    cs = np.array([1024])
    for c in cs:
        os.system(f'mkdir -p {cwd}/figures/{str(c)}')
        phase_average = Decompositions(f"{cwd}/{str(c)}", "fluid", length_scale=c).phase_average(20)
        for idx, loop in enumerate(phase_average):
            snap = AssignProps(loop)
            plot(snap, f'{cwd}/figures/{str(c)}/flow{idx}.png', c)
        del phase_average


if __name__ == "__main__":
    main()
    # cwd = os.getcwd()
    # # flow = FlowBase(f"{cwd}/data", "fluid", length_scale=384)
    # phase_average = Decompositions(f"{cwd}/data", "fluid", length_scale=384).phase_average(30)
    # print(np.shape(phase_average))
    # print(np.shape(flow.snaps))
