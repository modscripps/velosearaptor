#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module gadcp.adcp with general adcp functions"""

import os
import matplotlib.pyplot as plt
import gvpy as gv


def plot_raw_adcp(adcp, figsize=(17, 20)):
    """
    Plot raw RDI adcp dataset.

    Parameters
    ----------
    adcp : xarray.Dataset
        Raw RDI ADCP data read using gvpy.io.read_raw_rdi()
    figsize : tuple
        Provide figure size (default (17, 20))
    """

    fig1 = plt.figure(figsize=figsize)
    gs0 = fig1.add_gridspec(nrows=16, ncols=90)


    def plot_time_series(v, rowid, yincrease=False):
        tsax.append(fig1.add_subplot(gs0[rowid, 0:70]))
        if isinstance(v, list):
            for vi in v:
                adcp[vi].plot(ax=tsax[-1], yincrease=yincrease, label=vi)
                gv.plot.ysym()
            tsax[-1].legend()
        else:
            adcp[v].plot(ax=tsax[-1], yincrease=yincrease)
        return tsax


    tsax = []
    tsax = plot_time_series("pressure", 0)
    tsax = plot_time_series(["pitch", "roll"], 13, True)
    tsax = plot_time_series("heading", 14, True)
    tsax = plot_time_series("temperature", 15, True)


    def plot_beam_quantity(startrow, v, vmin, vmax, cmap, clabel):
        velcax = fig1.add_subplot(gs0[startrow + 1 : startrow + 5, 70])
        velax = []
        for i in range(4):
            velax.append(fig1.add_subplot(gs0[i + startrow + 1, 0:70]))
        for axi, (g, b) in zip(velax, adcp[v].groupby("beam")):
            h = b.plot(
                ax=axi, add_colorbar=False, vmin=vmin, vmax=vmax, cmap=cmap, yincrease=False
            )
            axi.text(
                0.01,
                0.1,
                "beam {}".format(g),
                transform=axi.transAxes,
                bbox=dict(fc="w", ec="0.3", alpha=0.7, boxstyle="Round, pad=0.2"),
            )
        plt.colorbar(h, cax=velcax, label=clabel)
        for axi in velax:
            axi.set(xlabel="", title="")
        for axi in velax:
            axi.tick_params(labelbottom=False)
        gv.plot.concise_date(velax[-1], show_offset=False)

        return velax


    def plot_time_mean_beam_quantity(startrow, v):
        meanax = fig1.add_subplot(gs0[startrow + 1 : startrow + 4, 80:90])
        cors = [vb.mean(dim="time") for (g, vb) in adcp[v].groupby("beam")]
        for b, ai in enumerate(cors):
            ai.plot(
                ax=meanax,
                y="z",
                label="beam{}".format(b + 1),
                yincrease=False,
                marker="o",
                linestyle="",
                alpha=0.8,
            )
        meanax.set(title="")
        meanax.legend()


    velax = plot_beam_quantity(0, "vel", -1.5, 1.5, "RdBu_r", "velocity [m/s]")
    corax = plot_beam_quantity(0 + 4, "cor", 0, 150, "viridis", "correlation")
    ampax = plot_beam_quantity(0 + 8, "amp", 50, 200, "magma", "amplitude")

    plot_time_mean_beam_quantity(4, "cor")
    plot_time_mean_beam_quantity(8, "amp")

    gv.plot.concise_date(tsax[-1])
    tsax[-1].tick_params(labelbottom=True)
    for axi in tsax[:-1]:
        axi.tick_params(labelbottom=False)
    for axi in tsax:
        axi.set(xlim=velax[-1].get_xlim(), xlabel="")

    infoax = fig1.add_subplot(gs0[0:3, 80:90])
    infoax.text(
        -0.2,
        0.7,
        "instrument: {}\nping type: {}\ncoordinate system: {}\nbin size: {}".format(
            adcp.attrs["sonar"], adcp.attrs["pingtype"], adcp.attrs["coordsystem"], adcp.attrs["cellsize"]
        ),
        transform=infoax.transAxes,
    )
    infoax.axis("off")


def plot_raw_adcp_auxillary(adcp, figsize=(12, 5)):
    """
    Plot raw RDI adcp dataset (auxillary data only).

    Parameters
    ----------
    adcp : xarray.Dataset
        Raw RDI ADCP data read using gvpy.io.read_raw_rdi()
    figsize : tuple
        Provide figure size (default (17, 20))
    """
    fig1 = plt.figure(figsize=figsize)
    gs0 = fig1.add_gridspec(nrows=4, ncols=15)


    def plot_time_series(v, rowid, yincrease=False):
        tsax.append(fig1.add_subplot(gs0[rowid, 0:7]))
        if isinstance(v, list):
            for vi in v:
                adcp[vi].plot(ax=tsax[-1], yincrease=yincrease, label=vi)
                gv.plot.ysym()
            tsax[-1].legend()
        else:
            adcp[v].plot(ax=tsax[-1], yincrease=yincrease)
        return tsax

    def plot_time_mean_beam_quantity(startcol, v):
            meanax = fig1.add_subplot(gs0[1:, startcol:startcol+3])
            cors = [vb.mean(dim="time") for (g, vb) in adcp[v].groupby("beam")]
            for b, ai in enumerate(cors):
                ai.plot(
                    ax=meanax,
                    y="z",
                    label="beam{}".format(b + 1),
                    yincrease=False,
                    marker="o",
                    linestyle="",
                    alpha=0.8,
                )
            meanax.set(title="")
            meanax.legend()

    tsax = []
    tsax = plot_time_series("pressure", 0)
    tsax = plot_time_series(["pitch", "roll"], 1, True)
    tsax = plot_time_series("heading", 2, True)
    tsax = plot_time_series("temperature", 3, True)

    plot_time_mean_beam_quantity(8, "cor")
    plot_time_mean_beam_quantity(12, "amp")

    gv.plot.concise_date(tsax[-1])
    tsax[-1].tick_params(labelbottom=True)
    for axi in tsax[:-1]:
        axi.tick_params(labelbottom=False)
    for axi in tsax:
        axi.set(xlim=(adcp.time[0].data, adcp.time[-1].data), xlabel="")

    infoax = fig1.add_subplot(gs0[0, 8:])
    infoax.text(
        0.1,
        0.4,
        "instrument: {}\nping type: {}\ncoordinate system: {}".format(
            adcp.attrs["sonar"], adcp.attrs["pingtype"], adcp.attrs["coordsystem"]
        ),
        transform=infoax.transAxes,
    )
    infoax.axis("off")
