#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module gadcp.adcp with adcp functions"""

import os

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import xarray as xr

import gvpy as gv
from pycurrents.adcp.rdiraw import extract_raw


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


class LADCPYoYoSplit(object):
    """
    Split raw ADCP files from yoyo or towyo time series.

    Parameters
    ----------
    downlooker : str
        Path to raw downlooker data
    uplooker : str
        Path to raw uplooker data

    Returns
    -------
    var : dtype
        description
    """
    def __init__(self, downlooker, uplooker):
        self.dn = gv.io.read_raw_rdi_uh(downlooker, auxillary_only=True)
        self.up = gv.io.read_raw_rdi_uh(uplooker, auxillary_only=True)
        self.dnraw = downlooker
        self.upraw = uplooker

        assert self.up.sysconfig.up == True
        assert self.dn.sysconfig.up == False

    # quick plotting routine
    def ts_plot(self, adcp, var, ax=None):
        if ax is None:
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 2))
        ax.plot(adcp.time, adcp[var])
        if var in ['pressure', 'XducerDepth']:
            ax.set(ylabel=var)
            gv.plot.ydecrease(ax)
        gv.plot.concise_date(ax)
        return ax
    # plot both pressure time series
    def plot_pressure(self):
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(10, 4), sharex=True)
        self.ts_plot(self.up, 'pressure', ax1)
        self.ts_plot(self.dn, 'pressure', ax2);

    def save_temp_files(self):
        # bring to same start/end times if appropriate
        dnt0 = self.dn.time[0]
        upt0 = self.up.time[0]
        self.deltat0 = np.timedelta64(upt0-dnt0, 's')
        print('difference in start time:', self.deltat0)
        # cut time series at the beginning if time deltat less than ten minutes.
        self.cut_criterion = np.timedelta64(600, 's')
        if np.abs(self.deltat0) < self.cut_criterion:
            print('bringing to same start time')
            self.cut_beginning = True
            if dnt0 > upt0:
                dnii0 = 0
                upii0 = np.argmax(self.up.time > dnt0)
            else:
                upii0 = 0
                dnii0 = np.argmax(self.dn.time > upt0)
        else:
            dnii0 = 0
            upii0 = 0

        dnt1 = self.dn.time[-1]
        upt1 = self.up.time[-1]
        self.deltat1 = np.timedelta64(dnt1-upt1, 's')
        print('difference in end time:', self.deltat1)
        # cut time series at the beginning if time deltat less than ten minutes.
        if np.abs(self.deltat1) < self.cut_criterion:
            print('bringing to same end time')
            if dnt1 > upt1:
                upii1 = np.argmax(self.up.time)
                dnii1 = np.amin(np.where(self.dn.time > upt1))
            else:
                dnii1 = np.argmax(self.dn.time)
                upii1 = np.amin(np.where(self.up.time > dnt1))
        else:
            upii1 = np.argmax(self.up.time)
            dnii1 = np.argmax(self.dn.time)

        # cut start/end if needed and save to current working directory.
        # we will delete this temporary file later on.
        print('\nsaving temporary files to current working directory\n')
        extract_raw(self.dnraw, 'wh', dnii0, dnii1, outfile='tmpdn.rdi', verbose=False);
        extract_raw(self.upraw, 'wh', upii0, upii1, outfile='tmpup.rdi', verbose=False);
        print('reading temporary files\n')
        self.dn = gv.io.read_raw_rdi_uh('tmpdn.rdi', auxillary_only=True)
        self.up = gv.io.read_raw_rdi_uh('tmpup.rdi', auxillary_only=True)


    # find split points
    def _lpfilt(self, data, pts):
        # design the Buterworth filter
        N  = 3    # Filter order
        Wn = 1/(pts*2) # Cutoff frequency
        B, A = sp.signal.butter(N, Wn, output='ba')
        # apply the filter
        lpdata = sp.signal.filtfilt(B,A,data)
        return lpdata

    def _split_point_function(self, dnup, dx, mindz):
        if dnup == 'dn':
            pressure = self.dn.pressure
        elif dnup == 'up':
            pressure = self.up.pressure
        # low pass filter pressure time series
        lptmp = self._lpfilt(pressure, 60)
        # find relative maxima and minima within +/-dx points
        tmpmax = sp.signal.argrelmax(lptmp, order=dx)[0]
        tmpmin = sp.signal.argrelmin(lptmp, order=dx)[0]

        allmaxmin = np.concatenate((tmpmax, tmpmin))
        allmaxminind = np.concatenate((np.ones_like(tmpmax), np.zeros_like(tmpmin)))

        allmaxminsortedind = np.argsort(allmaxmin)
        allmaxminsorted = allmaxmin[allmaxminsortedind]
        allmaxminsorteddep = pressure[allmaxminsorted]
        allmaxminsorteddepdiff = np.diff(allmaxminsorteddep)
        allmaxminsorteddepdiffabs = np.abs(allmaxminsorteddepdiff)

        pickmaxminsortedind1 = allmaxminsorteddepdiffabs > mindz
        pickmaxminsortedind = np.append(pickmaxminsortedind1, True)

        pickmaxminsorted = allmaxminsorted[pickmaxminsortedind]
        pickmaxminind = allmaxminind[allmaxminsortedind[pickmaxminsortedind]]

        shallow_indices = pickmaxminsorted[pickmaxminind == 0]
        deep_indices = pickmaxminsorted[pickmaxminind == 1]

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15, 3))
        ax.plot(lptmp, color='0.5')
        ax.plot(deep_indices, pressure[deep_indices], 'ro')
        ax.plot(shallow_indices, pressure[shallow_indices], 'bo')
        ax.set(title=dnup)
        gv.plot.ydecrease(ax=ax)

        return shallow_indices, deep_indices

    def find_split_points(self, dx=500, mindz=300):
        self.ind_shallow_dn, self.ind_deep_dn = self._split_point_function('dn', dx=dx, mindz=mindz)
        self.ind_shallow_up, self.ind_deep_up = self._split_point_function('up', dx=dx, mindz=mindz)

    def plot_indices(self):
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 4))
        ax.plot(self.ind_shallow_up - self.ind_shallow_dn, 'o')

    def split_files(self, experiment, cruise, cast, outdir):
        for adcp, indices, dnup in zip([self.dn, self.up], [self.ind_shallow_dn, self.ind_shallow_up], ['dn', 'up']):
            for ni, (ind0, ind1) in enumerate(zip(indices[:-1], indices[1:])):
                t0 = np.datetime64(adcp.time[ind0])
                t1 = np.datetime64(adcp.time[ind1])
                i0 = int(ind0)
                i1 = int(ind1)
                tbeg = np.datetime64(t0, 's').astype('str').replace('T', '_').replace(':', '')
                tend = np.datetime64(t1, 's').astype('str').replace('T', '_').replace(':', '')
                filename = '{}_{}_cast{:03g}_yo{:03g}_{}_{}_{}.rdi'.format(experiment, cruise, cast, ni, dnup, tbeg, tend)
                outfile = os.path.join(outdir, filename)
                rawfile = 'tmp{}.rdi'.format(dnup)
                extract_raw(rawfile, 'wh', i0, i1, outfile=outfile, verbose=False);
