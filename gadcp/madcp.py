#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module gadcp.madcp with functions for moored ADCPs"""

import os
import subprocess
import datetime
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import xarray as xr
from pathlib import Path

from pycurrents.adcp.rdiraw import extract_raw, Multiread
from pycurrents.system import Bunch
from pycurrents.adcp.mcm_avg import MCM, Pingavg

import gvpy as gv

from . import io


def proc(infile, lon, lat, end_pc, end_adcp):
    """Temporary function for processing ADCP raw data.

    Parameters
    ----------
    infile : str or Path
        ADCP data file.
    lon : float
        Mooring longitude
    lat : float
        Mooring latitude
    end_pc : tuple (int, int, int, int, int, int)
        Time at data download (UTC) as six element
        tuple (year, month, day, hour, minute, second).
    end_adcp : tuple (int, int, int, int, int, int)
        Time at data download (instrument).

    Returns
    -------
    TODO

    To Do
    -----
    - make time and depth grid parameters input parameters
    - allow for single ping processing/output, i.e. no time gridding
    - allow for external pressure time series input

    """
    print("reading configuration...")
    m = Multiread(infile, "wh")
    print("system configuration:", m.sysconfig)

    print("reading raw data...")
    raw = io.read_raw_rdi(infile)
    fig, ax = gv.plot.quickfig()
    raw.pressure.plot(ax=ax)
    ax.invert_yaxis()
    # find time away from surface
    ii = np.argwhere(raw.pressure.data > 50).squeeze()
    raw.pressure.isel(time=ii).plot(ax=ax, color="r")

    print("extract subsurface ping range")
    i0, i1 = ii[0], ii[-1]
    outfile = "adcp_cut.dat"
    inst = "wh"
    print("Extracting ping range %d to %d" % (i0, i1))
    data = extract_raw(infile, inst, i0, i1, outfile=outfile)

    # find start time in julian days in the raw time series
    t0 = raw.dday[0].data

    # read the cut time series
    mc = Multiread(outfile, inst)
    orientation = "up" if mc.sysconfig["up"] else "down"
    print("instrument orientation: ", orientation)

    # Parameters for time averaging and depth gridding
    outdir = "./"
    datadir = "./"
    fnamesdict = dict(adcp=[outfile])
    driftparamsdict = dict(adcp=dict(end_pc=end_pc,
                                     end_adcp=end_adcp,
                                     start_dday=t0))
    positionsdict = dict(adcp=(lon, lat))

    editparams = dict(max_e=0.2,          # absolute max e
                      max_e_deviation=2,  # max in terms of sigma
                      min_correlation=64) # 64 is RDI standard
    dgridparams = dict(dbot=1500, # int(self.p_median),
                       dtop=100,  # 50,
                       d_interval=16)
    tgridparams = dict(dt_hours = 1.0,  #  1.0/4,
                    t0 = 132,
                    #t1 = t1,
                    )

    print('time averaging and depth gridding')
    for key in fnamesdict.keys():
        mcm = MCM(fnamesdict[key], driftparamsdict[key], datadir=datadir)
        pa = Pingavg(mcm, lonlat=positionsdict[key],
                    dgridparams=dgridparams,
                    tgridparams=tgridparams,
                    editparams=editparams,
                    )
        pa.average_ensembles()
        npzname = key + '_hourly.npz'
        pa.save_npz(npzname, outdir=outdir)

    # load the generated file
    data = npz2nc(npzname)

    return m, data


def tmpload(file):
    mfd_prefix = "__mfd__"
    mfm_prefix = "__mfm__"
    fz = open(file, "rb")
    z = np.load(fz, mmap_mode=None, encoding='bytes', allow_pickle=True)
    out = Bunch()
    for k, v in z.items():
        if k.startswith(mfd_prefix):
            kma = k[len(mfd_prefix):]
            kmm = mfm_prefix + kma
            out[kma] = np.ma.array(v, mask=z[kmm], copy=False)
        elif k.startswith(mfm_prefix):
            continue
        else:
            if v.ndim == 0:
                v = v.item()
                if isinstance(v, dict):
                    v = Bunch(v)
            out[k] = v
    fz.close()

    return out


def npz2nc(npzfile):
    """Convert data structure from npz to xarray / netcdf file format.

    Parameters
    ----------
    npzfile : str
        ADCP data in npz file format.

    Returns
    -------
    out : xarray.Dataset
        ADCP data ready for saving to netcdf.
    """
    # load npz file
    dat = tmpload(npzfile)
    # identify variables
    k = dat.keys()
    varsint = []
    vars1d =  []
    vars2d = []
    for ki in k:
        try:
            tmp = dat[ki].shape
            if len(tmp)==1:
                vars1d.append(ki)
            elif len(tmp)==2:
                vars2d.append(ki)
        except:
            tmp = None
            varsint.append(ki)
        print(ki, tmp)

    # generate time vector
    base = datetime.datetime(dat.yearbase, 1, 1, 0, 0, 0)
    time = [base + datetime.timedelta(days=ti) for ti in dat.dday]
    adcptime = [np.datetime64(ti) for ti in time]
    # generate Dataset
    out = xr.Dataset({'pg': (['z', 'time'], dat.pg.T)},
                    coords={'time': (['time'], adcptime),
                            'z': (['z'], dat.dep)})
    for vari in vars2d:
        out[vari] = (['z', 'time'], dat[vari].T)
    for vari in vars1d:
        if vari not in ['dep', 'dday']:
            out[vari] = (['time'], dat[vari])

    return out
