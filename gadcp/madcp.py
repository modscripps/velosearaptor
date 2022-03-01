#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module gadcp.madcp with functions for moored ADCPs."""

import os
import subprocess
import datetime
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import xarray as xr
from pathlib import Path
import logging

from pycurrents.adcp.rdiraw import extract_raw, Multiread
from pycurrents.system import Bunch

from gadcp.mcm_avg import MCM, Pingavg

import gvpy as gv

logger = logging.getLogger(__name__)


def proc(
    infile,
    lon,
    lat,
    editparams,
    tgridparams,
    dgridparams,
    end_pc=None,
    end_adcp=None,
    meta_data=None,
    n_ensembles=None,
    ibad=None,
    pressure_scale_factor=1,
    logdir="log",
    verbose=False,
    plot_pressure=False,
):
    """Process moored ADCP raw data.

    Parameters
    ----------
    infile : str or Path
        ADCP data file.
    lon : float
        Mooring longitude
    lat : float
        Mooring latitude
    editparams : dict
        Raw data editing parameters.
    tgridparams : dict
        Time gridding parameters.
    dgridparams : dict
        Depth gridding parameters.
    end_pc : tuple (int, int, int, int, int, int), optional
        Time at data download (UTC) as six element
        tuple (year, month, day, hour, minute, second). Can be omitted if no
        time drift was recorded (not recommended).
    end_adcp : tuple (int, int, int, int, int, int), optional
        Time at data download (instrument).
    meta_data : dict, optional
        Provide additional metadata. All entries will be written into the attrs
        of the output `xr.Dataset`. If entries mooring and sn are provided,
        they are used to name the logfile.
    n_ensembles : None or int, optional
        Extract only this number of ensembles (after the instrument is in the
        water). Mostly for testing purposes. Defaults to None.
    ibad : int, optional
        The index of a beam to be excluded from the beam_to_xyz
        calculation. This is the zero-based index. Defaults to None.
    pressure_scale_factor : float, optional
        Scale factor for pressure time series. Defaults to 1 (no scaling).
    logdir : str, optional
        Directory for logfile. Defaults to 'log' and will create directory if
        it doesn't exist.
    verbose : bool, optional
        Print all log information to screen if True. Defaults to False (only warnings).
    plot_pressure : bool, optional
        Plot pressure time series and indicate subsurface pings if True. Defaults to False.

    Returns
    -------
    TODO

    To Do
    -----
    - allow for single ping processing/output, i.e. no time gridding
    - allow for external pressure time series input

    """
    # Parse metadata for naming the log.
    if meta_data is not None:
        if "sn" in meta_data and "mooring" in meta_data:
            sn = meta_data["sn"]
            mooring = meta_data["mooring"]
            has_mooring_id = True
        else:
            sn = None
            mooring = None
            has_mooring_id = False
    else:
        has_mooring_id = False

    # Set up a logger that will collect log info from this module and the
    # other pycurrent methods as well.
    logdir = Path(logdir)
    logdir.mkdir(exist_ok=True)
    if has_mooring_id:
        filename = f"{mooring}_{sn}.log"
    else:
        filename = "adcp_proc.log"
    # Delete any existing handlers. This may be bad style, but I kept adding handlers when developing this.
    logger.handlers = []
    logging.basicConfig(
        filename=logdir.joinpath(filename),
        filemode="w",
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
        level=logging.INFO,
        force=True,
    )
    ConsoleOutputHandler = logging.StreamHandler()
    ConsoleOutputHandler.setLevel(logging.WARNING)
    if verbose:
        ConsoleOutputHandler.setLevel(logging.INFO)
    logger.addHandler(ConsoleOutputHandler)
    datestr = gv.time.now_datestr()
    if has_mooring_id:
        logger.info(f"Processing {mooring} SN{sn} on {datestr}")
    else:
        logger.warning(
            "No meta data provided, logging to generic filename 'adcp_proc.log'"
        )
        logger.info(f"Processing on {datestr}")

    # Read configuration.
    m = Multiread(infile, "wh")
    # logger.info(f"system configuration: {m.sysconfig}")

    # Read raw data.
    raw = m.read(varlist=["VariableLeader"])
    raw.pressure = raw.VL["Pressure"] / 1000.0 * pressure_scale_factor
    raw.temperature = raw.VL["Temperature"] / 100.0

    # Find start time in julian days in the raw time series.
    t0 = raw.dday[0]

    # Find time away from surface.
    ii = np.argwhere(raw.pressure > 50).squeeze()

    if plot_pressure:
        fig, ax = gv.plot.quickfig(fgs=(6, 2.5))
        ax.plot(raw.dday, raw.pressure, label="all")
        ax.invert_yaxis()
        ax.set(xlabel="julian day", ylabel="pressure [dbar]")
        ax.plot(raw.dday[ii], raw.pressure[ii], color="r", label="subsurface")

    if n_ensembles is not None:
        if isinstance(infile, list):
            if len(infile) == 1:
                infile = infile[0]
            else:
                raise Exception(
                    "Can't have more than one file if using n_ensembles."
                )

        # Wondering if we really need to extract the subsurface range here? MCM
        # also looks for data in the water and discards the rest unless specified
        # otherwise. Currently, if we want to process less than the full time
        # series, we depend on the subset that is extracted here. This should be
        # easy to change in the future if desired. Writing a subset of the data to
        # disk does not seem to slow down the process very much.
        print("extracting subsurface ping range...")
        i0, i1 = ii[0], ii[-1]
        if n_ensembles is not None:
            i1 = i0 + n_ensembles
        cut_file = "adcp_cut.dat"
        inst = "wh"
        # print("extracted ping range %d to %d" % (i0, i1))
        _ = extract_raw(infile, inst, i0, i1, outfile=cut_file)
        fnames = [cut_file]
    else:
        fnames = infile

    # Parameters for time averaging and depth gridding
    outdir = "./"
    datadir = "./"
    driftparams = dict(end_pc=end_pc, end_adcp=end_adcp, start_dday=t0)
    positions = (lon, lat)

    # editparams = dict(
    #     max_e=0.2,  # absolute max e
    #     max_e_deviation=2,  # max in terms of sigma
    #     min_correlation=64,
    # )  # 64 is RDI standard
    # dgridparams = dict(
    #     dbot=1500, dtop=100, d_interval=16
    # )  # int(self.p_median),  # 50,
    # tgridparams = dict(
    #     # dt_hours=1.0,  #  1.0/4,
    #     # t0=132,
    #     # t1 = t1,
    #     burst_average=True,
    # )
    logger.info("processing settings")
    logger.info("-------------------")
    _log_params(editparams)
    _log_params(tgridparams)
    _log_params(dgridparams)
    logger.info("-------------------")

    # Time averaging and depth gridding.
    mcm = MCM(
        fnames,
        driftparams,
        datadir=datadir,
        lat=lat,
        ibad=ibad,
        pressure_scale_factor=pressure_scale_factor,
    )
    pa = Pingavg(
        mcm,
        lonlat=positions,
        dgridparams=dgridparams,
        tgridparams=tgridparams,
        editparams=editparams,
    )
    pa.average_ensembles()
    npzname = "adcp_proc.npz"
    pa.save_npz(npzname, outdir=outdir)

    if n_ensembles is not None and plot_pressure:
        ax.plot(
            mcm.tsdat.dday,
            mcm.tsdat.pressure,
            color="orange",
            label="n_ensembles",
        )
    if plot_pressure:
        ax.legend()

    # load the generated file as netcdf / xarray.Dataset
    data = npz2nc(npzname)

    # Add some more info.
    data.attrs["orientation"] = mcm.orientation
    data.attrs["magdec"] = pa.ave.magdec
    for att in ["max_e", "max_e_deviation", "min_correlation"]:
        data.attrs[att] = pa.ave.editparams[att]

    # Add meta data if provided.
    if meta_data is not None:
        for k, v in meta_data.items():
            data.attrs[k] = v
    data.attrs["proc time"] = np.datetime64("now").astype("str")

    # remove npz file
    os.remove(npzname)

    return m, mcm, pa, data


def npzload(file):
    mfd_prefix = "__mfd__"
    mfm_prefix = "__mfm__"
    fz = open(file, "rb")
    z = np.load(fz, mmap_mode=None, encoding="bytes", allow_pickle=True)
    out = Bunch()
    for k, v in z.items():
        if k.startswith(mfd_prefix):
            kma = k[len(mfd_prefix) :]
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
    dat = npzload(npzfile)
    # identify variables
    k = dat.keys()
    varsint = []
    vars1d = []
    vars2d = []
    for ki in k:
        try:
            tmp = dat[ki].shape
            if len(tmp) == 1:
                vars1d.append(ki)
            elif len(tmp) == 2:
                vars2d.append(ki)
        except:
            tmp = None
            varsint.append(ki)
        # print(ki, tmp)

    # generate time vector
    base = datetime.datetime(dat.yearbase, 1, 1, 0, 0, 0)
    time = [base + datetime.timedelta(days=ti) for ti in dat.dday]
    adcptime = [np.datetime64(ti) for ti in time]
    # generate Dataset
    out = xr.Dataset(
        {"pg": (["z", "time"], dat.pg.T)},
        coords={"time": (["time"], adcptime), "z": (["z"], dat.dep)},
    )
    for vari in vars2d:
        out[vari] = (["z", "time"], dat[vari].T)
    for vari in vars1d:
        if vari not in ["dep", "dday"]:
            out[vari] = (["time"], dat[vari])

    # add variable names and units for plotting
    out = add_meta_data(out)

    return out


def add_meta_data(ds):
    ds.u.attrs = dict(long_name="u", units="m/s")
    ds.v.attrs = dict(long_name="v", units="m/s")
    ds.w.attrs = dict(long_name="w", units="m/s")
    ds.e.attrs = dict(long_name="error velocity", units="m/s")
    ds.z.attrs = dict(long_name="depth", units="m")
    ds.temperature.attrs = dict(long_name="temperature", units="°C")
    ds.pressure.attrs = dict(long_name="pressure", units="dbar")
    return ds


def _log_params(pd):
    """Print ADCP processing parameter dict to logs.

    Parameters
    ----------
    pd : dict
        Parameter dict.
    """

    for k, v in pd.items():
        if k == "maskbins":
            logstr = (k, ":", np.flatnonzero(v))
            logger.info(logstr)
        else:
            logstr = (k, ":", v)
            logger.info(logstr)
