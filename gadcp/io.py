#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module gadcp.io with in/out functions. Mostly provides wrapper functions to UHs `Multiread`."""

import numpy as np
import xarray as xr

import gvpy as gv
from pycurrents.adcp.rdiraw import Multiread


def read_raw_rdi(file, auxillary_only=False):
    """
    Read raw RDI ADCP data file and return as xarray Dataset.

    Parameters
    ----------
    file : str or Path
        Path to raw data file.
    auxillary_only : bool
        Set to True to ignore 2d fields. (default False)

    Returns
    -------
    rdi : xarray.Dataset
        Raw ADCP data in xarray Dataset format.
    """

    m = Multiread(file, "wh")

    if auxillary_only:
        radcp = m.read(varlist=["FixedLeader"])
    else:
        if "BottomTrack" in m.available_varnames:
            radcp = m.read(
                varlist=[
                    "Velocity",
                    "PercentGood",
                    "Intensity",
                    "Correlation",
                    "BottomTrack",
                ]
            )
        else:
            radcp = m.read(
                varlist=["Velocity", "PercentGood", "Intensity", "Correlation"]
            )
    # convert time
    adcptime = gv.io.yday0_to_datetime64(radcp.yearbase, radcp.dday)

    jj = np.squeeze(radcp.dep.shape)
    assert radcp.nbins == jj
    ii = np.squeeze(radcp.dday.shape)
    assert radcp.nprofs == ii
    varsii = [
        "num_pings",
        "dday",
        "ens_num",
        "temperature",
        "heading",
        "pitch",
        "roll",
        "XducerDepth",
    ]

    out = xr.Dataset(
        data_vars={"dummy": (["z", "time"], np.ones((jj, ii)) * np.nan)}
    )

    # get 1d variables
    for v in varsii:
        out[v] = (["time"], radcp[v])
    # add pressure
    out["pressure"] = (["time"], radcp.VL["Pressure"] / 1000)

    # get 2d variables
    if auxillary_only is False:
        for v in ["vel", "cor", "amp"]:
            out[v] = (["beam", "z", "time"], np.transpose(radcp[v]))
        if "pg" in radcp:
            out["pg"] = (["beam", "z", "time"], np.transpose(radcp["pg"]))

    out.coords["time"] = (["time"], adcptime)
    out.coords["z"] = (["z"], radcp.dep)

    # bottom tracking
    if "bt_vel" in radcp.keys():
        out["bt_vel"] = (["time", "beam"], radcp.bt_vel)
        out["bt_depth"] = (["time", "beam"], radcp.bt_depth)
    out.coords["beam"] = np.array([1, 2, 3, 4])

    # drop dummy
    out = out.drop(["dummy"])

    # set a few attributes
    out.attrs["sonar"] = radcp.sonar.sonar
    out.attrs["coordsystem"] = radcp.trans.coordsystem
    out.attrs["pingtype"] = radcp.pingtype
    out.attrs["cellsize"] = radcp.CellSize

    return out


def read_raw_rdi_uh(file, auxillary_only=False):
    """
    Wrapper for UH's pycurrents.adcp.rdiraw.Multiread

    Parameters
    ----------
    file : str or Path
        Path to raw data file.

    Returns
    -------
    radcp : dict (Bunch)
        UH data structure with raw RDI data
    """

    m = Multiread(file, "wh")

    if auxillary_only:
        radcp = m.read(varlist=["FixedLeader"])
    else:
        if "BottomTrack" in m.available_varnames:
            radcp = m.read(
                varlist=[
                    "Velocity",
                    "PercentGood",
                    "Intensity",
                    "Correlation",
                    "BottomTrack",
                ]
            )
        else:
            radcp = m.read(
                varlist=["Velocity", "PercentGood", "Intensity", "Correlation"]
            )

    # convert time
    adcptime = gv.io.yday0_to_datetime64(radcp.yearbase, radcp.dday)
    radcp.time = adcptime

    # pressure and temperature
    radcp.pressure = radcp.VL["Pressure"] / 1000.0
    radcp.temperature = radcp.VL["Temperature"] / 100.0

    return radcp


def _is_number(s):
    """
    Check if string can be converted to a float.

    Parameters
    ----------
    s : str
        string

    Returns
    -------
    out : bool
        True if string can be converted to float, else False.
    """
    try:
        float(s)
        return True
    except ValueError:
        return False
