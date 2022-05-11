#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module gadcp.io with in/out functions. Mostly provides wrapper functions to UHs `Multiread`."""

import datetime
import numpy as np
import xarray as xr

from pycurrents.adcp.rdiraw import Multiread, extract_raw


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
    adcptime = yday0_to_datetime64(radcp.yearbase, radcp.dday)

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

    out = xr.Dataset(data_vars={"dummy": (["z", "time"], np.ones((jj, ii)) * np.nan)})

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
    adcptime = yday0_to_datetime64(radcp.yearbase, radcp.dday)
    radcp.time = adcptime

    # pressure and temperature
    radcp.pressure = radcp.VL["Pressure"] / 1000.0
    radcp.temperature = radcp.VL["Temperature"] / 100.0

    return radcp


def extract_raw_rdi(file, i0, i1, outfile, inst="wh"):
    """Extract ping range from raw RDI file and save to new raw file.

    Parameters
    ----------
    file : str or Path
        Input file.
    i0 : int
        Start ping.
    i1 : int
        End ping.
    outfile : str or Path
        Path and name of output file.
    inst : str
        One of ('wh','os','bb','ec'). Defaults to 'wh'.
    """
    data = extract_raw(file, inst, i0, i1, outfile=outfile)


def yday0_to_datetime64(baseyear, yday):
    """
    Convert year day (starting at yday 0) to numpy's datetime64 format.

    Parameters
    ----------
    baseyear : int
        Base year
    yday : float
        Year day

    Returns
    -------
    time : np.datetime64
        Time in numpy datetime64 format
    """
    base = datetime.datetime(baseyear, 1, 1, 0, 0, 0)
    time = [base + datetime.timedelta(days=ti) for ti in yday]
    # convert to numpy datetime64
    time64 = np.array([np.datetime64(ti, "ms") for ti in time])
    return time64


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


def cf_conventions():
    """Return dictionary with CF standard names and units.

    See https://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html
    and http://cfconventions.org/cf-conventions/cf-conventions.html
    for further information.

    Returns
    -------
    CF : dict
        CF standard names in units.
    """

    CF = dict(
        depth=dict(
            long_name="depth",
            standard_name="depth",
            units="m",
            positive="down",
        ),
        u=dict(
            long_name="u",
            standard_name="eastward_sea_water_velocity",
            units="m s-1",
            ancillary_variables="npings",
        ),
        v=dict(
            long_name="v",
            standard_name="northward_sea_water_velocity",
            units="m s-1",
            ancillary_variables="npings",
        ),
        w=dict(
            long_name="w",
            standard_name="upward_sea_water_velocity",
            units="m s-1",
            ancillary_variables="npings",
        ),
        e=dict(
            long_name="error velocity",
            standard_name="indicative_error_from_multibeam_acoustic_doppler_velocity_profiler_in_sea_water",
            units="m s-1",
            ancillary_variables="npings",
        ),
        pressure=dict(
            long_name="pressure",
            standard_name="sea_water_pressure",
            units="dbar",
            ancillary_variables="npings",
        ),
        temperature=dict(
            long_name="temperature",
            standard_name="sea_water_temperature",
            units="degrees_C",
            ancillary_variables="npings",
        ),
        xducer_depth=dict(
            long_name="transducer depth",
            standard_name="depth",
            units="m",
            positive="down",
        ),
        npings=dict(
            long_name="number of pings",
            standard_name="number_of_observations",
        ),
        pg=dict(
            long_name="percent good",
            standard_name="proportion_of_acceptable_signal_returns_from_acoustic_instrument_in_sea_water",
            units="percent",
            ancillary_variables="npings",
        ),
        amp=dict(
            long_name="amplitude",
            standard_name="signal_intensity_from_multibeam_acoustic_doppler_velocity_sensor_in_sea_water",
            ancillary_variables="npings",
        ),
    )
    return CF
