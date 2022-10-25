#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module velosearaptor.madcp with functions for moored ADCPs.

### General Use
Processing a raw ADCP file from a moored deployment is a two-step process.
- Instantiate a processing object with `ProcessADCP` or `ProcessADCPyml`. The
  former expects the path to the raw data and a number of dictionaries with
  processing parameters as input. The latter reads a .yml file containing the
  processing parameters. See the respective docstrings for more information.
- Process raw pings and possibly run a ping-averaging method on the data.
  Options here are `ProcessADCP.process_pings`,
  `ProcessADCP.average_ensembles`, and `ProcessADCP.burst_average_ensembles`.

### Notes
Some general notes for this module.

#### Depth Gridding
The depth vector for the ADCP raw data (in instrument coordinates) is
calculated in :meth:`pycurrents.rdiraw.FileBBWHOS` as <br>
`dep = np.arange(NCells) * CellSize + Bin1Dist`
The depth vector thus points to the center of each bin.

From the RDI manual *WorkHorse Monitor, Sentinel, Mariner, Quartermaster, and
Long Ranger ADCPs Commands and Output Data Format*:

> This [Bin1Dist] field contains the distance to the middle of the first depth
> cell (bin). This distance is a function of depth cell length (WS), the
> profiling mode (WM), the blank after transmit distance (WF), and speed of
> sound.

"""

import datetime
import logging
import os
import pathlib
from pathlib import Path
from shutil import which
from subprocess import PIPE, Popen  # for magdec
from warnings import warn

import gsw
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import xarray as xr
from pycurrents.adcp.rdiraw import Multiread
from pycurrents.adcp.transform import Transform, rdi_xyz_enu
from pycurrents.codas import to_date, to_day
from pycurrents.data import seawater
from pycurrents.num import interp1
from pycurrents.num.nptools import rangeslice
from pycurrents.system import Bunch
from tqdm import tqdm

from . import io

# Standard logging
logger = logging.getLogger(__name__)


class ProcessADCP:
    """Moored ADCP Processing.

    An instance of ProcessADCP is initialized by providing raw data and
    processing parameters. Once initialized, the instance method
    :meth:`average_ensembles` can be used to average over just a few or
    all pings.

    Magnetic declination is automatically calculated via a call to the command
    line tool
    [magdec](https://currents.soest.hawaii.edu/hgstage/geomag/file/tip) that
    must be installed. Once calculated, the magnetic declination is stored in
    `magdec` and automatically applied to the data.

    Prior to time averaging, data are edited based on correlation and error
    velocity thresholds. The variable `pg` is calculated based on the number of
    excluded pings, i.e. it is the numbur of good pings divided by the number
    of total pings within one depth bin and one time bin. It may be used to
    further filter the data.

    A pdf document on ADCP data collection and processing principles can be
    downloaded from RDI
    [here](https://www.comm-tec.com/Docs/Manuali/RDI/BBPRIME.pdf).

    Time-averaged data are grouped together in an `xarray.Dataset` in the
    instance attribute `ds` for easy access and convenient output to netcdf
    format.

    Parameters
    ----------
    raw_data : str or list or Path
        Location(s) of raw data.
    meta_data : dict
        Dictionary with meta data. At a minimum entries for `lon` and `lat` are
        needed. If `mooring` and `sn` provide mooring name and serial number
        then these will be used to name the log file produced during
        processing.
    driftparams : dict, optional
        Time drift parameters. See notes below.
    tgridparams : dict, optional
        Time gridding parameters. See notes below.
    dgridparams : dict, optional
        Depth gridding parameters. See notes below.
    editparams : dict, optional
        Editing parameters. See notes below.
    ibad : int, optional
        Mark beam with bad data (zero based). Defaults to None.
    logdir : str, optional
        Log file directory. Defaults to `log/`.
    verbose : bool, optional
        Output more processing info to screen.
    magdec : float, optional
        Magnetic declination in degrees.
    pressure : xr.DataArray, optional
        Pressure time series in dbar as xr.DataAarray with coordinate `time`.

    Attributes
    ----------
    files : list
        List pointing to raw data file(s).
    dday_start : float
        Start time of ADCP time series, determined either through `t0` in
        `tgridparams` or the start of the time series once at depth.
    dday_end : float
        End time of ADCP time series, determined either through `t1` in
        `tgridparams` or the end of the time series at depth.
    start_ddays : list
        Start times of averaging intervals
    dday_mid : float
        Time stamp for ping average as determined in :meth:`make_start_ddays`.
    dt : float
        Time that is inclusive of one average. For regular averaging, this is
        `dt_hours` in `tgridparams` converted to Julian days. For
        burst-averagint, this is a time interval that is inclusive of all pings
        in one bursts (but goes slightly beyond that into the time between
        bursts).
    time_drift_rate : float
        Clock drift calculated from `driftparams`.
    orientation : str
        Instrument orientation `up` or `down`.
    magdec : float
        Magnetic declination.
    m : pycurrents.adcp.rdiraw.Multiread
        Multiread instance.
    tsdat : Bunch
        Auxiliary data.
    raw : xarray.Dataset
        Raw ADCP data.
    ds : xarray.Dataset
        Time-averaged dataset added by running :meth:`average_ensembles`.

    Notes
    -----
    Various parameters are passed to the instance through dictionaries. Their
    specifics are described below. Gridding and editing parameters can be
    updated after creating the ProcessADCP instance via `parse_dgridparams`,
    `parse_tgridparams`, and  `parse_editparams`.

    **Time drift parameters**
    Provide clock drift parameters via `driftparams`. Accepted entries are
    - `end_adcp` : Time of ADCP at data download
    - `end_pc` : UTC time at data download

    The difference between `end_adcp` and `end_pc` is used to linearly correct
    for instrument clock drift.

    **Time gridding parameters**
    Provide time gridding parameters via `tgridparams`. Accepted entries are
    - `dt_hours` : Time grid interval. Defaults to 0.5h.
    - `t0` : Start time for gridding. Determined from data if not provided.
    - `t1` : End time for gridding. Determined from data if not provided.
    - burst_average : bool
        Set ensemble averaging to act on burst sampling scheme. Defaults to False.

    **Depth gridding parameters**
    Provide depth gridding parameters via `dgridparams`. Accepted entries are
    - `dtop` : Shallow depth in m.
    - `dbot` : Deep depth in m.
    - `dinterval` : Vertical grid size in m. Defaults to 5m.

    Values for `dbot` and `dtop` are generated if not provided.

    **Editing parameters**
    Provide editing parameters via `editparams`.
    - `max_e`=0.2,  # absolute max e
    - `max_e_deviation`=2,  # max in terms of sigma
    - `min_correlation`=64,  # 64 is RDI default
    - `maskbins` : Array with booleans indexing into the ADCP bins. Use the
      convenience method `generate_binmask`.
    - `pg_limit` : float or int or None.
            Percent good limit applied prior to interpolating to the universal
            depth grid in `burst_average_ensembles`.

    """

    # Default editing parameters.
    _editparams = dict(
        max_e=0.2,  # absolute max e
        max_e_deviation=2,  # max in terms of sigma
        min_correlation=64,  # 64 is RDI default
        maskbins=None,  # do not mask any bins
        pg_limit=50,  # percent good limit applied in `burst_average_ensembles`
    )

    def __init__(
        self,
        raw_data,
        meta_data,
        driftparams=None,
        tgridparams=None,
        dgridparams=None,
        editparams=None,
        ibad=None,
        logdir="log",
        verbose=False,
        plot=False,
        pressure_scale_factor=1,
        magdec=None,
        pressure=None,
    ):
        self.meta_data = Bunch(meta_data.copy())
        self.ibad = ibad
        self.logdir = logdir
        self.verbose = verbose

        self._magdec_provided = magdec
        self._magdec = magdec

        self._pressure_provided = pressure
        self._pressure_scale_factor = pressure_scale_factor

        self._raw = None
        self._default_dgridparams = None

        self.parse_file_locations(raw_data)
        self._initiate_data_reader()
        self._read_auxiliary_data()
        self._parse_meta_data()

        self._set_up_logger()

        self.parse_driftparams(driftparams)
        self._parse_sysconfig()
        self.parse_dgridparams(dgridparams)
        self.parse_tgridparams(tgridparams)
        self.parse_editparams(editparams)

        self.make_start_ddays()

        if plot:
            self.plot_pressure()

    def parse_file_locations(self, raw_data, min_file_size=1e4):
        """Parse input for raw data files.

        Input can either be a single file name as a str, a single file as a
        Path instance, a list of either of these, or a Path instance pointing to a
        directory with raw ADCP files. In the latter case, files that are
        smaller than a threshold will not be included in the processing.

        Outputs to attribute `files`. The output can be fed to Multiread instances.

        Parameters
        ----------
        raw_data : str or list or Path
            Location(s) of raw data.
        min_file_size : int
            Minimum size for file to be included. Defaults to 1e4 which
            corresponds to about 10kB and is a good value for excluding small
            files without any actual data.

        """

        def list_dir(dir, min_file_size):
            # List all raw files.
            all_raw_files = list(sorted(dir.glob("*.00*")))
            # only files larger than about 10kB
            files = [
                file.as_posix()
                for file in all_raw_files
                if file.stat().st_size > min_file_size
            ]
            return files

        input_type = type(raw_data)
        if input_type is list:
            if type(raw_data[0]) is str:
                self.files = raw_data
            elif type(raw_data[0]) is pathlib.PosixPath:
                self.files = [file.as_posix() for file in raw_data]
        elif input_type is str:
            if Path(raw_data).is_dir():
                self.files = list_dir(Path(raw_data), min_file_size)
            else:
                self.files = [raw_data]
        elif input_type is pathlib.PosixPath:
            if raw_data.is_dir():
                self.files = list_dir(raw_data, min_file_size)
            else:
                self.files = [raw_data.as_posix()]

    def _initiate_data_reader(self):
        """Initiate a Multiread data reader.

        Adds attribute `m`.

        """
        self.m = Multiread(self.files, sonar="wh", ibad=self.ibad)
        # Make some more meta data realily available by reading a single ping
        # from the raw data.
        ping = self.m.read(start=0, stop=1)
        self.meta_data.Bin1Dist = ping.FL.Bin1Dist / 100.0
        self.meta_data.NCells = ping.FL.NCells
        self.meta_data.CellSize = ping.FL.CellSize / 100.0

    def _generate_external_pressure_interpolator(self, dat):
        """Interpolate external pressure to pycurrents time vector.

        Parameters
        ----------
        dat : pycurrents.adcp.rdiraw.Bunch
            Data structure.

        Returns
        -------
        Interpolator

        """

        # We already have a function to convert the pycurrents time to
        # datetime64. Interpolate in this time domain.
        time = io.yday0_to_datetime64(dat.yearbase, dat.dday)
        p_interpolated = self._pressure_provided.interp(time=time).data

        # Beginning and end may have NaN's if the ADCP was started before
        # the external pressure sensor. Make sure this happens only when
        # outside the water and then replace with atmospheric pressure.
        if np.any(np.isnan(p_interpolated)):
            nan_mask = np.isnan(p_interpolated)
            i_nan = np.flatnonzero(nan_mask)

            first_few_good_median = np.median(p_interpolated[~nan_mask][:20])
            last_few_good_median = np.median(p_interpolated[~nan_mask][-20:])

        if np.isnan(p_interpolated[0]):
            if first_few_good_median < 1:
                i_divide = np.flatnonzero(np.diff(i_nan) - 1) + 1
                if i_divide.size == 0:
                    p_interpolated[nan_mask] = first_few_good_median
                else:
                    p_interpolated[0:i_divide] = first_few_good_median

        if np.isnan(p_interpolated[-1]):
            if last_few_good_median < 1:
                i_divide = np.flatnonzero(np.diff(i_nan) - 1) + 1
                if i_divide.size == 0:
                    p_interpolated[nan_mask] = last_few_good_median
                else:
                    p_interpolated[i_divide:] = last_few_good_median

        # Now generate an interpolation function that will take dday as input
        # for later per-ensemble interpolation.
        self._external_pressure_interpolator = sp.interpolate.interp1d(
            dat.dday,
            p_interpolated,
            bounds_error=False,
            fill_value="extrapolate",
        )

    def _external_pressure_to_dat(self, dat):
        """Interpolate external pressure to pycurrents time vector.

        Parameters
        ----------
        dat : pycurrents.adcp.rdiraw.Bunch
            Data structure.

        Returns
        -------
        array-like
            Pressure

        """
        return self._external_pressure_interpolator(dat.dday)

    def _scale_pycurrents_pressure(self, dat):
        # Initial pressure units: 10 Pa (about 1 mm or 0.001 decibar).
        # Converting to decibars.
        return dat.VL["Pressure"] / 1000.0 * self._pressure_scale_factor

    def _read_auxiliary_data(self):
        """Read auxiliary data.

        Adds attribute `tsdat`.

        """
        tsdat = self.m.read(varlist=["VariableLeader"])
        tsdat.temperature = tsdat.VL["Temperature"] / 100.0
        # Replace pressure if provided from external sensor
        if self._pressure_provided is not None:
            self._generate_external_pressure_interpolator(tsdat)
            tsdat.pressure = self._external_pressure_to_dat(tsdat)
        else:
            tsdat.pressure = self._scale_pycurrents_pressure(tsdat)
        self.tsdat = tsdat

    def _parse_meta_data(self):
        """Parse meta data.

        - Add essential meta data to attributes. Will throw a KeyError if no
          lon/lat provided.

        - Read serial number from raw data and complain if it does not match
          the meta data SN.
        """
        essential_meta_data = ["lon", "lat"]

        [
            self._safely_add_attribute_from_params(k, self.meta_data)
            for k in essential_meta_data
        ]

        # Check SN
        sn_internal = int.from_bytes(self.tsdat.FL.Inst_SN, "little")

        if "sn" not in self.meta_data:
            self.meta_data.sn = sn_internal

        # Check internal SN matches user set one
        if sn_internal != self.meta_data.sn:
            warn(
                f"Serial number in file, {sn_internal}, is different from that set by user, {self.meta_data.sn}. Keeping user value."
            )

    @property
    def default_dgridparams(self):
        """Determine default depth gridding parameters.

        The grid is centered on the  median depth of the ADCP (plus distance to
        the center of the first bin) to avoid unnecessary binning into
        neighboring depth cells. The default size of the depth bins mimicks the
        size of ADCP bins.
        """
        if self._default_dgridparams is None:
            # Only use pressure at depth, not on deck
            ii = np.flatnonzero(self.tsdat.pressure > 15)
            p = self.tsdat.pressure[ii]
            # Determine limits of the pressure distribution but leave out the
            # top 5 and bottom 2 percent of data points. This way we are hoping
            # to avoid any outliers and a possible pressure record of ascent
            # and/or descent.
            p_top, p_bot = np.round(
                seawater.depth2(np.percentile(p, [5, 98]), self.lat)
            )
            self.p_median = np.median(p)
            pdep_median = np.round(seawater.depth2(self.p_median, self.lat))
            n = self.meta_data.NCells + 2
            d_interval = self.meta_data.CellSize
            distance_to_first_bin = np.round(self.meta_data.Bin1Dist)
            distance_to_last_bin = distance_to_first_bin + n * d_interval

            if self.sysconfig["up"]:
                # Set minimum grid depth level. Anything shallower than 10m
                # will be garbage anyways so let's throw this out.
                dtop = p_top - distance_to_last_bin
                if dtop < 10:
                    dtop = 10
                dbot = pdep_median - distance_to_first_bin
                # Successively add bins until we reach maximum pressure. This
                # will take care of mooring knockdowns.
                while dbot < p_bot:
                    dbot += d_interval
            else:
                dtop = pdep_median + distance_to_first_bin
                while dtop > p_top:
                    dtop -= d_interval
                dbot = p_bot + distance_to_last_bin

            self._default_dgridparams = dict(
                dtop=dtop,
                dbot=dbot,
                d_interval=d_interval,
            )

        return self._default_dgridparams

    def parse_dgridparams(self, dgridparams):
        """Parse depth gridding parameters.

        See top level class notes for more info.

        Parameters
        ----------
        dgridparams : dict

        """
        self.dgridparams = Bunch(self.default_dgridparams)
        if dgridparams is not None:
            self.dgridparams.update_values(dgridparams, strict=True)
        else:
            logger.warning(
                "No depth gridding parameters provided, using default values."
            )
        # Default depth grid parameters are based on the median pressure to
        # avoid binning into neighboring grid cells as much as possible.
        # Therefore, we start assembling the depth grid from the bottom up for
        # an uplooker and from the top down for a downlooker.
        if self.orientation == "up":
            self.dgrid = np.arange(
                self.dgridparams.dbot,
                self.dgridparams.dtop,
                -self.dgridparams.d_interval,
                dtype=float,
            )
        elif self.orientation == "down":
            self.dgrid = np.arange(
                self.dgridparams.dtop,
                self.dgridparams.dbot,
                self.dgridparams.d_interval,
                dtype=float,
            )

    def parse_tgridparams(self, tgridparams):
        """Parse time gridding parameters.

        See top level class notes for more info.

        Parameters
        ----------
        tgridparams : dict

        """
        # Find time at depth to determine default time grid parameters.
        # Differentiate between time series only in the water and time series
        # including the overshoot on mooring deployment.
        p = self.tsdat.pressure
        if ~np.any(p < 10):
            t0 = self.dday[0]
            t1 = self.dday[-1]
        else:
            at_depth = np.nonzero(p > self.p_median)[0][0]
            t0 = self.dday[at_depth]
            in_water = np.nonzero(p > self.p_median / 2)[0][-1]
            t1 = self.dday[in_water]

        # Generate a set of default time gridding parameters and then update
        # from the input parameters provided.
        default_tgridparams = dict(dt_hours=0.5, t0=t0, t1=t1, burst_average=False)
        self.tgridparams = Bunch(default_tgridparams)
        if tgridparams is not None:
            self.tgridparams.update_values(tgridparams, strict=True)
        else:
            logger.warning(
                "No time gridding parameters provided, using default values."
            )

    def parse_editparams(self, editparams):
        """Parse editing parameters.

        See top level class notes for more info.

        Parameters
        ----------
        editparams : dict

        """
        self.editparams = Bunch(self._editparams)
        if editparams is not None:
            self.editparams.update_values(editparams, strict=True)
        else:
            logger.warning("No edit parameters provided, using default values.")

    def parse_driftparams(self, driftparams):
        """Parse time drift parameters.

        See top level class notes for more info.

        Parameters
        ----------
        driftparams : dict

        """
        driftparams = dict() if driftparams is None else driftparams
        self.driftparams = driftparams
        self.yearbase = self.m.yearbase
        t0 = self.tsdat.dday[0]
        self.t0 = t0
        t1_adcp = driftparams.get("end_adcp", None)
        if t1_adcp is not None:
            t1_pc = to_day(self.m.yearbase, *driftparams["end_pc"])
            t1_adcp = to_day(self.m.yearbase, *driftparams["end_adcp"])
            self.time_drift_rate = (t1_pc - t0) / (t1_adcp - t0)
        else:
            logger.warning(
                "No time drift parameters provided, not applying any clock correction."
            )
            self.time_drift_rate = 1

        self.dday = self._correct_dday(self.tsdat.dday)

    def _correct_dday(self, dday_orig):
        """Apply linear correction for clock drift.

        Parameters
        ----------
        dday_orig : array-like
            Origial time vector.

        Returns
        -------
        array-like
            Corrected time vector

        """

        return self.t0 + self.time_drift_rate * (dday_orig - self.t0)

    def _parse_sysconfig(self):
        # We have to get the up/down reading from sysconfig for a time when the
        # instrument was in the water. Look at pressure and determine ensembles
        # during time at depth.
        ii = np.flatnonzero(self.tsdat.pressure > 15)
        depth_ii = (ii[0] + ii[-1]) // 2
        try:
            self.tsdat.pressure[depth_ii] > 15
        except ValueError:
            print("could not determine ensemble index at depth for reading sysconfig")
        at_depth = self.m.read(
            varlist=["VariableLeader"], start=depth_ii, stop=depth_ii + 1
        )
        self.orientation = "up" if at_depth.sysconfig.up else "down"
        self.sysconfig = at_depth.sysconfig

    def make_start_ddays(self):
        """Generate time stamps for ping averaging.

        Notes
        -----
        If turning on burst averaging, other input values will be ignored and
        the averaging interval will be determined from the burst sampling
        scheme apparent in the ping pattern.

        """
        dday_start = self.tgridparams.t0
        dday_end = self.tgridparams.t1
        dt_hours = self.tgridparams.dt_hours
        is_burst_average = self.tgridparams.burst_average

        # Save whether we are averaging over bursts or not.
        self.is_burst_average = is_burst_average
        # Generate time stamps and stuff.
        if not is_burst_average:
            print("no burst average")
            self.dday_start = dday_start
            self.dday_end = dday_end
            self.dt = dt_hours / 24.0
            self.start_ddays = np.arange(dday_start, dday_end, dt_hours / 24.0)
            # Time stamps for the averages. Midpoints of averaging intervals.
            self.dday_mid = self.start_ddays + self.dt / 2
        else:
            dday_diff = np.diff(self.dday)
            # Determine ping interval within burst and time between bursts.
            burst_dt = np.median(dday_diff)
            print(f"time between pings within burst: {burst_dt * 24 * 60 * 60:1.1f} s")
            # It seems safe to assume that the time between bursts is at least
            # four times as long as the time between individual pings within a
            # burst.
            inter_burst_dt = np.median(dday_diff[dday_diff > burst_dt * 4])
            print(f"time between bursts: {inter_burst_dt * 24 * 60:1.1f} min")
            # Find starting points of all bursts.
            start_indices = np.flatnonzero(dday_diff > burst_dt * 4)
            # Increase index so we are at the end of the larger time differences.
            start_indices += 1
            # Include the very beginning.
            start_indices = np.insert(start_indices, 0, 0)
            self.start_ddays = self.dday[start_indices]
            self.dday_start = self.start_ddays[0]
            # Generate a dt that is inclusive of one burst.  We know the number
            # of pings in a burst from the difference between the
            # start_indices. Let's go a bit beyond the time needed (3 more ping
            # intervals chosen here).
            pings_per_burst = np.int32(np.median(np.diff(start_indices)))
            print(f"{pings_per_burst} pings per burst")
            print(f"{start_indices.shape[0]} bursts total")
            self.dt = burst_dt * pings_per_burst + burst_dt * 3

            self.dday_start = self.start_ddays[0]
            self.dday_end = dday_end

            # Time stamps in the middle of the burst
            self.dday_mid = self.start_ddays + pings_per_burst * burst_dt / 2

    def read_ensemble(self, iens):
        """Read ensembles (several individual pings grouped together).

        Parameters
        ----------
        iens : int
            Index into ensemble start times (they are generated in
            :meth:`make_start_ddays`).

        Returns
        -------
        dat : pycurrents.adcp.rdiraw.Bunch
            Dictionary with data.

        Raises
        ------
        ValueError
            If `iens` is too large to index into `start_ddays`.

        """

        if iens > len(self.start_ddays) - 1:
            raise ValueError("ens num %d is out of range" % iens)

        # Get indices within the interval.
        sl = rangeslice(
            self.dday, self.start_ddays[iens], self.start_ddays[iens] + self.dt
        )
        # Use the indices to extract data.
        dat = self.m.read(start=sl.start, stop=sl.stop)
        if dat is None:
            return None
        dat.dday_orig = dat.dday
        dat.dday = self._correct_dday(dat.dday_orig)
        if self._pressure_provided is not None:
            # Replace pressure if provided
            dat.pressure = self._external_pressure_to_dat(dat)
        else:
            dat.pressure = self._scale_pycurrents_pressure(dat)
        sign = -1 if self.orientation == "up" else 1
        pdepth = seawater.depth2(dat.pressure, self.lat)
        dat.depth = pdepth[:, np.newaxis] + sign * dat.dep
        return dat

    @property
    def magdec(self):
        """Magnetic declination.

        If not provided as input argument magdec is calculated using
        [magdec](https://currents.soest.hawaii.edu/hgstage/geomag/file/tip)
        (must be installed) based on `lon` and `lat`.

        """
        if self._magdec is None:
            if self.lat is None:
                logger.warning(
                    "No lon/lat provided, cannot calculate magnetic declination."
                )
                self._magdec = 0
            else:
                # Look for magdec executable
                magdec_found = True
                magdec_path = which("magdec")

                if magdec_path is None:
                    magdec_found = False
                    package_dir = os.path.dirname(__file__)

                if not magdec_found:
                    # Try this package directory
                    magdec_path = os.path.join(package_dir, "magdec")
                    magdec_found = os.path.isfile(magdec_path)

                if not magdec_found:
                    # Try the magdec installation directory
                    magdec_path = os.path.abspath(
                        os.path.join(package_dir, "../geomag/magdec")
                    )
                    magdec_found = os.path.isfile(magdec_path)

                if not magdec_found:
                    raise FileNotFoundError(
                        "Cannot find program magdec on the system path or paths within velosearaptor."
                    )

                logger.info(f"magdec found at {magdec_path}")

                n = len(self.start_ddays)
                dday_mid = self.start_ddays[n // 2]
                y, m, d = to_date(self.yearbase, dday_mid)[:3]

                output = Popen(
                    [
                        magdec_path,
                        str(self.lon),
                        str(self.lat),
                        str(y),
                        str(m),
                        str(d),
                    ],
                    stdout=PIPE,
                ).communicate()[0]
                output = output.strip()
                logger.info("magdec output is: %s", output)
                self._magdec = float(output.split()[0])
        elif self._magdec_provided is not None:
            logger.info(f"magdec {self._magdec_provided} provided.")
        return self._magdec

    @property
    def raw(self):
        """Raw ADCP data."""
        if self._raw is None:
            print("Reading raw data...")
            self._raw = io.read_raw_rdi(self.files)
            self._raw.coords["bin"] = (("z"), np.arange(self._raw.z.size))
        return self._raw

    def _edit(self, ens):
        """Apply editing to xyze."""
        ep = self.editparams
        cond = (ens.cor < ep.min_correlation).any(axis=-1)
        ens.xyze[cond] = np.ma.masked
        e = ens.xyze[:, :, 3]
        max_e = min(ep.max_e, e.std() * ep.max_e_deviation)
        ens.max_e_applied = max_e
        cond = np.abs(e) > max_e
        ens.xyze[cond] = np.ma.masked
        if ep.maskbins is not None:
            ens.xyze[:, ep.maskbins, :] = np.ma.masked

    def _to_enu(self, ens):
        """
        add enu
        GV: enu is east, north, up, errvel (optional) whereas xyz are
        instrument coordinates.
        """
        ens.enu = rdi_xyz_enu(
            ens.xyze,
            ens.heading + self.magdec,
            ens.pitch,
            ens.roll,
            orientation=self.orientation,
        )

    def _burst_average_depth(self, ens):
        """Depth-average within a burst.

        Average the depth vectors if doing burst-averages. Otherwise we just
        return all depth vectors of this ensemble.

        Parameters
        ----------
        ens : Bunch
            Ensemble data.

        Returns
        -------
        depth : array-like
            Depth vector for each ping.

        """
        if self.is_burst_average:
            depth_mean = ens.depth.mean(axis=0)
            depth = np.tile(depth_mean, (ens.dday.size, 1))
        else:
            depth = ens.depth
        return depth

    def _regrid_enu(self, ens, method="linear"):
        """Depth-grid enu velocities."""
        shape = (ens.dday.size, self.dgrid.size, ens.enu.shape[-1])
        enu_grid = np.ma.zeros(shape)
        enu_grid[:] = np.ma.masked
        # Average the depth vectors if doing burst-averages. Otherwise we just
        # return all depth vectors of this ensemble.
        # TO DO: If calculating averages over regular time intervals, we need
        # to low pass filter the pressure time series prior to creating the
        # depth vectors.
        depth = self._burst_average_depth(ens)
        for i in range(ens.dday.size):
            enu_grid[i] = interp1(
                depth[i], ens.enu[i], self.dgrid, axis=0, method=method
            )
        ens.enu_grid = enu_grid

    def _regrid_amp(self, ens, method="linear"):
        """Depth-grid amplitudes (averaged over all 4 beams)."""
        shape = (ens.dday.size, self.dgrid.size)
        amp_grid = np.ma.zeros(shape)
        amp_grid[:] = np.ma.masked
        depth = self._burst_average_depth(ens)
        for i in range(ens.dday.size):
            amp_grid[i] = interp1(
                depth[i],
                ens.amp[i].mean(axis=-1),
                self.dgrid,
                axis=0,
                method=method,
            )
        ens.amp_grid = amp_grid

    def _binmap_one_beam(self, ens, beam_number):
        """Binmap single ping data for a single beam by linear interpolation.

        Mapping is applied to velocity, amplitude and correlation.

        Currently only works for 4 beam ADCP.

        Parameters
        ----------
        ens : Bunch
            An ADCP dataset read by Multiread.
        beam_number : int
            An integer from 1 to 4 representing the beam number.

        Returns
        -------
        veli : ndarray
            Mapped velocity.
        ampi : ndarray
            Mapped amplitude.
        cori : ndarray
            Mapped correlation.

        """

        if beam_number not in [1, 2, 3, 4]:
            raise ValueError("Beam number must be 1, 2, 3 or 4.")

        tba = np.tan(np.deg2rad(ens.sysconfig.angle))  # Tangent of beam angle
        pitch = np.deg2rad(ens.pitch)
        roll = np.deg2rad(ens.roll)

        # The true bin distances
        if beam_number == 1:
            dep = (
                ens.dep[None, :]
                * ((np.cos(roll) - tba * np.sin(roll)) * np.cos(pitch))[:, None]
            )  # None adds a new axis.
        elif beam_number == 2:
            dep = (
                ens.dep[None, :]
                * ((np.cos(roll) + tba * np.sin(roll)) * np.cos(pitch))[:, None]
            )
        elif beam_number == 3:
            dep = (
                ens.dep[None, :]
                * ((np.cos(pitch) + tba * np.sin(pitch)) * np.cos(roll))[:, None]
            )
        elif beam_number == 4:
            dep = (
                ens.dep[None, :]
                * ((np.cos(pitch) - tba * np.sin(pitch)) * np.cos(roll))[:, None]
            )

        vel = ens.vel[..., beam_number - 1]
        amp = ens.amp[..., beam_number - 1]
        cor = ens.cor[..., beam_number - 1]

        # Calculate interpolating weights (this hogs RAM!)
        dz = np.diff(dep, axis=1)
        w = np.clip((ens.dep - dep[:, :-1, None]) / dz[:, :, None], 0, 1)

        # Determine data above or below the deepest bins
        above = (w == 1.0).all(axis=1)
        below = (w == 0.0).all(axis=1)
        bad = above | below

        # Calculate differences
        dvel = np.diff(vel, axis=1)
        damp = np.diff(amp, axis=1)
        dcor = np.diff(cor, axis=1)

        veli = vel[:, [0]] + np.sum(w * dvel[:, :, None], axis=1)
        veli[bad] = np.nan

        ampi = amp[:, [0]] + np.sum(w * damp[:, :, None], axis=1)
        ampi[bad] = np.nan

        cori = cor[:, [0]] + np.sum(w * dcor[:, :, None], axis=1)
        cori[bad] = np.nan

        return veli, ampi, cori

    def _binmap_all_beams(self, ens):
        """Binmap single ping data for all beams."""

        for beam_number in [1, 2, 3, 4]:
            veli, ampi, cori = self._binmap_one_beam(ens, beam_number)
            ens.vel[..., beam_number - 1] = veli
            ens.amp[..., beam_number - 1] = ampi
            ens.cor[..., beam_number - 1] = cori

            ens[f"vel{beam_number}"] = veli
            ens[f"amp{beam_number}"] = ampi
            ens[f"cor{beam_number}"] = cori

    def _calculate_xyze(self, ens, ibad=None):
        """Calculate xyze from along-beam data."""
        if ens.sysconfig.convex:
            geom = "convex"
        else:
            geom = "concave"

        trans = Transform(angle=ens.sysconfig.angle, geometry=geom)

        ens.xyze = trans.beam_to_xyz(ens.vel, ibad=ibad)

    def process_pings(self, start=None, stop=None, binmap=False, ens_size=50000):
        """Process single ping data without averaging.

        Adds results as dictionary under `ave` and as `xarray.Dataset` under `ds`.

        Writes processing parameters to the log file.

        Parameters
        ----------
        start : int, optional
            Start processing at this ping number.
        stop : int, optional
            Stop processing at this ping number.
        binmap : bool, optional
            Do binmapping of along-beam data.
        ens_size : int, optional
            Pings are processed in ensembles to reduce memory usage.
            This parameter sets how many pings are in an ensemble. The default is 50000.

        """
        if start is None:
            idx_start = np.searchsorted(self.dday, self.dday_start)
        if stop is None:
            idx_stop = np.searchsorted(self.dday, self.dday_end)

        ens_idxs = np.hstack((np.arange(idx_start, idx_stop, ens_size), idx_stop))
        write_idxs = ens_idxs - ens_idxs[0]  # Arrays we write to start at index 0
        npings = idx_stop - idx_start
        nens = ens_idxs.size - 1
        ndgrid = self.tsdat.dep.size

        logger.info("Processing all pings")
        logger.info(f"Binmapping is {binmap}")

        uvwe = np.ma.zeros((npings, ndgrid, 4), dtype=np.float32)

        pg = np.zeros((npings, ndgrid), dtype=np.int8)
        amp = np.ma.zeros((npings, ndgrid), dtype=np.float32)

        # temperature = np.ma.zeros((npings,), dtype=np.float32)
        # pressure = np.ma.zeros((npings,), dtype=np.float32)

        temperature = self.tsdat.temperature[idx_start:idx_stop]
        pressure = self.tsdat.pressure[idx_start:idx_stop]

        dday = self.dday[idx_start:idx_stop]

        # Loop over ensembles
        for i in tqdm(range(nens)):
            ens = self.m.read(start=ens_idxs[i], stop=ens_idxs[i + 1])
            idx0 = write_idxs[i]
            idx1 = write_idxs[i + 1]

            if ens is not None:
                # I pulled this out of read_ensembles because we need to calculate depth for _ave2nc to work.
                ens.dday_orig = ens.dday
                ens.dday = self._correct_dday(ens.dday_orig)
                if self._pressure_provided is not None:
                    # Replace pressure if provided
                    ens.pressure = self._external_pressure_to_dat(ens)
                else:
                    ens.pressure = self._scale_pycurrents_pressure(ens)
                sign = -1 if self.orientation == "up" else 1
                pdepth = seawater.depth2(ens.pressure, self.lat)
                ens.depth = pdepth[:, np.newaxis] + sign * ens.dep

                if binmap:
                    self._binmap_all_beams(ens)
                    # Now we have to recalculate xyze with the binmapped data.
                    self._calculate_xyze(ens)

                self._edit(ens)  # modifies xyze
                self._to_enu(ens)  # transform to earth coords (east, north, up)

            else:
                uvwe[idx0:idx1] = np.ma.masked
                amp[idx0:idx1] = np.ma.masked
                # pressure[idx0:idx1] = np.ma.masked
                # temperature[idx0:idx1] = np.ma.masked
                continue

            uvwe[idx0:idx1] = ens.enu

            # pgi = 100 * ens.enu_grid[..., 0].count(axis=0) // nprofs
            # pg[i] = pgi.astype(np.int8)
            amp[idx0:idx1] = ens.amp.mean(axis=-1)  # Average over beams... why?

        self.ave = Bunch(
            u=uvwe[..., 0],
            v=uvwe[..., 1],
            w=uvwe[..., 2],
            e=uvwe[..., 3],
            pg=pg,
            amp=amp,
            temperature=temperature,
            pressure=pressure,
            # npings=npings,
            dday=dday,
            yearbase=self.yearbase,
            dep=ens.dep,  # <<<<----- BAD!!! This a fudge because I don't want to calculated a depth vector. We use the depth vector from the last ensemble.
            editparams=self.editparams,
            tgridparams=self.tgridparams,
            # dgridparams=self.dgridparams,
            magdec=self.magdec,
            lon=self.lon,
            lat=self.lon,
        )

        self._ave2nc()
        self._add_meta_data_to_ds()
        self._log_processing_params()

    def average_ensembles(self, start=None, stop=None):
        """Time-averaging.

        Adds results as dictionary under `ave` and as `xarray.Dataset` under `ds`.

        Writes processing parameters to the log file.

        Parameters
        ----------
        start : int
            Range start for averaging. Index into start times of averaging
            intervals.
        stop : int
            Range start for averaging. Index into start times of averaging
            intervals.

        """

        nens_orig = len(self.start_ddays)
        indices_orig = np.arange(nens_orig)
        indices = indices_orig[start:stop]
        if start is None and stop is None:
            logger.info("Averaging all ensembles")
        else:
            logger.info(f"Averaging ensembles {indices[0]} to {indices[-1]}")
        nens = len(indices)
        ndgrid = len(self.dgrid)
        uvwe = np.ma.zeros((nens, ndgrid, 4), dtype=np.float32)
        uvwe_std = np.ma.zeros((nens, ndgrid, 4), dtype=np.float32)

        pg = np.zeros((nens, ndgrid), dtype=np.int8)
        amp = np.ma.zeros((nens, ndgrid), dtype=np.float32)

        temperature = np.ma.zeros((nens,), dtype=np.float32)
        pressure = np.ma.zeros((nens,), dtype=np.float32)
        pressure_std = np.ma.zeros((nens,), dtype=np.float32)
        pressure_max = np.ma.zeros((nens,), dtype=np.float32)

        npings = np.zeros((nens,), dtype=np.int16)

        dday = self.dday_mid[start:stop]

        for i, iens in enumerate(tqdm(indices)):
            ens = self.read_ensemble(iens)
            if ens is not None:
                self._edit(ens)  # modifies xyze
                self._to_enu(ens)  # transform to earth coords (east, north, up)
                self._regrid_enu(ens)
                self._regrid_amp(ens)

                nprofs = ens.enu_grid.shape[0]
            else:
                nprofs = 0
            npings[i] = nprofs
            if nprofs < 2:
                uvwe[i] = np.ma.masked
                uvwe_std[i] = np.ma.masked
                # (pg is not a masked array)
                amp[i] = np.ma.masked
                pressure[i] = np.ma.masked
                pressure_std[i] = np.ma.masked
                pressure_max[i] = np.ma.masked
                temperature[i] = np.ma.masked
                continue

            uvwe[i] = ens.enu_grid.mean(axis=0)
            uvwe_std[i] = ens.enu_grid.std(axis=0)

            pgi = 100 * ens.enu_grid[..., 0].count(axis=0) // nprofs
            pg[i] = pgi.astype(np.int8)
            amp[i] = ens.amp_grid.mean(axis=0)

            pressure[i] = ens.pressure.mean()
            pressure_std[i] = ens.pressure.std()
            pressure_max[i] = ens.pressure.max()
            temperature[i] = ens.temperature.mean()

        self.ave = Bunch(
            u=uvwe[..., 0],
            v=uvwe[..., 1],
            w=uvwe[..., 2],
            e=uvwe[..., 3],
            u_std=uvwe_std[..., 0],
            v_std=uvwe_std[..., 1],
            w_std=uvwe_std[..., 2],
            e_std=uvwe_std[..., 3],
            pg=pg,
            amp=amp,
            temperature=temperature,
            pressure=pressure,
            pressure_std=pressure_std,
            pressure_max=pressure_max,
            npings=npings,
            dday=dday,
            yearbase=self.yearbase,
            dep=self.dgrid,
            editparams=self.editparams,
            tgridparams=self.tgridparams,
            dgridparams=self.dgridparams,
            magdec=self.magdec,
            lon=self.lon,
            lat=self.lon,
        )

        self._ave2nc()
        self._add_meta_data_to_ds()
        self._log_processing_params()

    def burst_average_ensembles(self, start=None, stop=None, interpolate_bin=None):
        """Time-averaging prior to depth-gridding.

        Uses pre-defined editing parameters that can be updated with
        `parse_editparams`.

        Adds results as dictionary under `ave` and as `xarray.Dataset` under `ds`.

        Writes processing parameters to the log file.

        Parameters
        ----------
        start : int, optional
            Range start for averaging. Index into start times of averaging
            intervals. Defaults to None (start at beginning).
        stop : int, optional
            Range start for averaging. Index into start times of averaging
            intervals. Defaults to None (start at beginning).
        interpolate_bin : int or None, optional
            Interpolate over a single, previously masked, bin. Defaults to None (no interpolation).

        """
        pg_condition = self.editparams.pg_limit
        nens_orig = len(self.start_ddays)
        indices_orig = np.arange(nens_orig)
        indices = indices_orig[start:stop]
        if start is None and stop is None:
            logger.info("Averaging all ensembles")
        else:
            logger.info(f"Averaging ensembles {indices[0]} to {indices[-1]}")
        nens = len(indices)
        ndgrid = len(self.dgrid)
        uvwe = np.ma.zeros((nens, ndgrid, 4), dtype=np.float32)
        uvwe_std = np.ma.zeros((nens, ndgrid, 4), dtype=np.float32)

        pg = np.zeros((nens, ndgrid), dtype=np.int8)
        amp = np.ma.zeros((nens, ndgrid), dtype=np.float32)

        temperature = np.ma.zeros((nens,), dtype=np.float32)
        pressure = np.ma.zeros((nens,), dtype=np.float32)
        pressure_std = np.ma.zeros((nens,), dtype=np.float32)
        pressure_max = np.ma.zeros((nens,), dtype=np.float32)

        npings = np.zeros((nens,), dtype=np.int16)

        dday = self.dday_mid[start:stop]

        for i, iens in enumerate(tqdm(indices)):
            ens = self.read_ensemble(iens)
            if ens is not None:
                self._edit(ens)  # modifies xyze
                self._to_enu(ens)  # transform to earth coords (east, north, up)
                self._regrid_enu(ens)
                self._regrid_amp(ens)

                nprofs = ens.enu_grid.shape[0]
            else:
                nprofs = 0
            npings[i] = nprofs
            if nprofs < 2:
                uvwe[i] = np.ma.masked
                uvwe_std[i] = np.ma.masked
                # (pg is not a masked array)
                amp[i] = np.ma.masked
                pressure[i] = np.ma.masked
                pressure_std[i] = np.ma.masked
                pressure_max[i] = np.ma.masked
                temperature[i] = np.ma.masked
                continue

            # Depth vector for interpolation
            depth = self._burst_average_depth(ens)
            depth = depth[0, :]

            # Calculate burst-average on instrument-relative depth grid
            uvwe_inst = ens.enu.mean(axis=0)
            uvwe_std_inst = ens.enu.std(axis=0)

            pgi_inst = 100 * ens.enu[..., 0].count(axis=0) // nprofs

            if pg_condition is not None:
                pgi_index = pgi_inst < pg_condition
                uvwe_inst[pgi_index, :] = np.ma.masked
                uvwe_std_inst[pgi_index, :] = np.ma.masked

            if interpolate_bin is not None:
                zi = interpolate_bin
                neighbors = [zi - 2, zi - 1, zi + 1, zi + 2]
                dtmp = depth[neighbors]
                tmp = interp1(
                    dtmp,
                    uvwe_inst[neighbors, :],
                    depth[zi],
                    axis=0,
                    method="linear",
                )
                uvwe_inst[zi, :] = tmp

            # Interpolate burst-average to universal depth grid.
            uvwe_grid = interp1(depth, uvwe_inst, self.dgrid, axis=0, method="linear")
            uvwe_std_grid = interp1(
                depth, uvwe_std_inst, self.dgrid, axis=0, method="linear"
            )
            uvwe[i] = uvwe_grid
            uvwe_std[i] = uvwe_std_grid

            # Interpolate pg to universal depth grid. Not overly satisfying but
            # seems like that's what we need to do here.
            pgi_grid = interp1(depth, pgi_inst, self.dgrid, axis=0, method="linear")
            pg[i] = pgi_grid.astype(np.int8)

            # Not changed to averaging in instrument-relative coordinates first.
            amp[i] = ens.amp_grid.mean(axis=0)

            pressure[i] = ens.pressure.mean()
            pressure_std[i] = ens.pressure.std()
            pressure_max[i] = ens.pressure.max()
            temperature[i] = ens.temperature.mean()

        self.ave = Bunch(
            u=uvwe[..., 0],
            v=uvwe[..., 1],
            w=uvwe[..., 2],
            e=uvwe[..., 3],
            u_std=uvwe_std[..., 0],
            v_std=uvwe_std[..., 1],
            w_std=uvwe_std[..., 2],
            e_std=uvwe_std[..., 3],
            pg=pg,
            amp=amp,
            temperature=temperature,
            pressure=pressure,
            pressure_std=pressure_std,
            pressure_max=pressure_max,
            npings=npings,
            dday=dday,
            yearbase=self.yearbase,
            dep=self.dgrid,
            editparams=self.editparams,
            tgridparams=self.tgridparams,
            dgridparams=self.dgridparams,
            magdec=self.magdec,
            lon=self.lon,
            lat=self.lon,
        )

        self._ave2nc()
        self._add_meta_data_to_ds()
        self._log_processing_params()

    def _safely_add_attribute_from_params(self, key, d):
        """Add attribute from dictionary and throw an exception if the key is missing.

        Parameters
        ----------
        key : str
        d : dict
        """
        try:
            setattr(self, key, d[key])
        except KeyError as error:
            print(f"{error} is missing in input parameters")
            raise

    def _set_up_logger(self):
        """Set up logging to both a file and to screen.

        Default for screen logging is to show only warnings unless `verbose` is
        set to True.

        """

        # Parse metadata for naming the log.
        if self.meta_data is not None:
            if "sn" in self.meta_data and "mooring" in self.meta_data:
                sn = self.meta_data["sn"]
                mooring = self.meta_data["mooring"]
                has_mooring_id = True
            else:
                sn = None
                mooring = None
                has_mooring_id = False
        else:
            has_mooring_id = False

        # Set up a logger that will collect log info from this module and the
        # other pycurrent methods as well.
        logdir = Path(self.logdir)
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
        if self.verbose:
            ConsoleOutputHandler.setLevel(logging.INFO)
        logger.addHandler(ConsoleOutputHandler)

        # current date
        datestr = np.datetime64(datetime.datetime.now()).astype(datetime.datetime)
        strformat = "%Y-%m-%d"
        datestr = datestr.strftime(strformat)

        if has_mooring_id:
            logger.info(f"Processing {mooring} SN {sn} on {datestr}")
        else:
            logger.warning(
                "No meta data provided, logging to generic filename 'adcp_proc.log'"
            )
            logger.info(f"Processing on {datestr}")

    def _log_processing_params(self):
        """Write processing parameters to log file."""

        logger.info("processing settings")
        logger.info("-------------------")
        self._log_params(self.dgridparams)
        self._log_params(self.tgridparams)
        self._log_params(self.editparams)
        logger.info("-------------------")

    def _log_params(self, pd):
        """Print ADCP processing parameter dict to logs.

        Parameters
        ----------
        pd : dict
            Parameter dict.
        """

        for k, v in pd.items():
            if k == "maskbins" and v is not None:
                logstr = (k, ":", np.flatnonzero(v))
                logger.info(logstr)
            else:
                logstr = (k, ":", v)
                logger.info(logstr)

    def _ave2nc(self):
        """Convert data structure from ave to xarray.Dataset format."""
        # load npz file
        dat = self.ave.copy()
        dat = Bunch(dat)
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

        # Drop depth levels with all nan
        out = out.dropna(how="all", dim="z")

        # add variable names and units for plotting
        out = self._add_names_and_units(out)

        self.ds = out

    def _add_names_and_units(self, ds):
        """Add variable meta-data to Dataset.

        Parameters
        ----------
        ds : xarray.Dataset

        Returns
        -------
        ds : xarray.Dataset

        """

        ds.u.attrs = dict(long_name="u", units="m/s")
        ds.v.attrs = dict(long_name="v", units="m/s")
        ds.w.attrs = dict(long_name="w", units="m/s")
        ds.e.attrs = dict(long_name="error velocity", units="m/s")
        ds.z.attrs = dict(long_name="depth", units="m")
        ds.temperature.attrs = dict(long_name="temperature", units="°C")
        ds.pressure.attrs = dict(long_name="pressure", units="dbar")
        return ds

    def _add_meta_data_to_ds(self):
        # Add some more info.
        self.ds.attrs["orientation"] = self.orientation
        self.ds.attrs["magdec"] = self.magdec
        for att in ["max_e", "max_e_deviation", "min_correlation"]:
            self.ds.attrs[att] = self.editparams[att]

        # Add meta data if provided.
        if self.meta_data is not None:
            for k, v in self.meta_data.items():
                self.ds.attrs[k] = v
        self.ds.attrs["proc time"] = np.datetime64("now").astype("str")

        # Calculate transducer depth from pressure
        self.ds["xducer_depth"] = -gsw.z_from_p(self.ds.pressure, self.lat)
        self.ds.xducer_depth.attrs = dict(long_name="transducer depth", units="m")

    def plot_echo_stats(self):
        """Plot beam statistics (correlation and amplitude) from raw ADCP data."""
        r = self.raw

        fig, ax = plt.subplots(
            nrows=1,
            ncols=2,
            figsize=(5, r.bin.max().data * 0.3),
            constrained_layout=True,
            sharey=True,
        )
        r.cor.mean(dim="time").plot(
            hue="beam", y="bin", marker="o", linestyle="", ax=ax[0]
        )
        r.amp.mean(dim="time").plot(
            hue="beam", y="bin", marker="o", linestyle="", ax=ax[1]
        )
        ax[0].invert_yaxis()
        ax[1].set(ylabel="")
        ax[0].set_yticks(r.bin.data)
        for axi in ax:
            axi.grid(True)

    def plot_pressure(self):
        """Plot pressure time series and mark time at depth."""
        fig, ax = plt.subplots(
            nrows=1,
            ncols=1,
            figsize=(6, 2.5),
            constrained_layout=True,
        )
        self.raw.pressure.plot(ax=ax, label="all")
        self.raw.pressure.where(self.raw.pressure > 50).plot(ax=ax, label="subsurface")
        ax.invert_yaxis()
        ax.set(xlabel="", ylabel="pressure [dbar]")
        ax.legend()

    def generate_binmask(self, indices):
        binmask = self.raw.bin.data < 0
        binmask[indices] = True
        return binmask


class ProcessADCPyml(ProcessADCP):
    """Moored ADCP processing with parameters provided via .yml file.

    An example parameter file is included at
    [`notebooks/parameters.yml`](https://github.com/modscripps/velosearaptor/tree/main/notebooks/parameters.yml)"""

    def __init__(self, parameter_file, mooring, sn, **kwargs):
        p = io.parse_yaml_input(parameter_file, "MAVS2", 24606)
        super().__init__(
            p["data_dir"],
            meta_data=p["meta_data"],
            driftparams=p["driftparams"],
            tgridparams=p["tgridparams"],
            dgridparams=p["dgridparams"],
            editparams=p["editparams"],
            **kwargs,
        )
