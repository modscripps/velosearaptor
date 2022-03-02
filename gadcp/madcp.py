#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module gadcp.madcp with functions for moored ADCPs."""

import os
from subprocess import Popen, PIPE  # for magdec
import logging

import datetime
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import xarray as xr
import pathlib
from pathlib import Path
from tqdm.notebook import tqdm

from pycurrents.adcp.rdiraw import extract_raw, Multiread
from pycurrents.system import Bunch
from pycurrents.system import Bunch
from pycurrents.num.nptools import rangeslice
from pycurrents.num import interp1
from pycurrents.codas import to_date, to_day
from pycurrents.adcp.transform import rdi_xyz_enu
from pycurrents.file import npzfile
from pycurrents.data import seawater

# from gadcp.mcm_avg import MCM, Pingavg
from . import io

import gvpy as gv

# Standard logging
logger = logging.getLogger(__name__)


class ProcessADCP:
    """Docstring for ProcessADCP.

    Attributes
    ----------
    files : list
        List pointing to raw data file(s).

    """

    _editparams = dict(
        max_e=0.2,  # absolute max e
        max_e_deviation=2,  # max in terms of sigma
        min_correlation=64,  # 64 is RDI default
        maskbins=None,  # do not mask any bins
    )

    def __init__(
        self,
        raw_data,
        meta_data,
        driftparams,
        tgridparams,
        dgridparams,
        editparams,
        ibad=None,
        logdir="log",
        verbose=False,
        plot=False,
        pressure_scale_factor=1,
    ):
        """TODO: to be defined."""

        self.meta_data = meta_data.copy()
        self.driftparams = driftparams
        self.editparams = editparams
        self.ibad = ibad
        self.logdir = logdir
        self.verbose = verbose

        self.pressure_scale_factor = pressure_scale_factor
        self._magdec = None
        self._raw = None

        self._set_up_logger()
        self.parse_file_locations(raw_data)
        self.parse_meta_data()
        self.generate_data_reader()
        self.parse_dgridparams(dgridparams)
        self.parse_tgridparams(tgridparams)
        self.parse_editparams(editparams)
        # self._log_processing_params()

        self.make_start_ddays()

        if plot:
            self.plot_pressure()

    def parse_file_locations(self, raw_data, min_file_size=1e4):
        """Parse input for raw data files.

        Input can either be a single file name as a str, a single file as a
        Path object, a list of either of these, or a Path object pointing to a
        directory with raw ADCP files. In the latter case, files that are
        smaller than a threshold will not be included in the processing.

        Adds attribute `files`.

        Parameters
        ----------
        raw_data : str or list or Path
            Location(s) of raw data.
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

    def parse_meta_data(self):
        """Add essential meta data to attributes and remove them from the meta_data dict."""
        essential_meta_data = ["lon", "lat"]
        [
            self._safely_add_attribute_from_params(k, self.meta_data)
            for k in essential_meta_data
        ]

    def parse_dgridparams(self, dgridparams):
        self.p_median = np.median(self.tsdat.pressure)
        if self.sysconfig["up"]:
            default_dgridparams = dict(
                dbot=int(self.p_median), dtop=10, d_interval=5
            )
        else:
            default_dgridparams = dict(
                dtop=int(self.p_median),
                dbot=int(self.p_median) + 1000,
                d_interval=5,
            )
        self.dgridparams = Bunch(default_dgridparams)
        if dgridparams is not None:
            self.dgridparams.update_values(dgridparams, strict=True)
        else:
            logger.warning(
                "No depth gridding parameters provided, using default values."
            )
        self.dgrid = np.arange(
            self.dgridparams.dtop,
            self.dgridparams.dbot,
            self.dgridparams.d_interval,
            dtype=float,
        )

    def parse_tgridparams(self, tgridparams):
        # Find time at depth.
        p = self.tsdat.pressure
        at_depth = np.nonzero(p > self.p_median)[0][0]
        t0 = (2 + np.ceil(self.dday[at_depth] * 24)) / 24.0
        in_water = np.nonzero(p > self.p_median / 2)[0][-1]
        t1 = (np.floor(self.dday[in_water] * 24) - 2) / 24.0

        # Generate a set of default time gridding parameters and then update
        # from the input parameters provided.
        default_tgridparams = dict(
            dt_hours=0.5, t0=t0, t1=t1, burst_average=False
        )
        self.tgridparams = Bunch(default_tgridparams)
        if tgridparams is not None:
            self.tgridparams.update_values(tgridparams, strict=True)
        else:
            logger.warning(
                "No time gridding parameters provided, using default values."
            )

    def parse_editparams(self, editparams):
        self.editparams = Bunch(self._editparams)
        if editparams is not None:
            self.editparams.update_values(editparams, strict=True)
        else:
            logger.warning("No edit parameters provided, using default values.")

    def generate_data_reader(self):
        self.m = Multiread(self.files, sonar="wh", ibad=self.ibad)
        tsdat = self.m.read(varlist=["VariableLeader"])
        # Initial units: 10 Pa (about 1 mm or 0.001 decibar).
        # Converting to decibars.
        tsdat.pressure = (
            tsdat.VL["Pressure"] / 1000.0 * self.pressure_scale_factor
        )
        tsdat.temperature = tsdat.VL["Temperature"] / 100.0
        self.tsdat = tsdat

        self.yearbase = self.m.yearbase
        t0 = tsdat.dday[0]
        t1_adcp = self.driftparams.get("end_adcp", None)
        if t1_adcp is not None:
            t1_pc = to_day(self.m.yearbase, *self.driftparams["end_pc"])
            t1_adcp = to_day(self.m.yearbase, *self.driftparams["end_adcp"])
            self.rate = (t1_pc - t0) / (t1_adcp - t0)
        else:
            self.rate = 1
        self.t0 = t0
        self.dday = self.correct_dday(tsdat.dday)

        # We have to get the up/down reading from sysconfig for
        # a time when the instrument was in the water, so we
        # use the middle of the deployment.
        imid = len(tsdat.dday) // 2
        middle = self.m.read(
            varlist=["VariableLeader"], start=imid, stop=imid + 1
        )
        self.orientation = "up" if middle.sysconfig.up else "down"
        self.sysconfig = middle.sysconfig

    def correct_dday(self, orig_dday):
        return self.t0 + self.rate * (orig_dday - self.t0)

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
        burst_average = self.tgridparams.burst_average

        # Save whether we are averaging over bursts or not.
        self.burst_average = burst_average
        # Generate time stamps and stuff.
        if not burst_average:
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
            print(
                f"time between pings within burst: {burst_dt * 24 * 60 * 60:1.1f} s"
            )
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
        dat.dday = self.correct_dday(dat.dday_orig)
        dat.pressure = dat.VL["Pressure"] / 1000.0 * self.pressure_scale_factor
        sign = -1 if self.orientation == "up" else 1
        pdepth = seawater.depth2(dat.pressure, self.lat)
        dat.depth = pdepth[:, np.newaxis] + sign * dat.dep
        return dat

    @property
    def magdec(self):
        if self._magdec is None:
            if self.lat is None:
                logger.warning("No magnetic declination is available; using 0")
                self._magdec = 0
            else:
                n = len(self.start_ddays)
                dday_mid = self.start_ddays[n // 2]
                y, m, d = to_date(self.yearbase, dday_mid)[:3]
                output = Popen(
                    [
                        "magdec",
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
        return self._magdec

    @property
    def raw(self):
        if self._raw is None:
            self._raw = io.read_raw_rdi(self.files)
            self._raw.coords["bin"] = (("z"), np.arange(self._raw.z.size))
        return self._raw

    def to_enu(self, ens):
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

    def edit(self, ens):
        """
        Apply editing to xyze
        """
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

    def burst_average_depth(self, ens):
        # Average the depth vectors if doing burst-averages. Otherwise we just
        # return all depth vectors of this ensemble.
        if self.burst_average:
            depth_mean = ens.depth.mean(axis=0)
            depth = np.tile(depth_mean, (ens.dday.size, 1))
        else:
            depth = ens.depth
        return depth

    def regrid_enu(self, ens, method="linear"):
        """
        add enu_grid
        """
        shape = (ens.dday.size, self.dgrid.size, ens.enu.shape[-1])
        enu_grid = np.ma.zeros(shape)
        enu_grid[:] = np.ma.masked
        # Average the depth vectors if doing burst-averages. Otherwise we just
        # return all depth vectors of this ensemble.
        # TO DO: If calculating averages over regular time intervals, we need
        # to low pass filter the pressure time series prior to creating the
        # depth vectors.
        depth = self.burst_average_depth(ens)
        for i in range(ens.dday.size):
            enu_grid[i] = interp1(
                depth[i], ens.enu[i], self.dgrid, axis=0, method=method
            )
        ens.enu_grid = enu_grid

    def regrid_amp(self, ens, method="linear"):
        """
        add amp_grid (averaged over all 4 beams)
        """
        shape = (ens.dday.size, self.dgrid.size)
        amp_grid = np.ma.zeros(shape)
        amp_grid[:] = np.ma.masked
        depth = self.burst_average_depth(ens)
        for i in range(ens.dday.size):
            amp_grid[i] = interp1(
                depth[i],
                ens.amp[i].mean(axis=-1),
                self.dgrid,
                axis=0,
                method=method,
            )
        ens.amp_grid = amp_grid

    def average_ensembles(self, start=None, stop=None):
        nens_orig = len(self.start_ddays)
        indices_orig = np.arange(nens_orig)
        indices = indices_orig[start:stop]
        if start is None and stop is None:
            logger.info('Averaging all ensembles')
        else:
            logger.info(f'Averaging ensembles {indices[0]} to {indices[-1]}')
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

        # Midpoints of averaging intervals. Now calculating the time stamps for
        # the averages directly in MCM.make_start_days() where we deal with the
        # rest of the time vector. There we now also correctly calculate time
        # stamps for burst averages.
        # dday = self.start_ddays[start:stop] + self.tgridparams.dt_hours / 48.0
        dday = self.dday_mid[start:stop]

        for i, iens in enumerate(tqdm(indices)):
            ens = self.read_ensemble(iens)
            if ens is not None:
                self.edit(ens)  # modifies xyze
                self.to_enu(ens)  # transform to earth coords (east, north, up)
                self.regrid_enu(ens)
                self.regrid_amp(ens)

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

        self._log_processing_params()

    def _safely_add_attribute_from_params(self, key, d):
        """Add attribute from dictionary and throw an exception if the key is missing.

        Parameters
        ----------
        key : str
        d : dict
        """
        try:
            setattr(self, key, d.pop(key))
        except KeyError:
            print(f"{key} is missing in input parameters")

    def _set_up_logger(self):
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
        datestr = gv.time.now_datestr()
        if has_mooring_id:
            logger.info(f"Processing {mooring} SN {sn} on {datestr}")
        else:
            logger.warning(
                "No meta data provided, logging to generic filename 'adcp_proc.log'"
            )
            logger.info(f"Processing on {datestr}")

    def _log_processing_params(self):
        logger.info("processing settings")
        logger.info("-------------------")
        _log_params(self.dgridparams)
        _log_params(self.tgridparams)
        _log_params(self.editparams)
        logger.info("-------------------")

    def _ave2nc(self):
        """Convert data structure from ave to xarray / netcdf file format."""
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
        ds.u.attrs = dict(long_name="u", units="m/s")
        ds.v.attrs = dict(long_name="v", units="m/s")
        ds.w.attrs = dict(long_name="w", units="m/s")
        ds.e.attrs = dict(long_name="error velocity", units="m/s")
        ds.z.attrs = dict(long_name="depth", units="m")
        ds.temperature.attrs = dict(long_name="temperature", units="°C")
        ds.pressure.attrs = dict(long_name="pressure", units="dbar")
        return ds

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
            gv.plot.axstyle(axi)

    def plot_pressure(self):
        fig, ax = gv.plot.quickfig(fgs=(6, 2.5))
        self.raw.pressure.plot(ax=ax, label="all")
        self.raw.pressure.where(self.raw.pressure > 50).plot(
            ax=ax, label="subsurface"
        )
        ax.invert_yaxis()
        ax.set(xlabel="", ylabel="pressure [dbar]")
        ax.legend()


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
