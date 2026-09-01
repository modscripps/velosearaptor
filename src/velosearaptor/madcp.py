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
import pathlib
from pathlib import Path
from shutil import which
from subprocess import PIPE, Popen  # for magdec
from typing import ClassVar
from warnings import warn

import gsw
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import xarray as xr
from pycurrents.adcp.rdiraw import FileBBWHOS, Multiread
from pycurrents.adcp.transform import Transform, rdi_xyz_enu
from pycurrents.codas import to_date, to_day
from pycurrents.data import seawater
from pycurrents.num import interp1
from pycurrents.num.nptools import rangeslice
from pycurrents.system import Bunch
from tqdm import tqdm

from . import __version__, io, tools

# Standard logging
logger = logging.getLogger(__name__)


def _find_magdec():
    """Locate the magdec executable.

    Looks on the system path first, then inside the package directory, then
    for a `geomag/magdec` built by `install_magdec.sh` in any directory above
    the package (the repository root of a source checkout).

    Returns
    -------
    str or None
        Path to the magdec executable, or None if it cannot be found.
    """
    magdec_path = which("magdec")
    if magdec_path is not None:
        return magdec_path

    package_dir = Path(__file__).resolve().parent

    candidate = package_dir / "magdec"
    if candidate.is_file():
        return str(candidate)

    # Walk up from the package directory looking for geomag/magdec. This finds
    # the repository root regardless of how deeply the package is nested, so it
    # keeps working across layout changes (e.g. the move to src/).
    for parent in package_dir.parents:
        candidate = parent / "geomag" / "magdec"
        if candidate.is_file():
            return str(candidate)

    return None


def _masked_bins(maskbins):
    """Bin indices `maskbins` masks, whichever form it was given in.

    `maskbins` is documented as a boolean array indexing into the bins (see
    `ProcessADCP.generate_binmask`), but `_edit` applies it as a plain numpy
    index, so an integer list of bin numbers masks exactly the same bins and
    is the natural thing to write in a parameter file. `np.flatnonzero` on
    such a list silently returns *positions* rather than bin numbers, and
    drops bin 0 because it is falsy. Dispatch on dtype instead.

    Returns an empty array for `None`, so callers can treat "no bins" once.
    """
    if maskbins is None:
        return np.array([], dtype=int)
    mb = np.asarray(maskbins)
    if mb.dtype == bool:
        return np.flatnonzero(mb)
    return mb.astype(int).ravel()


# The editing criteria, as flags (issue #30). Each returns a boolean array
# marking the cells it rejects, and none of them touches the data. Applying
# the flags is the caller's business, which is what lets `process_pings` run
# the two beam-space criteria per chunk and the error velocity threshold once
# over the whole record.
#
# Module level and explicitly parameterized, instead of methods reading
# `self`, so that a criterion can be evaluated without constructing a
# `ProcessADCP`. No stage in this module reads mutable pipeline state, so this
# costs nothing today and is what a composable pipeline would need (issue
# #15).


def _correlation_flag(cor, min_correlation):
    """Beams whose correlation is below `min_correlation`.

    The full beam axis is kept. Reducing it to a per-cell verdict is
    `_correlation_cell_flag`, which is where the decisions live that issues
    #18 and #30 want to change.
    """
    return np.asarray(cor) < min_correlation


def _no_data_flag(vel):
    """Beams carrying no velocity at a cell.

    Whichever cause put it there: the instrument rejected the beam, or
    binmapping had no source bin to fill the cell from (issue #78). The two
    are not separable, because `_binmap_all_beams` overwrites the beam data
    and some raw-masked cells come back finite from it.

    Editing does not cause this and does not mask for it; the cells are
    already masked when `_edit_masks` is entered. The flag exists so that
    validity can be assembled from flags alone, which is what percent good is
    counted from (issue #30).
    """
    return np.ma.getmaskarray(vel)


def _cell_flag(flag_beam, ibad=None):
    """Reduce a per-beam flag to a per-cell verdict.

    One failing beam rejects the whole cell, which is harsher than the
    instrument's own convention for the correlation threshold (issue #30).
    Beams excluded via `ibad` do not enter the velocity solution, so what they
    do must not reject cells either (issue #90); `Transform.beam_to_xyz` fills
    the excluded beam from the other three before transforming, so its mask
    never reaches `xyze` and neither may its flag.

    Shared by every beam-space criterion, so softening the whole-cell rule for
    three-beam solutions (issue #18) is one change, here.
    """
    kept = flag_beam if ibad is None else np.delete(flag_beam, ibad, axis=-1)
    return kept.any(axis=-1)


def _maskbins_flag(shape, maskbins):
    """Cells in the bins the user declared bad, for every ping."""
    flag = np.zeros(shape, dtype=bool)
    flag[:, _masked_bins(maskbins)] = True
    return flag


def _max_e_flag(e, max_e_applied):
    """Cells whose error velocity exceeds the applied threshold.

    Cells already masked in `e` come back False, so the flag carries only what
    this test rejects. They stay masked either way. A threshold that could not
    be computed, NaN for a fully masked ensemble, rejects nothing.
    """
    return np.ma.filled(np.abs(e) > max_e_applied, False)


# The editing criteria, in the order they are applied. That order is also the
# precedence that makes the published counts disjoint (issue #30).
QC_CRITERIA = ("nodata", "cor", "maskbins", "max_e")


def _attributed_flags(flag_no_data, flag_cor, flag_maskbins):
    """The beam-space criteria, made disjoint so each cell has one reason.

    The raw flags overlap, and the overlap is not a nuisance, it is a
    misattribution. `_correlation_flag` compares the fill data underneath
    already-masked cells, and `_mask_binmapped` writes 0 there, which is below
    any correlation threshold. So `flag_cor` fires on essentially every cell
    `flag_no_data` already claimed: 100% of them under binmapping, which is
    29% of everything `flag_cor` reports. Publishing the raw flags as
    per-criterion counts would report that as correlation failure (issue #30).

    Precedence is the order the criteria are applied. `flag_max_e` needs no
    exclusion: it is computed on an error velocity these three have already
    masked, and `_max_e_flag` fills masked cells with False.
    """
    return (
        flag_no_data,
        flag_cor & ~flag_no_data,
        flag_maskbins & ~(flag_no_data | flag_cor),
    )


class ProcessADCP:
    """Moored ADCP Processing.

    An instance of ProcessADCP is initialized by providing raw data and
    processing parameters. Once initialized, the instance method
    :meth:`average_ensembles` can be used to average over just a few or
    all pings.

    Magnetic declination is automatically calculated via a call to the command
    line tool
    [magdec](https://currents.soest.hawaii.edu/git/Oceanography_Tools/geomag) that
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
        Supply external pressure time series in dbar as xr.DataAarray with
        coordinate `time`.
    use_raw_pressure : bool, optional
        Pressure time series is low-pass filtered at 30 minutes or 50 ping
        cutoff period, whatever is shorter. Set this to True to use raw
        pressure. Defaults to False. Does not matter for burst-averaging
        (pressure is averaged over full burst) or when supplying external
        pressure time series.
    pressure_scale_factor : float, optional
        Can be used to scale a (likely broken) pressure time series. Defaults
        to 1 (no scaling).
    depth_offset : float, optional
        Constant offset in meters added to the published depth axis to correct
        a constant bias in the pressure sensor. Positive means the instrument
        was deeper than the pressure record says. Defaults to 0 (no offset).

        The offset is applied only to the output: to the `depth` coordinate of
        the averaging methods and to `xducer_depth` on all three paths. All
        gridding and interpolation happens in raw pressure-derived depth
        exactly as it does without an offset, so `u`, `v`, `w` and `pg` are
        bit-identical to an unoffset run and the published axis is the
        unoffset axis translated by a constant. This is the same result as
        shifting the finished file by hand, which supplying a corrected
        external pressure series cannot reproduce: `depth(p)` is nonlinear, so
        a constant shift in dbar is not a constant shift in meters, and
        re-deriving per-ping depths re-samples the interpolation onto the grid.

        Because nothing upstream of the output moves, `dtop` and `dbot` in
        `dgridparams` are read in the **raw**, uncorrected frame. Requesting
        `dtop=100` with `depth_offset=5` publishes an axis starting at 105.
        The offset does not apply to `z`, the transducer-relative distance
        published by `process_pings`, and it does not apply to the in-water
        pressure thresholds, which stay in pressure space.

        `pressure` always holds the unmodified measurement. When the offset is
        nonzero the equivalent corrected pressure is stored alongside it as
        `pressure_corrected`, so that recomputing depth from the file
        reproduces the published depth axis.

        Note that `pressure_scale_factor` is applied twice on the low-pass
        pressure path (a known bug), so combining the two parameters is not
        currently reliable at `pressure_scale_factor != 1`.


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
            depth grid in `burst_average_ensembles`. Not applied in
            `average_ensembles` as the user can filter based on pg later.

    """

    # Default editing parameters. Copied into a Bunch by `parse_editparams`,
    # so this is only ever read, never mutated.
    _editparams: ClassVar[dict] = {
        "max_e": 0.2,  # absolute max e
        "max_e_deviation": 2,  # max in terms of sigma
        "min_correlation": 64,  # 64 is RDI default
        "maskbins": None,  # do not mask any bins
        "pg_limit": 50,  # percent good limit applied in `burst_average_ensembles`
    }

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
        magdec=None,
        pressure=None,
        use_raw_pressure=False,
        pressure_scale_factor=1,
        depth_offset=0,
    ):
        self.meta_data = Bunch(meta_data.copy())
        self.ibad = ibad
        self.logdir = logdir
        self.verbose = verbose

        self._magdec_provided = magdec
        self._magdec = magdec

        self._pressure_provided = pressure
        self._pressure_scale_factor = pressure_scale_factor
        self._use_raw_pressure = use_raw_pressure
        self._depth_offset = float(depth_offset)

        self._raw = None
        self._default_dgridparams = None

        self.parse_file_locations(raw_data)
        self._initiate_data_reader()
        self._read_auxiliary_data()
        self._parse_meta_data()

        self._set_up_logger()

        self.parse_driftparams(driftparams)
        self._ensure_monotonic_dday()
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
            all_raw_files = sorted(dir.glob("*.00*"))
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
        # Auto-detect sonar type from the first file so that non-Workhorse
        # instruments (e.g. Sentinel V) are handled correctly.
        probe = FileBBWHOS(self.files[0], sonar=None, trim=False)
        sonar = str(probe.sonar)
        probe.close()
        self.m = Multiread(self.files, sonar=sonar, ibad=self.ibad)
        # Make some more meta data readily available by reading a single ping
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
        # the external pressure sensor. Patch a leading and/or a trailing NaN
        # run with atmospheric pressure, but only if the instrument was on
        # deck at the time. Any other NaN is an error: pressure sets the depth
        # of every velocity sample, so a gap that cannot be filled with
        # atmospheric pressure must fail loudly.
        i_good = np.flatnonzero(~np.isnan(p_interpolated))
        if i_good.size == 0:
            raise ValueError(
                "External pressure is NaN over the whole ADCP time range. "
                "Check that the external pressure record overlaps the ADCP "
                "record in time."
            )

        if i_good.size < p_interpolated.size:
            first_few_good_median = np.median(p_interpolated[i_good][:20])
            last_few_good_median = np.median(p_interpolated[i_good][-20:])

            # Lengths of the leading and the trailing NaN run (zero if the
            # record does not start / end with a NaN).
            n_leading = i_good[0]
            n_trailing = p_interpolated.size - 1 - i_good[-1]

            if n_leading > 0:
                if first_few_good_median < 1:
                    p_interpolated[:n_leading] = first_few_good_median
                else:
                    raise ValueError(
                        f"External pressure is NaN for the first {n_leading} "
                        "ping(s) of the ADCP record, but the first valid "
                        f"pressure ({first_few_good_median:.2f} dbar) shows "
                        "the instrument was not on deck (< 1 dbar). Such a "
                        "gap cannot be filled with atmospheric pressure; "
                        "provide external pressure covering the start of the "
                        "ADCP record."
                    )

            if n_trailing > 0:
                if last_few_good_median < 1:
                    p_interpolated[p_interpolated.size - n_trailing :] = (
                        last_few_good_median
                    )
                else:
                    raise ValueError(
                        f"External pressure is NaN for the last {n_trailing} "
                        "ping(s) of the ADCP record, but the last valid "
                        f"pressure ({last_few_good_median:.2f} dbar) shows "
                        "the instrument was not on deck (< 1 dbar). Such a "
                        "gap cannot be filled with atmospheric pressure; "
                        "provide external pressure covering the end of the "
                        "ADCP record."
                    )

            # Whatever NaN is left sits between good data.
            i_nan = np.flatnonzero(np.isnan(p_interpolated))
            if i_nan.size > 0:
                n_gaps = np.flatnonzero(np.diff(i_nan) > 1).size + 1
                raise ValueError(
                    f"External pressure has {n_gaps} interior NaN gap(s) "
                    f"covering {i_nan.size} ping(s) of the ADCP record "
                    f"(first at ping index {i_nan[0]}, last at {i_nan[-1]}). "
                    "Interior gaps cannot be filled with atmospheric "
                    "pressure; interpolate or otherwise fill the external "
                    "pressure record before passing it in."
                )

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

    @property
    def pressure_lp(self):
        """Low-pass filtered pressure time series"""
        if "pressure_lp" not in self.tsdat:
            self.tsdat.pressure_lp = self._lowpassfilter_pressure()
        return self.tsdat.pressure_lp

    def _lowpassfilter_pressure(self):
        t64 = io.yday0_to_datetime64(self.tsdat.yearbase, self.tsdat.dday)
        # `scipy.signal.filtfilt` has no notion of a time base; it treats every
        # sample as 1/fs apart. Derive fs from the *mean* sampling interval over
        # the record so the cutoff below is expressed in wall-clock time. Using
        # the modal interval instead puts a burst-sampled record on a fictitious
        # time base in which the 30 minute cap can never bind, and the filter
        # then removes a large part of a genuine mooring knockdown.
        dt = tools.timedelta64_to_s(np.diff(t64))
        mean_dt = tools.timedelta64_to_s(t64[-1] - t64[0]) / (t64.size - 1)
        # Same non-uniformity test as the burst detection in
        # :meth:`make_start_ddays`: a gap longer than four times the median ping
        # interval is a break between bursts.
        median_dt = np.median(dt)
        if np.any(dt > 4 * median_dt):
            warn(
                "Pressure record is not uniformly sampled (median ping interval "
                f"{median_dt:.1f} s, mean {mean_dt:.1f} s, longest gap "
                f"{dt.max():.1f} s). The low-pass filter treats all pings as "
                f"{mean_dt:.1f} s apart, so structure within a burst is not "
                "resolved. Pass use_raw_pressure=True to skip the filter.",
                stacklevel=2,
            )
        fs = 1 / mean_dt
        # Make cutoff period either 30 minutes or 50 data points,
        # whatever is shorter.
        cutoff = min(50 * mean_dt, 1800)
        if cutoff <= 2 * mean_dt:
            raise ValueError(
                f"Mean ping interval is {mean_dt:.1f} s, so a {cutoff:.0f} s "
                "cutoff period is at or below the Nyquist period of the record "
                "and cannot be filtered. Pass use_raw_pressure=True, or supply "
                "external pressure."
            )
        # `tsdat.pressure` already carries `pressure_scale_factor`, applied in
        # `_scale_pycurrents_pressure`. Do not apply it a second time here.
        return tools.lowpassfilter(self.tsdat.pressure, 1 / cutoff, fs)

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
                dtop = max(dtop, 10)
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

            self._default_dgridparams = {
                "dtop": dtop,
                "dbot": dbot,
                "d_interval": d_interval,
            }

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

        # `dtop` and `dbot` are in the uncorrected, pressure-derived frame, so
        # with a depth offset the published axis is not the one that was asked
        # for. Say so here rather than let it be discovered in the output.
        depth_offset = getattr(self, "_depth_offset", 0.0)
        if depth_offset:
            logger.info(
                "Depth grid spans %g to %g m in uncorrected depth; "
                "depth_offset of %g m puts the published axis at %g to %g m.",
                self.dgrid.min(),
                self.dgrid.max(),
                depth_offset,
                self.dgrid.min() + depth_offset,
                self.dgrid.max() + depth_offset,
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
            # in_water = np.nonzero(p > self.p_median / 2)[0][-1]
            # t1 = self.dday[in_water]
            at_depth_last = np.nonzero(p > self.p_median / 1.05)[0][-1]
            t1 = self.dday[at_depth_last]

        # Generate a set of default time gridding parameters and then update
        # from the input parameters provided.
        default_tgridparams = {
            "dt_hours": 0.5,
            "t0": t0,
            "t1": t1,
            "burst_average": False,
        }
        self.tgridparams = Bunch(default_tgridparams)
        if tgridparams is not None:
            # Convert time to dday if provided as str. The conversion must be
            # relative to the yearbase of the raw data - reading the year off
            # the requested time instead puts any window in a later calendar
            # year a whole year off.
            for time in ["t0", "t1"]:
                if time in tgridparams and isinstance(tgridparams[time], str):
                    t64 = np.datetime64(tgridparams[time])
                    tgridparams[time] = (
                        t64 - np.datetime64(f"{self.yearbase}-01-01")
                    ) / np.timedelta64(1, "D")
            # update parameters in processing object
            self.tgridparams.update_values(tgridparams, strict=True)
            data_t0, data_t1 = self.dday[0], self.dday[-1]
            req_t0, req_t1 = self.tgridparams.t0, self.tgridparams.t1
            # An empty or reversed window selects no data, and it passes the
            # overlap test below, which compares each endpoint against the data
            # range on its own.
            empty = req_t0 >= req_t1
            # A requested window that does not overlap the data at all yields
            # an empty record - fail loudly instead.
            no_overlap = req_t1 < data_t0 or req_t0 > data_t1
            if empty or no_overlap:
                times = io.yday0_to_datetime64(
                    self.yearbase,
                    [float(req_t0), float(req_t1), float(data_t0), float(data_t1)],
                ).astype("datetime64[s]")
                if empty:
                    raise ValueError(
                        f"Requested time range starts at {times[0]} and ends "
                        f"at {times[1]}, so it selects no data. `t1` must lie "
                        "after `t0`."
                    )
                raise ValueError(
                    f"Requested time range {times[0]} to {times[1]} does not "
                    f"overlap the data range {times[2]} to {times[3]}."
                )
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
        driftparams = {} if driftparams is None else driftparams
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

    def _uncorrect_dday(self, dday_corrected):
        """Invert :meth:`_correct_dday`.

        Used to carry a repair made on ``self.dday`` back into the raw
        ``self.tsdat.dday`` it was derived from.

        Parameters
        ----------
        dday_corrected : array-like
            Drift-corrected time vector.

        Returns
        -------
        array-like
            The corresponding raw time vector.

        """

        return self.t0 + (dday_corrected - self.t0) / self.time_drift_rate

    def _ensure_monotonic_dday(self):
        """Detect and fix non-monotonic time vectors.

        Backward and repeated time steps are located and each one is classified
        and repaired on its own, walking the record from start to end.  Four
        cases:

        - Repeated timestamps (a short run of identical values): the first
          occurrence is kept and the copies are spread linearly towards the
          next larger timestamp.  Neighbor averaging cannot repair these, so
          they are handled separately from backward steps.
        - Isolated non-monotonic pings (a short run of <= ``_max_interp_pings``
          pings below the pre-jump value): replaced by linear interpolation
          from the nearest valid neighbors.  The outlier may be the ping
          *before* the backward step rather than the ping after it — a forward
          clock spike — in which case the spike is the ping that gets moved and
          the good timestamps following it are left alone.  Whichever
          hypothesis moves fewer pings is the one applied.  The two are equally
          costly when exactly one ping sits below the jump, and that tie is
          resolved in favour of repairing the ping *after* the jump, which is
          what the pre-existing single-bad-ping behaviour did.  Array length is
          preserved so ``self.m`` index alignment is maintained.
        - Segment overlap (backward jump where all remaining pings stay below
          the pre-jump value): the array is truncated at the jump.  Only tail
          removal is performed so index alignment is maintained.
        - Ambiguous (longer non-monotonic stretch that eventually recovers, or
          a long run of repeated timestamps): raises ``ValueError`` so the
          caller can inspect and decide.

        The result is verified to be strictly increasing before returning; if
        the repair did not achieve that, ``ValueError`` is raised rather than
        passing a partly fabricated, still non-monotonic time vector on to the
        rest of the processing.

        Called automatically from ``__init__`` after ``parse_driftparams``
        sets ``self.dday``.
        """
        _max_interp_pings = 3

        dday = self.dday
        n = dday.size
        if n < 2 or np.all(np.diff(dday) > 0):
            return

        # `self.dday` is a drift-corrected copy of `self.tsdat.dday`, not the
        # same array, so every value repaired below has to be carried back into
        # the raw record at the end. Track which pings move rather than
        # rewriting the whole array, so pings the repair never touched stay
        # bit-identical.
        repaired = np.zeros(n, dtype=bool)

        positive = np.diff(dday)
        positive = positive[positive > 0]
        if positive.size == 0:
            raise ValueError(
                "Non-monotonic time: the time vector contains no forward step "
                "at all. Cannot auto-fix — inspect the time vector and handle "
                "manually."
            )
        median_dt = np.median(positive)

        def _fill_between(i0, i1, count):
            """Spread `count` pings linearly between dday[i0] and dday[i1]."""
            if i1 < n and dday[i1] > dday[i0]:
                step = (dday[i1] - dday[i0]) / (count + 1)
            else:
                step = median_dt
            dday[i0 + 1 : i0 + 1 + count] = dday[i0] + step * np.arange(1, count + 1)
            repaired[i0 + 1 : i0 + 1 + count] = True

        i = 1
        truncate_at = None
        while i < n:
            if dday[i] > dday[i - 1]:
                i += 1
                continue

            if dday[i] == dday[i - 1]:
                # --- repeated timestamps ---
                j = i
                while j < n and dday[j] == dday[i - 1]:
                    j += 1
                ndup = j - i
                if ndup > _max_interp_pings:
                    raise ValueError(
                        f"Non-monotonic time: {ndup + 1} pings share the "
                        f"timestamp at ping {i - 1}. Cannot auto-fix — inspect "
                        f"the time vector and handle manually."
                    )
                logger.warning(
                    "Non-monotonic time: spreading %d repeated timestamp(s) "
                    "at indices %s.",
                    ndup,
                    np.arange(i, j),
                )
                _fill_between(i - 1, j, ndup)
                i = j
                continue

            # --- strictly backward step at i ---
            k = i
            while k < n and dday[k] < dday[i - 1]:
                k += 1
            n_overlap = k - i

            # Which side of the jump is the outlier?  If dropping the single
            # ping before the jump restores order — the pings after it are
            # increasing and start above dday[i - 2] — then dday[i - 1] is a
            # forward clock spike.  Choose whichever hypothesis moves fewer
            # pings; on a tie the ping after the jump is the one repaired.
            spike_before = (i == 1 or dday[i] > dday[i - 2]) and np.all(
                np.diff(dday[i:k]) > 0
            )

            if spike_before and n_overlap > 1:
                logger.warning(
                    "Non-monotonic time: interpolating forward clock spike at "
                    "index %d.",
                    i - 1,
                )
                if i == 1:
                    dday[0] = dday[1] - median_dt
                    repaired[0] = True
                else:
                    dday[i - 1] = (dday[i - 2] + dday[i]) / 2
                    repaired[i - 1] = True
                i += 1
            elif n_overlap <= _max_interp_pings:
                logger.warning(
                    "Non-monotonic time: interpolating %d isolated "
                    "non-monotonic ping(s) at indices %s.",
                    n_overlap,
                    np.arange(i, k),
                )
                _fill_between(i - 1, k, n_overlap)
                i = k
            elif k >= n:
                # --- segment overlap: all remaining pings below pre-jump max ---
                logger.warning(
                    "Non-monotonic time: backward jump at ping %d with %d "
                    "overlapping pings. Truncating to first %d pings.",
                    i - 1,
                    n - i,
                    i,
                )
                truncate_at = i
                break
            else:
                # --- ambiguous: long non-monotonic stretch that recovers ---
                raise ValueError(
                    f"Non-monotonic time: backward jump at ping {i - 1} "
                    f"with {n_overlap} non-monotonic pings that eventually "
                    f"recover. Cannot auto-fix — inspect the time vector "
                    f"and handle manually."
                )

        if truncate_at is not None:
            keep = slice(None, truncate_at)
            self.dday = self.dday[keep]
            self.tsdat.dday = self.tsdat.dday[keep]
            self.tsdat.pressure = self.tsdat.pressure[keep]
            self.tsdat.temperature = self.tsdat.temperature[keep]
            self.tsdat.ens_num = self.tsdat.ens_num[keep]

        if not np.all(np.diff(self.dday) > 0):
            still_bad = np.flatnonzero(np.diff(self.dday) <= 0)
            raise ValueError(
                f"Non-monotonic time: the repair did not produce a strictly "
                f"increasing time vector; {still_bad.size} bad step(s) remain, "
                f"first at ping {still_bad[0]}. Inspect the time vector and "
                f"handle manually."
            )

        # Carry the repaired values back into the raw time base. Truncation
        # already sliced both arrays, so only indices still in range apply.
        moved = np.flatnonzero(repaired[: self.dday.size])
        if moved.size > 0:
            self.tsdat.dday[moved] = self._uncorrect_dday(self.dday[moved])

    def _parse_sysconfig(self):
        # We have to get the up/down reading from sysconfig for a time when the
        # instrument was in the water. Look at pressure and determine ensembles
        # during time at depth.
        ii = np.flatnonzero(self.tsdat.pressure > 15)
        if ii.size == 0:
            raise ValueError(
                "Could not determine ensemble index at depth for reading "
                "sysconfig: no pressure record exceeds 15 dbar."
            )
        depth_ii = (ii[0] + ii[-1]) // 2
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
            # Find starting points of all bursts.
            start_indices = np.flatnonzero(dday_diff > burst_dt * 4)
            if start_indices.size == 0:
                # No gap anywhere is longer than four ping intervals, so there
                # are no bursts to find. Without this guard the median below is
                # taken over an empty selection and returns NaN, and the ping
                # count per burst comes out as INT32_MIN, which then fails much
                # later with an unrelated message (issue #101).
                raise ValueError(
                    "This data is not burst sampled: no gap between pings is "
                    "longer than four times the median ping interval of "
                    f"{burst_dt * 24 * 60 * 60:1.1f} s. Set "
                    "tgridparams['burst_average'] to False and give a "
                    "dt_hours averaging interval instead."
                )
            inter_burst_dt = np.median(dday_diff[dday_diff > burst_dt * 4])
            print(f"time between bursts: {inter_burst_dt * 24 * 60:1.1f} min")
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
            raise ValueError(f"ens num {iens} is out of range")

        # Get indices within the interval
        sl = rangeslice(
            self.dday, self.start_ddays[iens], self.start_ddays[iens] + self.dt
        )
        # Use the indices to extract data
        dat = self.m.read(start=sl.start, stop=sl.stop)
        if dat is None:
            return None
        dat.dday_orig = dat.dday
        dat.dday = self._correct_dday(dat.dday_orig)
        if self._pressure_provided is not None:
            # Replace pressure if provided
            dat.pressure = self._external_pressure_to_dat(dat)
        elif self.is_burst_average or self._use_raw_pressure:
            # Use raw pressure
            dat.pressure = self._scale_pycurrents_pressure(dat)
        else:
            # Use low-pass filtered pressure
            mask = (self.tsdat.ens_num >= dat.ens_num[0]) & (
                self.tsdat.ens_num <= dat.ens_num[-1]
            )
            dat.pressure = self.pressure_lp[mask]

        sign = -1 if self.orientation == "up" else 1
        pdepth = seawater.depth2(dat.pressure, self.lat)
        dat.depth = pdepth[:, np.newaxis] + sign * dat.dep
        return dat

    @property
    def magdec(self):
        """Magnetic declination.

        If not provided as input argument magdec is calculated using
        [magdec](https://currents.soest.hawaii.edu/git/Oceanography_Tools/geomag)
        (must be installed) based on `lon` and `lat`.

        """
        if self._magdec is None:
            if self.lat is None:
                logger.warning(
                    "No lon/lat provided, cannot calculate magnetic declination."
                )
                self._magdec = 0
            else:
                magdec_path = _find_magdec()

                if magdec_path is None:
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

    def _edit_masks(self, ens):
        """Mask cells rejected by correlation and by `maskbins`.

        Split out of `_edit` so that `process_pings` can run this per chunk
        and still compute the error velocity threshold once over the whole
        record (issue #100).

        Both criteria are beam-space quantities, so both flags are computed
        before anything is written to `xyze`. `flag_cor_beam`, `flag_cor` and
        `flag_maskbins` are left on `ens` for the caller (issue #30).

        `flag_no_data_beam` and its cell reduction `flag_no_data` record the
        invalidity that arrives here rather than being caused here, and
        nothing is masked for them. Without it the flags do not add up to the
        mask on `xyze`: on the binmapped per-ping path nearly a third of every
        invalid cell was already invalid before any criterion looked at it.

        `ens.valid` is the complement of the flags raised so far, and is what
        percent good is counted from. `_edit` narrows it with `flag_max_e`;
        `process_pings` applies that flag itself, once over the whole record.

        `reason_nodata`, `reason_cor` and `reason_maskbins` are the same
        criteria made disjoint by precedence, so that each rejected cell is
        attributed to exactly one of them. They are what the published
        `nbad_*` counts are built from; `ens.valid` is unaffected, since the
        union is the same either way.
        """
        ep = self.editparams
        ens.flag_no_data_beam = _no_data_flag(ens.vel)
        ens.flag_no_data = _cell_flag(ens.flag_no_data_beam, self.ibad)
        ens.flag_cor_beam = _correlation_flag(ens.cor, ep.min_correlation)
        ens.flag_cor = _cell_flag(ens.flag_cor_beam, self.ibad)
        ens.flag_maskbins = _maskbins_flag(ens.flag_cor.shape, ep.maskbins)
        ens.valid = ~(ens.flag_no_data | ens.flag_cor | ens.flag_maskbins)
        (
            ens.reason_nodata,
            ens.reason_cor,
            ens.reason_maskbins,
        ) = _attributed_flags(ens.flag_no_data, ens.flag_cor, ens.flag_maskbins)

        ens.xyze[ens.flag_cor] = np.ma.masked
        # Before the standard deviation in `_edit`, not after: bins the user
        # has declared bad must not set the threshold that decides which of
        # the kept bins survive. When the masked bins are the noisy ones their
        # contribution pushes `max_e_deviation * sigma` past `max_e` and the
        # adaptive criterion switches off entirely (issue #100).
        ens.xyze[ens.flag_maskbins] = np.ma.masked

    def _adaptive_max_e(self, e):
        """The error velocity threshold to apply to the samples in `e`.

        `e` is whatever survived `_edit_masks`. Accumulate in float64: the
        per-ping path hands this a float32 array of the whole record, where a
        naive float32 sum of squares loses real precision.
        """
        ep = self.editparams
        max_e = min(ep.max_e, np.ma.std(e, dtype=np.float64) * ep.max_e_deviation)
        # Return a plain float (NaN if the threshold could not be computed,
        # e.g. for a fully masked ensemble) so it can be propagated to the
        # output dataset.
        return float(max_e) if np.isfinite(max_e) else np.nan

    def _edit(self, ens):
        """Apply editing to xyze.

        Leaves the flags of `_edit_masks` and `flag_max_e` on `ens`, and
        `ens.valid`, their complement. `ens.valid` is `~getmaskarray(xyze)`
        per cell once this returns, which is what makes percent good
        computable from the flags (issue #30).
        """
        self._edit_masks(ens)
        e = ens.xyze[:, :, 3]
        ens.max_e_applied = self._adaptive_max_e(e)
        ens.flag_max_e = _max_e_flag(e, ens.max_e_applied)
        ens.valid &= ~ens.flag_max_e
        # Already disjoint from the other three, so it is the flag unchanged.
        ens.reason_max_e = ens.flag_max_e
        ens.xyze[ens.flag_max_e] = np.ma.masked

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

    def _regrid_enu_amp(self, ens, method="linear"):
        """Depth-grid ENU velocities and amplitudes in a single pass.

        Combines the work of _regrid_enu and _regrid_amp, calling interp1
        once per ping on a concatenated (nbins, 6) array instead of twice
        on separate arrays. Output arrays use NaN instead of masked values.

        The sixth column is validity, masked exactly where the cell is
        invalid, so it rides along at no extra interpolation. `valid_grid` is
        where it came back finite, and `average_ensembles` counts percent good
        from it (issue #30). It is the only change this makes to a shared
        array, so the column must leave `enu_grid` and `amp_grid` bitwise
        alone, which `tests/test_qc_flags.py` pins.
        """
        npings = ens.dday.size
        ndgrid = self.dgrid.size
        ncols_enu = ens.enu.shape[-1]

        enu_grid = np.full((npings, ndgrid, ncols_enu), np.nan)
        amp_grid = np.full((npings, ndgrid), np.nan)
        valid_grid = np.zeros((npings, ndgrid), dtype=bool)

        depth = self._burst_average_depth(ens)

        for i in range(npings):
            amp_col = ens.amp[i].mean(axis=-1, keepdims=True)
            valid_col = np.ma.masked_array(
                np.ones((ens.valid.shape[1], 1)), mask=~ens.valid[i][:, np.newaxis]
            )
            combined = np.ma.concatenate([ens.enu[i], amp_col, valid_col], axis=-1)
            result = interp1(depth[i], combined, self.dgrid, axis=0, method=method)
            result_filled = np.ma.filled(result, np.nan)
            enu_grid[i] = result_filled[:, :ncols_enu]
            amp_grid[i] = result_filled[:, ncols_enu]
            valid_grid[i] = np.isfinite(result_filled[:, ncols_enu + 1])

        ens.enu_grid = enu_grid
        ens.amp_grid = amp_grid
        ens.valid_grid = valid_grid

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

    @staticmethod
    def _mask_binmapped(x):
        """Mask cells a binmapped field carries no data for (issue #78).

        `_binmap_one_beam` marks cells outside the mapped beam range with
        NaN. `ens.vel` is a masked array, so assigning NaN into it writes
        data without setting the mask; `ens.amp` and `ens.cor` are integer
        typed, so their NaN does not survive the assignment at all and is
        cast to whatever the platform makes of `uint8(nan)` (0 on x86).
        Convert the invalid cells into masked cells here, so that the mask,
        not a NaN or a cast artifact, is what travels into `xyze`.
        """
        bad = ~np.isfinite(np.ma.filled(x, np.nan))
        return np.ma.masked_array(np.where(bad, 0, np.ma.filled(x, 0)), mask=bad)

    def _binmap_all_beams(self, ens):
        """Binmap single ping data for all beams."""

        # `ens.amp` and `ens.cor` come off the instrument as plain integer
        # arrays. Promote all three fields to masked arrays so that the
        # cells binmapping cannot fill can be masked below instead of being
        # written as NaN (issue #78).
        ens.vel = np.ma.masked_array(ens.vel)
        ens.amp = np.ma.masked_array(ens.amp)
        ens.cor = np.ma.masked_array(ens.cor)

        for beam_number in [1, 2, 3, 4]:
            veli, ampi, cori = (
                self._mask_binmapped(x) for x in self._binmap_one_beam(ens, beam_number)
            )
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
            Index of the first ping to process. This is a **ping index** into
            the raw ping sequence, not a time and not an averaging interval.
            Defaults to None, i.e. the first ping at or after `tgridparams.t0`.
        stop : int, optional
            Index of the ping to stop before, again a **ping index**; the ping
            at `stop` is not processed, so `stop - start` pings are returned.
            Defaults to None, i.e. the first ping after `tgridparams.t1`.
        binmap : bool, optional
            Do binmapping of along-beam data.
        ens_size : int, optional
            Pings are processed in ensembles to reduce memory usage.
            This parameter sets how many pings are in an ensemble. The default is 50000.
            It does not affect the results: the error velocity threshold is
            computed once over the whole record, not per ensemble.

        Notes
        -----
        `start` and `stop` here do **not** mean what they mean in
        :meth:`average_ensembles` and :meth:`burst_average_ensembles`. Those
        index into the averaging intervals set up by :meth:`make_start_ddays`;
        this method does no averaging and has no intervals, so its `start` and
        `stop` can only index individual pings.

        """
        # Validate the explicit indices. Unchecked they fail far from here:
        # `start` past the record gives an empty `ens_idxs`, and `start > stop`
        # a negative `npings` that dies inside `np.ma.zeros`. That is the same
        # class of problem as the `UnboundLocalError` fixed below (issue #101).
        npings_total = self.dday.size
        for name, value in (("start", start), ("stop", stop)):
            if value is None:
                continue
            if value < 0:
                raise ValueError(
                    f"{name}={value} is negative. Ping indices are not wrapped "
                    "Python-style here; give an index between 0 and "
                    f"{npings_total}."
                )
            if value > npings_total:
                raise ValueError(
                    f"{name}={value} is past the end of the record, which has "
                    f"{npings_total} pings. Give an index between 0 and "
                    f"{npings_total}."
                )

        if start is None:
            idx_start = np.searchsorted(self.dday, self.dday_start)
        else:
            idx_start = start
        if stop is None:
            idx_stop = np.searchsorted(self.dday, self.dday_end)
        else:
            idx_stop = stop

        if idx_start >= idx_stop:
            raise ValueError(
                f"start={idx_start} and stop={idx_stop} select no pings. "
                "`stop` is exclusive, so it must lie after `start`. Ping "
                "indices count forward from the start of the record."
            )

        ens_idxs = np.hstack((np.arange(idx_start, idx_stop, ens_size), idx_stop))
        write_idxs = ens_idxs - ens_idxs[0]  # Arrays we write to start at index 0
        npings = idx_stop - idx_start
        nens = ens_idxs.size - 1
        ndgrid = self.tsdat.dep.size

        logger.info("Processing all pings")
        logger.info(f"Binmapping is {binmap}")

        # Kept for the output attributes: binmapping changes the velocities
        # and nothing in the file otherwise records whether it ran.
        self._binmap = binmap

        uvwe = np.ma.zeros((npings, ndgrid, 4), dtype=np.float32)

        # Percent good comes from here, not from what survives in `uvwe`.
        # False, not True: a chunk the reader cannot return skips the loop
        # body below without raising a single flag, and a True default would
        # publish it as 100 % good (issue #30).
        valid = np.zeros((npings, ndgrid), dtype=bool)

        # Percent good says how many pings survived. These say why the rest
        # did not, one reason per rejected cell (issue #30). `int8` because
        # there is one ping per cell on this path, so every count is 0 or 1,
        # and this is the largest product the package writes.
        nbad = {name: np.zeros((npings, ndgrid), dtype=np.int8) for name in QC_CRITERIA}

        pg = np.zeros((npings, ndgrid), dtype=np.int8)
        amp = np.ma.zeros((npings, ndgrid), dtype=np.float32)

        # temperature = np.ma.zeros((npings,), dtype=np.float32)
        # pressure = np.ma.zeros((npings,), dtype=np.float32)

        temperature = self.tsdat.temperature[idx_start:idx_stop]
        pressure = self.tsdat.pressure[idx_start:idx_stop]
        max_e_applied = np.full((npings,), np.nan, dtype=np.float32)

        dday = self.dday[idx_start:idx_stop]

        # Loop over ensembles
        for i in tqdm(range(nens)):
            ens = self.m.read(start=ens_idxs[i], stop=ens_idxs[i + 1])
            idx0 = write_idxs[i]
            idx1 = write_idxs[i + 1]

            if ens is not None:
                # I pulled this out of read_ensembles because we need to
                # calculate depth for _ave2nc to work.
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
                    self._calculate_xyze(ens, ibad=self.ibad)

                # Only the masks here. The error velocity threshold is
                # applied after the loop, over the whole record — see below.
                self._edit_masks(ens)  # modifies xyze
                self._to_enu(ens)  # transform to earth coords (east, north, up)

            else:
                uvwe[idx0:idx1] = np.ma.masked
                amp[idx0:idx1] = np.ma.masked
                # pressure[idx0:idx1] = np.ma.masked
                # temperature[idx0:idx1] = np.ma.masked
                continue

            uvwe[idx0:idx1] = ens.enu
            amp[idx0:idx1] = ens.amp.mean(axis=-1)  # Average over beams... why?
            valid[idx0:idx1] = ens.valid
            for name in ("nodata", "cor", "maskbins"):
                nbad[name][idx0:idx1] = ens[f"reason_{name}"]

        # Apply the error velocity threshold to the record as a whole. Running
        # `_edit` inside the loop made `ens_size` the window the standard
        # deviation was estimated over, so the same file processed with
        # different `ens_size` came out with different velocities even though
        # the parameter is documented as a memory knob (issue #100).
        #
        # No second pass and no extra array: `rdi_xyz_enu` carries the error
        # velocity through untouched, so `uvwe[..., 3]` holds exactly what
        # `_edit` would have seen in `xyze`, and masking a cell of `enu` in all
        # four components is equivalent to masking it before the rotation.
        # Both claims are pinned by tests/test_adaptive_max_e.py.
        e = uvwe[..., 3]
        max_e = self._adaptive_max_e(e)
        # The record-wide counterpart of `_edit`'s `flag_max_e`. The two
        # beam-space flags were raised per chunk inside the loop, so on this
        # path the three do not share a lifecycle (issue #30).
        flag_max_e = _max_e_flag(e, max_e)
        uvwe[flag_max_e] = np.ma.masked
        valid &= ~flag_max_e
        # The record-wide criterion, written after the loop where it is
        # raised. It is already disjoint from the other three.
        nbad["max_e"][:] = flag_max_e
        max_e_applied[:] = max_e

        # Percent good is binary for single-ping data: 100 where the ping
        # survived editing, 0 where it was edited out. Counted from the flags
        # the criteria raised rather than read back out of `uvwe`, so it does
        # not depend on the masking `_edit` performs as a side effect
        # (issue #30). `ndgrid` here is the instrument bin count, so nothing
        # is gridded between the flags and the count.
        pg[:] = 100 * valid.astype(np.int8)

        self.ave = Bunch(
            u=uvwe[..., 0],
            v=uvwe[..., 1],
            w=uvwe[..., 2],
            e=uvwe[..., 3],
            pg=pg,
            nbad_nodata=nbad["nodata"],
            nbad_cor=nbad["cor"],
            nbad_maskbins=nbad["maskbins"],
            nbad_max_e=nbad["max_e"],
            amp=amp,
            temperature=temperature,
            pressure=pressure,
            max_e_applied=max_e_applied,
            # npings=npings,
            dday=dday,
            yearbase=self.yearbase,
            dep=ens.dep,  # <<<<----- BAD!!! This a fudge because I don't want to calculated a depth vector. We use the depth vector from the last ensemble.
            editparams=self.editparams,
            tgridparams=self.tgridparams,
            # dgridparams=self.dgridparams,
            magdec=self.magdec,
            lon=self.lon,
            lat=self.lat,
        )

        self._processing_method = "process_pings"
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
        uvwe = np.full((nens, ndgrid, 4), np.nan, dtype=np.float32)
        uvwe_std = np.full((nens, ndgrid, 4), np.nan, dtype=np.float32)

        pg = np.zeros((nens, ndgrid), dtype=np.int8)
        ngood = np.zeros((nens, ndgrid), dtype=np.int32)
        amp = np.full((nens, ndgrid), np.nan, dtype=np.float32)

        temperature = np.full((nens,), np.nan, dtype=np.float32)
        pressure = np.full((nens,), np.nan, dtype=np.float32)
        pressure_std = np.full((nens,), np.nan, dtype=np.float32)
        pressure_max = np.full((nens,), np.nan, dtype=np.float32)

        npings = np.zeros((nens,), dtype=np.int32)
        max_e_applied = np.full((nens,), np.nan, dtype=np.float32)

        dday = self.dday_mid[start:stop]

        for i, iens in enumerate(tqdm(indices)):
            ens = self.read_ensemble(iens)
            if ens is not None:
                self._edit(ens)  # modifies xyze
                max_e_applied[i] = ens.max_e_applied
                self._to_enu(ens)  # transform to earth coords (east, north, up)
                self._regrid_enu_amp(ens)

                nprofs = ens.enu_grid.shape[0]
            else:
                nprofs = 0
            npings[i] = nprofs
            if nprofs < 2:
                continue

            with np.errstate(all="ignore"):
                uvwe[i] = np.nanmean(ens.enu_grid, axis=0)
                uvwe_std[i] = np.nanstd(ens.enu_grid, axis=0)
                amp[i] = np.nanmean(ens.amp_grid, axis=0)

            # Counted from the validity gridded alongside the velocities in
            # `_regrid_enu_amp`, not from what came back finite there
            # (issue #30). The two are the same array today; the point is that
            # the count no longer depends on `_edit` masking as a side effect.
            ngood[i] = np.sum(ens.valid_grid, axis=0)
            # Widen before multiplying: ngood is int32 storage, and
            # 100 * ngood wraps once ngood > 21474836, which would turn the
            # best bins negative (issue #94). ngood itself is int32 rather than
            # int16 so that averaging intervals longer than 32767 pings neither
            # raise nor wrap (issue #101).
            pgi = 100 * ngood[i].astype(np.int64) // nprofs
            pg[i] = pgi.astype(np.int8)

            p = np.ma.filled(ens.pressure, np.nan)
            t = np.ma.filled(ens.temperature, np.nan)
            with np.errstate(all="ignore"):
                pressure[i] = np.nanmean(p)
                pressure_std[i] = np.nanstd(p)
                pressure_max[i] = np.nanmax(p)
                temperature[i] = np.nanmean(t)

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
            ngood=ngood,
            amp=amp,
            temperature=temperature,
            pressure=pressure,
            pressure_std=pressure_std,
            pressure_max=pressure_max,
            npings=npings,
            max_e_applied=max_e_applied,
            dday=dday,
            yearbase=self.yearbase,
            dep=self.dgrid,
            editparams=self.editparams,
            tgridparams=self.tgridparams,
            dgridparams=self.dgridparams,
            magdec=self.magdec,
            lon=self.lon,
            lat=self.lat,
        )

        self._processing_method = "average_ensembles"
        self._ave2nc()
        self._add_meta_data_to_ds()
        self._log_processing_params()

    @staticmethod
    def _interpolation_neighbors(zi, nbins):
        """Bin indices used to interpolate over bin `zi`.

        Two bins on either side where they exist, clipped to the profile.
        Without the clipping, `zi < 2` indexes negative bins (which wrap to the
        far end of the profile and break the monotonicity `interp1` requires)
        and `zi > nbins - 3` runs off the end (issue #101).

        The clipping costs no accuracy. The interpolation is piecewise linear,
        so only the two bins bracketing `zi` enter the result and the outer
        pair never contributes; a clipped window such as `[0, 2, 3]` for
        `zi = 1` gives bitwise what the full window would.

        Parameters
        ----------
        zi : int
            Index of the bin to interpolate over.
        nbins : int
            Number of bins in the profile.

        Returns
        -------
        list of int
            Neighbouring bin indices, at least one on either side of `zi`.

        """
        if not 0 <= zi < nbins:
            raise ValueError(
                f"interpolate_bin={zi} is not a bin of this profile; it must "
                f"be between 0 and {nbins - 1}."
            )
        neighbors = [i for i in (zi - 2, zi - 1, zi + 1, zi + 2) if 0 <= i < nbins]
        if min(neighbors) > zi or max(neighbors) < zi:
            # The first and last bin have a neighbour on one side only, so the
            # target lies outside any window we could build. `interp1` masks an
            # out-of-range target rather than extrapolating to it, so clipping
            # here would leave the bin masked and do nothing at all. Refuse,
            # rather than silently ignore the request.
            raise ValueError(
                f"interpolate_bin={zi} has no neighbouring bin on both sides; "
                f"it must be between 1 and {nbins - 2} for this profile."
            )
        return neighbors

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
            Interpolate over a single, previously masked, bin. Index into the
            instrument-relative bins, so it must be between 1 and `nbins - 2`;
            the first and last bin have no neighbour on one side and raise.
            The interpolated bin keeps its own `pg` and `ngood`, which is zero
            for a masked bin, so the interpolation stays visible in the output.
            Defaults to None (no interpolation).

        """
        pg_condition = self.editparams.pg_limit
        # Kept for the output attributes: an interpolated bin is not measured
        # data, so which bin was filled in has to be recoverable from the file.
        self._interpolate_bin = interpolate_bin
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

        # `pg` and `ngood` are float and start at NaN so that a grid cell the
        # instrument never sampled stays distinguishable from one where every
        # ping was rejected, which reads 0. For `pg`, float64 keeps the
        # published dtype, which `_ave2nc` already promotes when it masks `pg`
        # on `amp`. `ngood` is masked by nothing, so its off-grid cells reached
        # the file as "0 good pings"; giving it NaN changes its published dtype
        # from int32 to float64 (issue #82).
        pg = np.full((nens, ndgrid), np.nan, dtype=np.float64)
        ngood = np.full((nens, ndgrid), np.nan, dtype=np.float64)

        # Same NaN initialization and the same reason as `pg` and `ngood`
        # above: an interval with fewer than two pings is skipped below
        # without writing, and a grid depth off the profile must not read
        # "0 pings rejected" (issue #82).
        nbad = {
            name: np.full((nens, ndgrid), np.nan, dtype=np.float64)
            for name in QC_CRITERIA
        }
        amp = np.ma.zeros((nens, ndgrid), dtype=np.float32)

        temperature = np.ma.zeros((nens,), dtype=np.float32)
        pressure = np.ma.zeros((nens,), dtype=np.float32)
        pressure_std = np.ma.zeros((nens,), dtype=np.float32)
        pressure_max = np.ma.zeros((nens,), dtype=np.float32)

        npings = np.zeros((nens,), dtype=np.int32)
        max_e_applied = np.full((nens,), np.nan, dtype=np.float32)

        dday = self.dday_mid[start:stop]

        for i, iens in enumerate(tqdm(indices)):
            ens = self.read_ensemble(iens)
            if ens is not None:
                self._edit(ens)  # modifies xyze
                max_e_applied[i] = ens.max_e_applied
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
                # (pg is not a masked array; it is NaN-filled already)
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

            # Counted from the flags rather than from what survived in
            # `ens.enu` (issue #30). Every ping of a burst sits on one tiled
            # depth vector, so the count is exact on the instrument bins and
            # only the count is interpolated onto the depth grid below.
            #
            # The sum of a boolean array is a platform int, as `count` was, so
            # `100 * ngood_inst` is wide enough not to wrap. Do not "tidy"
            # this to a narrow integer: at more than 327 good pings an int16
            # product overflows and pg goes negative, which is what issue #94
            # was in `average_ensembles`. The float storage arrays above are
            # only written after pg is computed.
            ngood_inst = np.sum(ens.valid, axis=0)
            pgi_inst = 100 * ngood_inst // nprofs

            if pg_condition is not None:
                pgi_index = pgi_inst < pg_condition
                uvwe_inst[pgi_index, :] = np.ma.masked
                uvwe_std_inst[pgi_index, :] = np.ma.masked

            if interpolate_bin is not None:
                zi = interpolate_bin
                # Note that pg and ngood are deliberately left alone for the
                # interpolated bin, so the interpolation stays visible in the
                # product and cannot be mistaken for measured data.
                neighbors = self._interpolation_neighbors(zi, depth.size)
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

            # One `interp1` call for every count, the way `_regrid_enu_amp`
            # does it for velocity and amplitude. `interp1` takes a 2-D
            # `y_old` and casts it to float64 internally.
            #
            # `np.floor` reproduces the truncation the old `.astype(np.int8)`
            # and `.astype(np.int32)` performed on the way in, so the
            # published values do not move and these stay integer-valued.
            # Drop it only as a deliberate change.
            #
            # Depths outside the instrument's profile come back masked. Fill
            # them with NaN rather than casting them to 0, which would claim
            # every ping at that depth was bad.
            counts_inst = np.column_stack(
                [pgi_inst, ngood_inst]
                + [np.sum(ens[f"reason_{name}"], axis=0) for name in QC_CRITERIA]
            )
            counts_grid = np.floor(
                np.ma.filled(
                    interp1(depth, counts_inst, self.dgrid, axis=0, method="linear"),
                    np.nan,
                )
            )
            pg[i] = counts_grid[:, 0]
            ngood[i] = counts_grid[:, 1]
            for j, name in enumerate(QC_CRITERIA):
                nbad[name][i] = counts_grid[:, 2 + j]

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
            ngood=ngood,
            nbad_nodata=nbad["nodata"],
            nbad_cor=nbad["cor"],
            nbad_maskbins=nbad["maskbins"],
            nbad_max_e=nbad["max_e"],
            amp=amp,
            temperature=temperature,
            pressure=pressure,
            pressure_std=pressure_std,
            pressure_max=pressure_max,
            npings=npings,
            max_e_applied=max_e_applied,
            dday=dday,
            yearbase=self.yearbase,
            dep=self.dgrid,
            editparams=self.editparams,
            tgridparams=self.tgridparams,
            dgridparams=self.dgridparams,
            magdec=self.magdec,
            lon=self.lon,
            lat=self.lat,
        )

        self._processing_method = "burst_average_ensembles"
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
        datestr = datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%d")

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
                logstr = (k, ":", _masked_bins(v))
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
            except AttributeError:
                # Scalars have no .shape; treat them as integer-valued entries.
                tmp = None
                varsint.append(ki)
            # print(ki, tmp)

        # generate time vector
        # Naive on purpose, see velosearaptor.io.yday0_to_datetime64.
        base = datetime.datetime(dat.yearbase, 1, 1, 0, 0, 0)  # noqa: DTZ001
        time = [base + datetime.timedelta(days=ti) for ti in dat.dday]
        adcptime = [np.datetime64(ti, "ns") for ti in time]
        # generate Dataset
        out = xr.Dataset(
            {"pg": (["depth", "time"], dat.pg.T)},
            coords={"time": (["time"], adcptime), "depth": (["depth"], dat.dep)},
        )
        for vari in vars2d:
            out[vari] = (["depth", "time"], dat[vari].T)
        for vari in vars1d:
            if vari not in ["dep", "dday"]:
                out[vari] = (["time"], dat[vari])

        # Percent good is currently defined everywhere. Set to NaN where we
        # don't have any amplitude data (i.e. no data).
        out["pg"] = out.pg.where(~np.isnan(out.amp), other=np.nan)

        # The counts say how many pings each criterion rejected. Off the
        # instrument's profile no ping was measured at all, so 0 would claim
        # something that was never true (issue #82).
        #
        # Not on `process_pings`, where `z` is the instrument bin axis: every
        # bin is sampled whenever the ping was read, the only NaN cell is an
        # unreadable chunk, and `u` and `pg` already report that. Masking
        # there would promote the counts from int8 to float64 and more than
        # double the largest product this package writes.
        if getattr(self, "_processing_method", None) != "process_pings":
            for name in QC_CRITERIA:
                var = f"nbad_{name}"
                if var in out:
                    out[var] = out[var].where(~np.isnan(out.amp), other=np.nan)

        # Drop depth levels that carry no velocity. Keying the drop on the
        # whole Dataset would keep levels alive on the strength of amp alone:
        # editing masks velocities only, so amp can be finite in levels with
        # zero finite velocity samples, and such levels should not define the
        # product's depth axis.
        has_velocity = np.isfinite(out.u).any(dim="time")
        out = out.isel(depth=has_velocity.values)

        # Drop pressure_std, pressure_max, and e_std
        dropvars = ["pressure_std", "pressure_max", "e_std"]
        for var in dropvars:
            if var in out:
                out = out.drop(var)

        # Change u/v/w std to standard error. Scale by the number of pings
        # that actually entered the mean in each bin; fall back to the total
        # ping count if the per-bin count is not available.
        if any(f"{var}_std" in out for var in ["u", "v", "w"]):
            if "ngood" in out:
                n_samples = out["ngood"].where(out["ngood"] > 0)
            else:
                n_samples = out["npings"]
            for var in ["u", "v", "w"]:
                if f"{var}_std" in out:
                    out = out.rename({f"{var}_std": f"{var}_error"})
                    out[f"{var}_error"] = out[f"{var}_error"] / np.sqrt(n_samples)

        # Calculate transducer depth from pressure
        out["xducer_depth"] = -gsw.z_from_p(out.pressure, self.lat)

        # Apply the constant depth offset (issue #92). This is the only place
        # it enters. Gridding and interpolation upstream run entirely in raw
        # pressure-derived depth, so translating the axis here leaves the
        # velocity field bit-identical to an unoffset run by construction,
        # which is what the issue asks for and what re-deriving depths from a
        # corrected pressure series cannot deliver.
        #
        # `getattr` because tests/test_dropna.py and tests/test_std_error.py
        # call this method on a `ProcessADCP.__new__` shell that never runs
        # `__init__`.
        depth_offset = getattr(self, "_depth_offset", 0.0)
        if depth_offset:
            out["xducer_depth"] = out["xducer_depth"] + depth_offset
            # `depth` is water depth and moves with the instrument. The
            # `process_pings` axis is a distance from the transducer, renamed
            # to `z` below, and an offset on the instrument's depth says
            # nothing about it.
            if getattr(self, "_processing_method", None) != "process_pings":
                out = out.assign_coords(depth=out.depth.values + depth_offset)
            # The measurement stays in `pressure`. This is what a sensor at
            # the corrected depth would have read, so that recomputing depth
            # from the file reproduces the published axis.
            out["pressure_corrected"] = (
                ["time"],
                gsw.p_from_z(-out.xducer_depth.values, self.lat),
            )

        # `process_pings` does not regrid, so its vertical axis is the
        # distance from the transducer to the center of each bin, not water
        # depth, and it runs upwards for an uplooker. Publish it as `z`, the
        # name the raw dataset already uses for the same quantity, so that
        # `depth` always means water depth. Bin depth stays recoverable per
        # ping as xducer_depth +/- z. `_add_names_and_units` runs after the
        # rename and therefore attaches the `z` attributes, not the `depth`
        # ones.
        if getattr(self, "_processing_method", None) == "process_pings":
            out = out.rename({"depth": "z"})

        # add variable names and units for plotting
        out = self._add_names_and_units(out)
        out = self._add_depth_offset_comments(out, depth_offset)
        out = self._add_pg_comment(out, getattr(self, "_processing_method", None))
        out = self._drop_absent_ancillary_variables(out)

        self.ds = out

    @staticmethod
    def _add_depth_offset_comments(ds, depth_offset):
        """Record the applied depth offset on the variables it moved.

        Injected after `_add_names_and_units`, which assigns the static
        `io.cf_conventions` entries wholesale and would otherwise overwrite
        this. Nothing is written when no offset was applied, so a file never
        carries a correction it did not receive.
        """
        if not depth_offset:
            return ds
        note = (
            f"A constant depth_offset of {depth_offset} m has been added to "
            "this variable to correct a constant bias in the pressure sensor. "
            "Positive means the instrument was deeper than the pressure "
            "record says. The offset was applied to the output only; all "
            "gridding used the uncorrected pressure-derived depth, so the "
            "velocities are unchanged by it and depth gridding parameters "
            "were interpreted in the uncorrected frame."
        )
        for name in ["depth", "xducer_depth"]:
            if name in ds.variables:
                attrs = dict(ds[name].attrs)
                existing = attrs.get("comment")
                attrs["comment"] = f"{existing} {note}" if existing else note
                ds[name].attrs = attrs
        return ds

    # What `pg` counts differs on every path, and the name it shares with the
    # instrument's own field means something else again (issue #82). The text
    # is injected here because `io.cf_conventions` returns one static dict and
    # `_add_names_and_units` assigns it wholesale, so it cannot be stored
    # per path.
    _PG_COMMENTS: ClassVar[dict] = {
        "process_pings": (
            "Fraction of pings with a valid velocity at this bin, computed "
            "per ping by `process_pings`. It is binary, 0 or 100: 100 means "
            "the ping survived editing at this bin, 0 means it was edited "
            "out. No time averaging has happened on this path, so this says "
            "nothing about a sample count. The two averaging paths have to "
            "choose whether to count before or after depth gridding, and "
            "they choose differently; with no averaging window there is "
            "nothing to choose here."
        ),
        "average_ensembles": (
            "Fraction of the pings in each averaging interval with a valid "
            "velocity at this grid depth, computed by `average_ensembles` "
            "after every ping was interpolated onto the universal depth "
            "grid. The count is taken after gridding because the transducer "
            "moves within an averaging interval, by several meters against "
            "a 4 m bin on a typical mooring, so the pings in one interval "
            "share no common bin grid to count on. "
            "`burst_average_ensembles` counts on instrument bins before "
            "gridding, because within a burst the transducer does not move. "
            "`interp1` widens an edited bin, so both grid intervals "
            "touching one come out invalid and this reads lower than the "
            "fraction of pings that passed editing. It also absorbs mooring "
            "knockdown: a grid depth the instrument did not reach for part "
            "of the interval counts those pings as bad. `ngood` carries the "
            "sample count this is derived from."
        ),
        "burst_average_ensembles": (
            "Fraction of the pings in each burst with a valid velocity at "
            "this bin, computed by `burst_average_ensembles` on "
            "instrument-relative bins and then interpolated onto the "
            "universal depth grid. The count is taken before gridding "
            "because every ping in a burst is assigned the burst mean "
            "transducer depth, so the pings share one bin grid and counting "
            "in bin space is exact. `average_ensembles` counts after "
            "gridding instead, because its intervals are long enough for "
            "the transducer to move and no common bin grid exists. "
            "It is NaN at depths outside the "
            "instrument's profile, which is why it is float. Velocities were "
            "screened against the `pg_limit` attribute. A bin filled in by "
            "`interpolate_bin` keeps its own value of 0, so an interpolated "
            "bin stays visible in the product. `ngood` carries the sample "
            "count this is derived from."
        ),
    }

    _PG_RDI_NOTE = (
        "This is a velosearaptor-defined quantity. It shares its name with "
        "RDI's own four PercentGood fields, which "
        "`velosearaptor.io.read_raw_rdi` reads and no processing step uses, "
        "and it means something different from them."
    )

    @classmethod
    def _add_pg_comment(cls, ds, method):
        """Say which of the three quantities called `pg` this one is."""
        if "pg" not in ds.variables or method not in cls._PG_COMMENTS:
            return ds
        attrs = dict(ds.pg.attrs)
        attrs["comment"] = f"{cls._PG_COMMENTS[method]} {cls._PG_RDI_NOTE}"
        ds["pg"].attrs = attrs
        return ds

    @staticmethod
    def _drop_absent_ancillary_variables(ds):
        """Drop `ancillary_variables` names this file does not carry.

        `io.cf_conventions` declares `npings` on twelve entries, and
        `process_pings` writes no `npings`, so nine variables pointed at
        something absent: a CF checker fails, and a reader looking for the
        sample count behind `pg` or the error estimates finds nothing
        (issue #102). Done generically rather than per variable, so a
        variable added later cannot reintroduce it.
        """
        for v in ds.variables:
            declared = ds[v].attrs.get("ancillary_variables")
            if not declared:
                continue
            kept = [name for name in declared.split() if name in ds.variables]
            attrs = dict(ds[v].attrs)
            if kept:
                attrs["ancillary_variables"] = " ".join(kept)
            else:
                del attrs["ancillary_variables"]
            ds[v].attrs = attrs
        return ds

    def _add_names_and_units(self, ds):
        """Add variable attributes based on CF conventions.

        Parameters
        ----------
        ds : xarray.Dataset

        Returns
        -------
        ds : xarray.Dataset

        """
        CF = io.cf_conventions()
        for v in ds.variables:
            if v in CF:
                ds[v].attrs = CF[v]
        return ds

    @property
    def _pressure_source(self):
        """Which pressure record set the depths in this run.

        The four sources give different `xducer_depth` and, on the averaging
        paths, a different depth grid, and the choice is not recoverable from
        the numbers in the file.
        """
        if self._pressure_provided is not None:
            return "external"
        # `process_pings` scales the raw record per ping and never reaches the
        # low pass, which lives in `read_ensemble`.
        if getattr(self, "_processing_method", None) == "process_pings":
            return "raw"
        if self.is_burst_average or self._use_raw_pressure:
            return "raw"
        return "lowpass"

    def _yday_to_isoformat(self, yday):
        """Year day as a calendar time string for the output attributes."""
        return str(io.yday0_to_datetime64(self.yearbase, float(yday)))

    def _add_meta_data_to_ds(self):
        # Add some more info.
        method = getattr(self, "_processing_method", None)
        self.ds.attrs["velosearaptor_version"] = __version__
        self.ds.attrs["processing_method"] = "unknown" if method is None else method
        self.ds.attrs["orientation"] = self.orientation
        self.ds.attrs["magdec"] = self.magdec
        for att in ["max_e", "max_e_deviation", "min_correlation"]:
            value = self.editparams[att]
            # netcdf attributes cannot hold None.
            self.ds.attrs[att] = "none" if value is None else value
        # `pg_limit` screens percent good in `burst_average_ensembles` alone,
        # as the class docstring says. Writing the number into every file told
        # a reader the product had been screened at that threshold when it had
        # not (issue #102).
        if method == "burst_average_ensembles":
            value = self.editparams["pg_limit"]
            self.ds.attrs["pg_limit"] = "none" if value is None else value
        else:
            self.ds.attrs["pg_limit"] = "not applied"
        masked = _masked_bins(self.editparams["maskbins"])
        # netcdf attributes cannot hold an empty array either.
        self.ds.attrs["maskbins"] = "none" if masked.size == 0 else masked
        # Recorded even when zero: a reader needs to know the depth frame,
        # and "no offset applied" is part of that.
        self.ds.attrs["depth_offset"] = getattr(self, "_depth_offset", 0.0)

        # Where the depths came from, and how they were scaled.
        self.ds.attrs["pressure_source"] = self._pressure_source
        self.ds.attrs["pressure_scale_factor"] = self._pressure_scale_factor

        # Clock correction. `driftparams` is not written to the log file at
        # all, so this was previously unrecoverable.
        applied = self.driftparams.get("end_adcp") is not None
        self.ds.attrs["clock_correction"] = "applied" if applied else "none"
        self.ds.attrs["time_drift_rate"] = self.time_drift_rate

        self.ds.attrs["ibad"] = "none" if self.ibad is None else self.ibad

        # The requested time window, in calendar time. It is held internally
        # as a year day against `yearbase`, which is not a form a reader of
        # the netCDF file can be expected to decode.
        self.ds.attrs["t0"] = self._yday_to_isoformat(self.tgridparams.t0)
        self.ds.attrs["t1"] = self._yday_to_isoformat(self.tgridparams.t1)

        # Parameters that only mean something on the path that ran. Writing
        # them everywhere is the same fault as the old `pg_limit`.
        if method in ["average_ensembles", "burst_average_ensembles"]:
            # `make_start_ddays` reads `dt_hours` only when not burst
            # averaging; a burst gets its interval from the ping pattern. So
            # record the interval that actually governed, and `dt_hours`
            # itself only where it did something.
            self.ds.attrs["averaging_interval_hours"] = self.dt * 24
            if method == "average_ensembles":
                self.ds.attrs["dt_hours"] = self.tgridparams.dt_hours
            for att in ["dtop", "dbot", "d_interval"]:
                self.ds.attrs[att] = self.dgridparams[att]
        if method == "process_pings":
            # netcdf attributes cannot hold a boolean either.
            self.ds.attrs["binmap"] = str(bool(getattr(self, "_binmap", False)))
        if method == "burst_average_ensembles":
            value = getattr(self, "_interpolate_bin", None)
            self.ds.attrs["interpolate_bin"] = "none" if value is None else value

        # Add meta data if provided.
        if self.meta_data is not None:
            for k, v in self.meta_data.items():
                self.ds.attrs[k] = v
        self.ds.attrs["proc time"] = np.datetime64("now").astype("str")

    def plot_echo_stats(self):
        """Plot beam statistics (correlation and amplitude) from raw ADCP data."""
        r = self.raw

        _, ax = plt.subplots(
            nrows=1,
            ncols=2,
            figsize=(5, r.bin.max().data * 0.15),
            constrained_layout=True,
            sharey=True,
        )
        r.cor.mean(dim="time").plot(
            hue="beam", y="bin", marker="o", markersize=5, linestyle="", ax=ax[0]
        )
        r.amp.mean(dim="time").plot(
            hue="beam",
            y="bin",
            marker="o",
            markersize=5,
            linestyle="",
            ax=ax[1],
            add_legend=False,
        )
        ax[0].invert_yaxis()
        ax[1].set(ylabel="")
        ax[0].yaxis.set_major_locator(mpl.ticker.MultipleLocator(base=2))
        for axi in ax:
            axi.grid(True)

    def plot_pressure(self):
        """Plot pressure time series and mark time at depth."""
        _, ax = plt.subplots(
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
    [`notebooks/parameters.yml`](https://github.com/modscripps/velosearaptor/tree/main/notebooks/parameters.yml)
    """

    def __init__(self, parameter_file, mooring, sn, **kwargs):
        p = io.parse_yaml_input(parameter_file, mooring, sn)
        super().__init__(
            p["data_dir"],
            meta_data=p["meta_data"],
            driftparams=p["driftparams"],
            tgridparams=p["tgridparams"],
            dgridparams=p["dgridparams"],
            editparams=p["editparams"],
            **kwargs,
        )
