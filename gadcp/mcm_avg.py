"""
Edit and average the data from moored ADCPs.
"""

import logging
import os
from subprocess import PIPE, Popen  # for magdec

import numpy as np
from pycurrents.adcp.rdiraw import Multiread
from pycurrents.adcp.transform import rdi_xyz_enu
from pycurrents.codas import to_date, to_day
from pycurrents.data import seawater
from pycurrents.file import npzfile
from pycurrents.num import interp1
from pycurrents.num.nptools import rangeslice
from pycurrents.system import Bunch

# for the xyz_to_enu hotfix
# from pycurrents.adcp._transform import _heading_rotate, _heading_rotate_m
# from pycurrents.adcp._transform import _hpr_rotate, _hpr_rotate_m
# from pycurrents.adcp.transform import _process_vel, _process_attitude

# Standard logging
logger = logging.getLogger(__name__)


class MCM:
    """
    Data reader for moored ADCP.

    Includes clock correction, ensemble time range generation,
    and ability to read an ensemble at a time. It uses
    Multiread, so it can handle multi-file datasets.

    The clock drift correction is based on the assumption that the
    clock was set correctly when the instrument was deployed, and
    has drifted linearly until data collection was stopped.

    """

    def __init__(
        self,
        fnames,
        driftparams,
        datadir="./",
        sonar="wh",
        ibad=None,
        lat=30,
        pressure_scale_factor=1,
    ):
        """
        Generate MCM object for moored ADCP data.

        Parameters
        ----------
        fnames : list
            List with paths to raw data files.
        driftparams : dict
            Time drift information.
        datadir : str
            Raw data directory.
        sonar : str, optional
            ADCP type. Defaults to 'wh'. See Multiread docs for more info.
        ibad : int, optional
            The index of a beam to be excluded from the beam_to_xyz
            calculation. This is the zero-based index. Defaults to None.
        lat : float, optional
            Latitude in degrees for calculating depth from pressure. Defaults
            to 30.
        pressure_scale_factor : float, optional
            Scale factor for pressure time series. Defaults to 1 (no scaling).
            This was introduced for a malfunctioning pressure sensor and should
            not be necessary in most cases.
        """
        self.fnames = fnames
        self.driftparams = driftparams
        self.lat = lat
        self.pressure_scale_factor = pressure_scale_factor

        self.m = Multiread(self.fnames, sonar, ibad=ibad)
        tsdat = self.m.read(varlist=["VariableLeader"])
        # Initial units: 10 Pa (about 1 mm or 0.001 decibar).
        # Converting to decibars.
        tsdat.pressure = tsdat.VL["Pressure"] / 1000.0 * pressure_scale_factor
        self.tsdat = tsdat

        self.yearbase = self.m.yearbase

        t0 = driftparams.get("start_dday", None)
        t1_adcp = driftparams.get("end_adcp", None)
        if t0 is None:
            t0 = tsdat.dday[0]
        if t1_adcp is not None:
            t1_pc = to_day(self.m.yearbase, *driftparams["end_pc"])
            t1_adcp = to_day(self.m.yearbase, *driftparams["end_adcp"])
            self.rate = (t1_pc - t0) / (t1_adcp - t0)
        else:
            self.rate = 1
        self.t0 = t0
        self.dday = self.correct_dday(tsdat.dday)

        # We have to get the up/down reading from sysconfig for
        # a time when the instrument was in the water, so we
        # use the middle of the deployment.
        imid = len(tsdat.dday) // 2
        middle = self.m.read(varlist=["VariableLeader"], start=imid, stop=imid + 1)
        self.orientation = "up" if middle.sysconfig.up else "down"

    def correct_dday(self, orig_dday):
        return self.t0 + self.rate * (orig_dday - self.t0)

    def make_start_ddays(self, dday_start, dday_end, dt_hours, burst_average=False):
        """
        Generate time stamps for time averaging in Pingavg.average_ensembles().

        Parameters
        ----------
        dday_start : float
            Start time stamp.
        dday_end : float
            End time stamp.
        dt_hours : float
            Averaging interval in hours.
        burst_average : bool, optional
            Average over bursts. Defaults to False. If True, other parameters will be ignored.

        Notes
        -----
        If turning on burst averaging, other input values will be ignored and
        the averaging interval will be determined from the burst sampling
        scheme apparent in the ping pattern.
        """
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
            print(f"processing {start_indices.shape[0]} bursts")
            self.dt = burst_dt * pings_per_burst + burst_dt * 3

            self.dday_start = self.start_ddays[0]
            self.dday_end = dday_end

            # Time stamps in the middle of the burst
            self.dday_mid = self.start_ddays + pings_per_burst * burst_dt / 2

    def read_ensemble(self, iens):
        if iens > len(self.start_ddays) - 1:
            raise ValueError("ens num %d is out of range" % iens)

        # get indices within dday
        sl = rangeslice(
            self.dday, self.start_ddays[iens], self.start_ddays[iens] + self.dt
        )
        # use the indices to extract data
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


class Pingavg:
    """
    Edit raw moored ADCP data, and average on a uniform time grid or average
    overs bursts.
    """

    _editparams = dict(
        max_e=0.2,  # absolute max e
        max_e_deviation=2,  # max in terms of sigma
        min_correlation=64,  # 64 is RDI default
        maskbins=None,  # do not mask any bins
    )

    def __init__(
        self,
        mcm,
        lonlat=None,
        editparams=None,
        tgridparams=None,
        dgridparams=None,
    ):

        self.mcm = mcm
        self.yearbase = mcm.yearbase
        self.lonlat = lonlat
        self._magdec = None

        self.editparams = Bunch(self._editparams)
        if editparams is not None:
            self.editparams.update_values(editparams, strict=True)

        self.p_median = np.median(mcm.tsdat.pressure)
        default_dgridparams = dict(dbot=int(self.p_median), dtop=10, d_interval=1)
        self.dgridparams = Bunch(default_dgridparams)
        if dgridparams is not None:
            self.dgridparams.update_values(dgridparams, strict=True)
        self.dgrid = np.arange(
            self.dgridparams.dtop,
            self.dgridparams.dbot,
            self.dgridparams.d_interval,
            dtype=float,
        )

        # Find time at depth.
        p = mcm.tsdat.pressure
        at_depth = np.nonzero(p > self.p_median)[0][0]
        t0 = (2 + np.ceil(mcm.dday[at_depth] * 24)) / 24.0
        in_water = np.nonzero(p > self.p_median / 2)[0][-1]
        t1 = (np.floor(mcm.dday[in_water] * 24) - 2) / 24.0

        # Generate a set of default time gridding parameters and then update
        # from the input parameters provided.
        default_tgridparams = dict(dt_hours=0.5, t0=t0, t1=t1, burst_average=False)
        self.tgridparams = Bunch(default_tgridparams)
        if tgridparams is not None:
            self.tgridparams.update_values(tgridparams, strict=True)

        # Generate a time vector.
        self.mcm.make_start_ddays(
            self.tgridparams.t0,
            self.tgridparams.t1,
            self.tgridparams.dt_hours,
            self.tgridparams.burst_average,
        )
        self.start_ddays = self.mcm.start_ddays
        self.dday_mid = self.mcm.dday_mid
        self.burst_average = self.mcm.burst_average

    # The following is modified from ladcp.py.
    @property
    def magdec(self):
        if self._magdec is None:
            if self.lonlat is None:
                logger.warning("No magnetic declination is available; using 0")
                self._magdec = 0
            else:
                lonlat = self.lonlat
                n = len(self.start_ddays)
                dday_mid = self.start_ddays[n // 2]
                y, m, d = to_date(self.yearbase, dday_mid)[:3]
                output = Popen(
                    [
                        "magdec",
                        str(lonlat[0]),
                        str(lonlat[1]),
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
            orientation=self.mcm.orientation,
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

        for i, iens in enumerate(indices):
            ens = self.mcm.read_ensemble(iens)
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

        self.ave = Bunch(  # uvwe=uvwe,
            # uvwe_std=uvwe_std,
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
            lon=self.lonlat[0],
            lat=self.lonlat[1],
        )

    def save_npz(self, fname, outdir="./"):
        fpath = os.path.join(outdir, fname)
        npzfile.savez(fpath, **self.ave)


def rdi_xyz_enu_tmp(vel, heading, pitch, roll, orientation="down", gimbal=False):
    """
    GV: I used this as a hotfix for a bug in the pycurrents package.
    The bug should be fixed now and this should not be needed anymore.
    I am leaving this in here for now but go back to using rdi_xyz_enu()
    above in Pingavg.
    Also leaving the imports at the top needed for this function.

    Transform a velocity or sequence from xyz to enu.

    vel is a sequence or array of entries, [u, v, w, ...],
    where u, v and w are optionally followed by e; that is, vel
    can have 3, or 4 entries, and only the first three will
    be modified on output.  vel may be 1, 2, or 3-D.

    heading is a single value in compass degrees, or a 1-D sequence
    that must match the first dimension of vel if vel is 2-D or 3-D.

    pitch and roll have the same constraints as heading.
    pitch is tilt1; roll is tilt2;

    gimbal = False (default) adjusts the raw pitch (tilt1) for the
    roll (tilt2), as appropriate for the internal sensor.

    The data ordering convention is that if vel is 3-D, the indices
    are time, depth, component; heading is assumed to be a
    constant or a time series, but not varying with depth.

    This is a wrapper around a cython function.
    """
    vel, velshape = _process_vel(vel, min_comps=3)
    heading = _process_attitude(heading, "heading", vel, velshape)
    pitch = _process_attitude(pitch, "pitch", vel, velshape)
    roll = _process_attitude(roll, "roll", vel, velshape)

    if (
        np.ma.is_masked(vel)
        or np.ma.is_masked(heading)
        or np.ma.is_masked(pitch)
        or np.ma.is_masked(roll)
    ):
        # If any input has masked values, use the fully masked function.
        velmask = np.ma.getmaskarray(vel).astype(np.int8)
        hmask = np.ma.mask_or(
            np.ma.getmaskarray(roll), np.ma.getmaskarray(pitch), shrink=False
        )
        hmask = np.ma.mask_or(hmask, np.ma.getmaskarray(heading), shrink=False)
        hmask = hmask.astype(np.int8)  # cython 11.2 can't handle bool
        velr, outmask = _hpr_rotate_m(
            vel, heading, pitch, roll, velmask, hmask, orientation, gimbal
        )
        velr = np.ma.array(velr, mask=outmask, copy=False)
    else:
        # If either input is a masked array with mask False,
        # strip off the mask, do the calculation,
        # and convert the output back into a masked array.
        maskout = np.ma.isMaskedArray(vel)
        if maskout:
            vel = vel.view(np.ndarray)
        if np.ma.isMaskedArray(heading):
            heading = heading.view(np.ndarray)  # or heading.data?
        if np.ma.isMaskedArray(pitch):
            pitch = pitch.view(np.ndarray)
        if np.ma.isMaskedArray(roll):
            roll = roll.view(np.ndarray)
        velr = _hpr_rotate(vel, heading, pitch, roll, orientation, gimbal)
        if maskout:
            velr = np.ma.asarray(velr)

    return velr.reshape(velshape)
