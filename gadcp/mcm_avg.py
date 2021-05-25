"""
Edit and average the data from moored ADCPs.

This was developed initially for the wh300 upward-looking
instruments at the top of the MIXET moorings.

To use it, one must provide information and a small script
like this (taken from MIXET):

#####################


from pycurrents.adcp.mcm_avg import MCM, Pingavg

outdir = "./output"
datadir = "../wh300_data"

fnamesdict = dict(M45N=["mcm_4_5N.dat"],
                  M15N=['_RDI_000_16470.000', '_RDI_512_16470.000'],
                  M05N=['_RDI_000_16678.000', '_RDI_512_16678.000'],
                 )

driftparamsdict = dict(M45N=dict(end_pc=(2012, 11, 16, 1, 21, 48),
                                 end_adcp=(2012, 11, 16, 1, 27, 12),
                                 start_dday=None,
                                 ),
                       M15N=dict(end_pc=(2012, 11, 14, 13, 54, 15),
                                 end_adcp=(2012, 11, 14, 13, 52, 54)),
                       M05N=dict(end_pc=(2012, 11, 13, 21, 00, 22),
                                 end_adcp=(2012, 11, 13, 21, 00, 19)),
                      )

#positions from Scott:
#M1    4.5N    04  30.09 N       156  00.70 E
#M2    3N       03  13.082 N     155  59.956 E
#M3    1.5N    01  30.128 N     156  00.424 E
#M4     .5N     00  33.847 N     156  00.009 E
#M5    .5S      00   34.495 S     156  00.06 E


positionsdict = dict(M45N=(156.0117, 4.5015),
                     M15N=(156.0071, 1.5021),
                     M05N=(156.0002, 0.5651),
                     )


editparams = dict(max_e=0.2,          # absolute max e
                  max_e_deviation=2,  # max in terms of sigma
                  )

dgridparams = dict(#dbot=int(self.p_median),
                   dtop=10,
                   d_interval=1,
                   )

tgridparams = dict(dt_hours = 1.0,
                   #t0 = t0,
                   #t1 = t1,
                   )


for key in fnamesdict.keys():
    mcm = MCM(fnamesdict[key], driftparamsdict[key], datadir=datadir)
    pa = Pingavg(mcm, lonlat=positionsdict[key],
                 dgridparams=dgridparams,
                 tgridparams=tgridparams,
                 )
    pa.average_ensembles()
    fname = key + '_hourly.npz'
    pa.save_npz(fname, outdir=outdir)

### end MIXET example
"""

import os
from subprocess import Popen, PIPE  # for magdec
import logging

import numpy as np

from pycurrents.system import Bunch
from pycurrents.num.nptools import rangeslice
from pycurrents.num import interp1
from pycurrents.adcp.rdiraw import Multiread
from pycurrents.codas import to_date, to_day
from pycurrents.adcp.transform import rdi_xyz_enu
from pycurrents.file import npzfile
from pycurrents.data import seawater

# Standard logging
L = logging.getLogger(__name__)


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
        """
        # self.fnames = [os.path.join(datadir, f) for f in fnames]
        self.fnames = fnames
        self.driftparams = driftparams
        self.lat = lat
        self.pressure_scale_factor = pressure_scale_factor

        self.m = Multiread(self.fnames, sonar, ibad=ibad)
        tsdat = self.m.read(varlist=["VariableLeader"])
        # initial units: 10 Pa (about 1 mm or 0.001 decibar)
        tsdat.pressure = tsdat.VL["Pressure"] / 1000.0 * pressure_scale_factor  # in decibars
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
        middle = self.m.read(
            varlist=["VariableLeader"], start=imid, stop=imid + 1
        )
        self.orientation = "up" if middle.sysconfig.up else "down"

    def correct_dday(self, orig_dday):
        return self.t0 + self.rate * (orig_dday - self.t0)

    def make_start_ddays(
        self, dday_start, dday_end, dt_hours, burst_average=False
    ):
        """
        Generate time stamps for time averaging in Pingavg.average_ensembles().

        Parameters
        ----------
        dday_start : float

        dday_end : float

        dt_hours :

        burst_average : bool, optional
            Average over bursts. Defaults to False.

        Notes
        -----
        If turning on burst averaging, other input values will be ignored.
        """
        if not burst_average:
            print("no burst average")
            self.dday_start = dday_start
            self.dday_end = dday_end
            self.dt = dt_hours / 24.0
            self.start_ddays = np.arange(dday_start, dday_end, dt_hours / 24.0)
        else:
            # determine burst length and time between bursts
            dday_diff = np.diff(self.dday)
            burst_dt = np.median(dday_diff)
            print(burst_dt * 24 * 60)
            # It seems safe to assume that the time between bursts is at least
            # four times as long as the time between individual bursts within a
            # burst.
            inter_burst_dt = np.median(dday_diff[dday_diff > burst_dt * 4])
            print(inter_burst_dt * 24 * 60)
            # find starting points of all bursts
            start_indices = np.flatnonzero(dday_diff > burst_dt * 4)
            start_indices += 1
            start_indices = np.insert(start_indices, 0, 0)
            self.start_ddays = self.dday[start_indices]
            self.dday_start = self.start_ddays[0]
            # generate a dt that is inclusive of one burst
            # we know the number of pings in a burst from the difference
            # between the start_indices:
            pings_per_burst = np.int32(np.median(np.diff(start_indices)))
            print(f"{pings_per_burst} pings per burst")
            print(f"{start_indices.shape[0]} bursts")
            self.dt = burst_dt * pings_per_burst + burst_dt * 3

            # keeping this in here for now to run a test, need to delete when
            # above fully works.
            # self.dday_start = dday_start
            self.dday_end = dday_end
            # self.dt = dt_hours / 24.0
            # self.start_ddays = np.arange(dday_start, dday_end, dt_hours / 24.0)

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
    Edit raw moored ADCP data, and average on a uniform time grid.
    """

    _editparams = dict(
        max_e=0.2,  # absolute max e
        max_e_deviation=2,  # max in terms of sigma
        min_correlation=64,  # 64 is RDI default
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
        default_dgridparams = dict(
            dbot=int(self.p_median), dtop=10, d_interval=1
        )
        self.dgridparams = Bunch(default_dgridparams)
        if dgridparams is not None:
            self.dgridparams.update_values(dgridparams, strict=True)
        self.dgrid = np.arange(
            self.dgridparams.dtop,
            self.dgridparams.dbot,
            self.dgridparams.d_interval,
            dtype=float,
        )

        p = mcm.tsdat.pressure
        at_depth = np.nonzero(p > self.p_median)[0][0]
        t0 = (2 + np.ceil(mcm.dday[at_depth] * 24)) / 24.0
        in_water = np.nonzero(p > self.p_median / 2)[0][-1]
        t1 = (np.floor(mcm.dday[in_water] * 24) - 2) / 24.0

        default_tgridparams = dict(
            dt_hours=0.5, t0=t0, t1=t1, burst_average=False
        )

        self.tgridparams = Bunch(default_tgridparams)
        if tgridparams is not None:
            self.tgridparams.update_values(tgridparams, strict=True)

        self.mcm.make_start_ddays(
            self.tgridparams.t0,
            self.tgridparams.t1,
            self.tgridparams.dt_hours,
            self.tgridparams.burst_average,
        )
        self.start_ddays = self.mcm.start_ddays

    # The following is modified from ladcp.py.
    @property
    def magdec(self):
        if self._magdec is None:
            if self.lonlat is None:
                L.info("No magnetic declination is available; using 0")
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
                L.info("magdec output is: %s", output)
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

    def regrid_enu(self, ens, method="linear"):
        """
        add enu_grid
        """
        shape = (ens.dday.size, self.dgrid.size, ens.enu.shape[-1])
        enu_grid = np.ma.zeros(shape)
        enu_grid[:] = np.ma.masked
        for i in range(ens.dday.size):
            enu_grid[i] = interp1(
                ens.depth[i], ens.enu[i], self.dgrid, axis=0, method=method
            )
        ens.enu_grid = enu_grid

    def regrid_amp(self, ens, method="linear"):
        """
        add amp_grid (averaged over all 4 beams)
        """
        shape = (ens.dday.size, self.dgrid.size)
        amp_grid = np.ma.zeros(shape)
        amp_grid[:] = np.ma.masked
        for i in range(ens.dday.size):
            amp_grid[i] = interp1(
                ens.depth[i],
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

        # midpoints of averaging intervals (dividing by 48 because of midpoint,
        # so half the value):
        dday = self.start_ddays[start:stop] + self.tgridparams.dt_hours / 48.0

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
