#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module gadcp.adcp with general adcp functions"""


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
        if var in ["pressure", "XducerDepth"]:
            ax.set(ylabel=var)
            gv.plot.ydecrease(ax)
        gv.plot.concise_date(ax)
        return ax

    # plot both pressure time series
    def plot_pressure(self):
        fig, (ax1, ax2) = plt.subplots(
            nrows=2, ncols=1, figsize=(10, 4), sharex=True
        )
        self.ts_plot(self.up, "pressure", ax1)
        self.ts_plot(self.dn, "pressure", ax2)

    def save_temp_files(self):
        # bring to same start/end times if appropriate
        dnt0 = self.dn.time[0]
        upt0 = self.up.time[0]
        self.deltat0 = np.timedelta64(upt0 - dnt0, "s")
        print("difference in start time:", self.deltat0)
        # cut time series at the beginning if time deltat less than ten minutes.
        self.cut_criterion = np.timedelta64(600, "s")
        if np.abs(self.deltat0) < self.cut_criterion:
            print("bringing to same start time")
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
        self.deltat1 = np.timedelta64(dnt1 - upt1, "s")
        print("difference in end time:", self.deltat1)
        # cut time series at the beginning if time deltat less than ten minutes.
        if np.abs(self.deltat1) < self.cut_criterion:
            print("bringing to same end time")
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
        print("\nsaving temporary files to current working directory\n")
        extract_raw(
            self.dnraw, "wh", dnii0, dnii1, outfile="tmpdn.rdi", verbose=False
        )
        extract_raw(
            self.upraw, "wh", upii0, upii1, outfile="tmpup.rdi", verbose=False
        )
        print("reading temporary files\n")
        self.dn = gv.io.read_raw_rdi_uh("tmpdn.rdi", auxillary_only=True)
        self.up = gv.io.read_raw_rdi_uh("tmpup.rdi", auxillary_only=True)

    # find split points
    def _lpfilt(self, data, pts):
        # design the Buterworth filter
        N = 3  # Filter order
        Wn = 1 / (pts * 2)  # Cutoff frequency
        B, A = sp.signal.butter(N, Wn, output="ba")
        # apply the filter
        lpdata = sp.signal.filtfilt(B, A, data)
        return lpdata

    def _split_point_function(self, dnup, dx, mindz):
        if dnup == "dn":
            pressure = self.dn.pressure
        elif dnup == "up":
            pressure = self.up.pressure
        # low pass filter pressure time series
        lptmp = self._lpfilt(pressure, 60)
        # find relative maxima and minima within +/-dx points
        tmpmax = sp.signal.argrelmax(lptmp, order=dx)[0]
        tmpmin = sp.signal.argrelmin(lptmp, order=dx)[0]

        allmaxmin = np.concatenate((tmpmax, tmpmin))
        allmaxminind = np.concatenate(
            (np.ones_like(tmpmax), np.zeros_like(tmpmin))
        )

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
        ax.plot(lptmp, color="0.5")
        ax.plot(deep_indices, pressure[deep_indices], "ro")
        ax.plot(shallow_indices, pressure[shallow_indices], "bo")
        ax.set(title=dnup)
        gv.plot.ydecrease(ax=ax)

        return shallow_indices, deep_indices

    def find_split_points(self, dx=500, mindz=300):
        self.ind_shallow_dn, self.ind_deep_dn = self._split_point_function(
            "dn", dx=dx, mindz=mindz
        )
        self.ind_shallow_up, self.ind_deep_up = self._split_point_function(
            "up", dx=dx, mindz=mindz
        )

    def plot_indices(self):
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 4))
        ax.plot(self.ind_shallow_up - self.ind_shallow_dn, "o")

    def split_files(self, experiment, cruise, cast, outdir):
        for adcp, indices, dnup in zip(
            [self.dn, self.up],
            [self.ind_shallow_dn, self.ind_shallow_up],
            ["dn", "up"],
        ):
            for ni, (ind0, ind1) in enumerate(zip(indices[:-1], indices[1:])):
                t0 = np.datetime64(adcp.time[ind0])
                t1 = np.datetime64(adcp.time[ind1])
                i0 = int(ind0)
                i1 = int(ind1)
                tbeg = (
                    np.datetime64(t0, "s")
                    .astype("str")
                    .replace("T", "_")
                    .replace(":", "")
                )
                tend = (
                    np.datetime64(t1, "s")
                    .astype("str")
                    .replace("T", "_")
                    .replace(":", "")
                )
                filename = "{}_{}_cast{:03g}_yo{:03g}_{}_{}_{}.rdi".format(
                    experiment, cruise, cast, ni, dnup, tbeg, tend
                )
                outfile = os.path.join(outdir, filename)
                rawfile = "tmp{}.rdi".format(dnup)
                extract_raw(
                    rawfile, "wh", i0, i1, outfile=outfile, verbose=False
                )
