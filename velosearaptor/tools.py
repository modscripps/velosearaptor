#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module velosearaptor.tools with various helper tools."""

import numpy as np
import scipy


def lowpassfilter(x, lowcut, fs, order=3):
    """Low-pass filter a signal using a butterworth filter.

    Parameters
    ----------
    x : array-like
        Time series.

    lowcut : float
        Cut-off frequency in units of fs.

    fs : float
        Sampling frequency.

    order : int
        Filter order.

    Returns
    -------
    lpx : array-like
        Low-pass filtered time series.

    Notes
    -----
    For example, if sampling four times per hour, fs=4. A cut-off period of 24
    hours is then expressed as lowcut=1/24.
    """
    b, a = _butter_lowpass(lowcut, fs, order=order)
    lpx = scipy.signal.filtfilt(b, a, x)
    return lpx


def _butter_lowpass(lowcut, fs, order=3):
    nyq = 0.5 * fs
    low = lowcut / nyq
    b, a = scipy.signal.butter(order, low, btype="lowpass")
    return b, a


def timedelta64_to_s(td64):
    """Convert np.timedelta64 with any time unit to float [s].

    Parameters
    ----------
    td64 : np.timedelta64
        Time delta.

    Returns
    -------
    np.float64
        Time delta in seconds.
    """
    assert td64.dtype.str[:3] == "<m8"
    return td64.astype("<m8[ns]").astype("float") / 1e9


def dominant_period_in_s(dt64):
    """Return the dominant period [s] in a time vector.

    Parameters
    ----------
    dt64 : np.datetime64
        Time vector.

    Returns
    -------
    period : np.float64
        Dominant period in seconds.
    """

    res = scipy.stats.mode(np.diff(dt64))
    return timedelta64_to_s(res.mode)
