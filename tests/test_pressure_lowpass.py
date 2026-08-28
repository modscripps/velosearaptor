"""Tests for the pressure low-pass filter cutoff (issue #99).

`_lowpassfilter_pressure` hands `scipy.signal.filtfilt` a single sampling
frequency, so every ping is treated as `1/fs` apart. Deriving `fs` from the
*modal* ping interval means that on burst-sampled records the 1800 s cap is
expressed in a fictitious time base and can never bind.

All tests here are synthetic. The bundled raw files cannot separate this
effect from genuine pressure-sensor noise (see the notes in
`plans/qc-and-depth-offset.md`, item 9), so `_lowpassfilter_pressure` is driven
directly on a constructed time base.
"""

import warnings

import numpy as np
import pytest
from pycurrents.system import Bunch

from velosearaptor import madcp, tools

YEARBASE = 2021

# Synthetic burst schedule: 10 pings 1 s apart, once per hour.
PINGS_PER_BURST = 10
PING_DT = 1.0
BURST_PERIOD = 3600.0
N_BURSTS = 100

# Knockdown: 20 dbar for one hour, i.e. exactly one burst.
BASE_PRESSURE = 1000.0
KNOCKDOWN = 20.0
KNOCKDOWN_START = 50 * BURST_PERIOD
KNOCKDOWN_DURATION = 3600.0
NOISE_STD = 1.0


def uniform_seconds(n=1000, dt=10.0):
    """Seconds since record start for a uniformly sampled record."""
    return np.arange(n) * dt


def burst_seconds():
    """Seconds since record start for the synthetic burst schedule."""
    starts = np.arange(N_BURSTS) * BURST_PERIOD
    offsets = np.arange(PINGS_PER_BURST) * PING_DT
    return (starts[:, np.newaxis] + offsets[np.newaxis, :]).ravel()


def mean_dt(seconds):
    return (seconds[-1] - seconds[0]) / (seconds.size - 1)


def make_processor(seconds, pressure):
    """Build a bare `ProcessADCP` carrying only what the low pass needs.

    Bypasses `__init__` on purpose: reading a raw file would drag in a real,
    unhelpful time base.
    """
    p = madcp.ProcessADCP.__new__(madcp.ProcessADCP)
    p._pressure_scale_factor = 1.0
    p.tsdat = Bunch(
        yearbase=YEARBASE,
        dday=seconds / 86400.0,
        pressure=pressure,
    )
    return p


def record_filter_args(monkeypatch):
    """Capture the `(lowcut, fs)` handed to `tools.lowpassfilter`."""
    calls = []
    real = tools.lowpassfilter

    def spy(x, lowcut, fs, order=3):
        calls.append({"lowcut": lowcut, "fs": fs, "cutoff": 1.0 / lowcut})
        return real(x, lowcut, fs, order=order)

    monkeypatch.setattr(madcp.tools, "lowpassfilter", spy)
    return calls


def synthetic_pressure(seconds, seed=0, noise_std=NOISE_STD):
    """Base pressure plus a one-hour knockdown plus white noise."""
    knock = np.where(
        (seconds >= KNOCKDOWN_START) & (seconds < KNOCKDOWN_START + KNOCKDOWN_DURATION),
        KNOCKDOWN,
        0.0,
    )
    truth = BASE_PRESSURE + knock
    rng = np.random.default_rng(seed)
    return truth + rng.normal(0.0, noise_std, seconds.size), truth


def ignore_nonuniform_warning():
    """Context manager silencing the non-uniformity warning."""
    ctx = warnings.catch_warnings()
    ctx.__enter__()
    warnings.filterwarnings("ignore", message=".*not uniformly sampled.*")
    return ctx


class TestUniformRecord:
    """Regression guards: nothing about uniform records may change."""

    @pytest.mark.parametrize(
        ("dt", "expected_cutoff"),
        [
            (10.0, 500.0),  # 50 pings is shorter than 30 minutes
            (60.0, 1800.0),  # 30 minutes is shorter than 50 pings
        ],
    )
    def test_cutoff_unchanged(self, monkeypatch, dt, expected_cutoff):
        seconds = uniform_seconds(dt=dt)
        pressure, _ = synthetic_pressure(seconds)
        calls = record_filter_args(monkeypatch)

        make_processor(seconds, pressure)._lowpassfilter_pressure()

        assert len(calls) == 1
        assert calls[0]["fs"] == pytest.approx(1.0 / dt)
        assert calls[0]["cutoff"] == pytest.approx(expected_cutoff)

    def test_no_nonuniformity_warning(self):
        seconds = uniform_seconds()
        pressure, _ = synthetic_pressure(seconds)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            make_processor(seconds, pressure)._lowpassfilter_pressure()


class TestBurstRecord:
    """The burst schedule from issue #99."""

    def test_schedule_is_what_we_think_it_is(self):
        seconds = burst_seconds()
        diffs = np.diff(seconds)
        values, counts = np.unique(diffs, return_counts=True)
        modal = values[np.argmax(counts)]
        assert modal == pytest.approx(PING_DT)
        # The modal interval is off from the mean by more than two orders of
        # magnitude; that gap is the whole bug.
        assert mean_dt(seconds) == pytest.approx(356.77, abs=0.01)

    def test_cutoff_capped_in_wall_clock_time(self, monkeypatch):
        seconds = burst_seconds()
        pressure, _ = synthetic_pressure(seconds)
        calls = record_filter_args(monkeypatch)

        ctx = ignore_nonuniform_warning()
        try:
            make_processor(seconds, pressure)._lowpassfilter_pressure()
        finally:
            ctx.__exit__(None, None, None)

        assert len(calls) == 1
        # fs must describe the mean rate over the record, not the modal one.
        assert calls[0]["fs"] == pytest.approx(1.0 / mean_dt(seconds))
        # The cap is 30 minutes of wall-clock time and it must bind here:
        # 50 mean intervals span about 5 hours.
        assert calls[0]["cutoff"] == pytest.approx(1800.0)

    def test_knockdown_amplitude_recovered(self):
        seconds = burst_seconds()
        pressure, truth = synthetic_pressure(seconds)

        ctx = ignore_nonuniform_warning()
        try:
            lp = make_processor(seconds, pressure)._lowpassfilter_pressure()
        finally:
            ctx.__exit__(None, None, None)

        window = (seconds >= KNOCKDOWN_START) & (
            seconds < KNOCKDOWN_START + KNOCKDOWN_DURATION
        )
        peak = (lp[window] - BASE_PRESSURE).max()
        # Tolerance 3 dbar. A third-order Butterworth run through filtfilt
        # overshoots a one-hour boxcar by ~1.6 dbar with the cutoff at 1800 s
        # (measured with the noise switched off); 1 dbar of white noise adds up
        # to ~1 dbar more at the peak. Today's value is 8.2 dbar, so the guard
        # is nowhere near vacuous.
        assert peak == pytest.approx(KNOCKDOWN, abs=3.0)
        # And the filter must not be adding error relative to the input noise.
        assert np.sqrt(np.mean((lp - truth) ** 2)) < NOISE_STD

    def test_nonuniformity_warning_fires(self):
        seconds = burst_seconds()
        pressure, _ = synthetic_pressure(seconds)
        with pytest.warns(UserWarning, match="not uniformly sampled"):
            make_processor(seconds, pressure)._lowpassfilter_pressure()
