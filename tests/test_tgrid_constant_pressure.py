"""Time-grid defaults with an external pressure that is constant at depth.

An ADCP with a failed pressure sensor is processed with a constant pressure
at depth from mooring geometry. The default time window starts at the first
ping deeper than the median pressure. For a constant record no ping is
deeper than the median, and `parse_tgridparams` raised an IndexError even
when the caller passed `t0` and `t1`, which make the default unnecessary.
"""

import pathlib

import numpy as np
import pytest
import xarray as xr

from velosearaptor import io
from velosearaptor.madcp import ProcessADCP

ADCPFILE = pathlib.Path(__file__).parent / "data/binmap_16670013.000"
META_DATA = {"mooring": "Test", "project": "Test", "lon": 0, "lat": 0}
DECK = 0.3  # dbar
AT_DEPTH = 1174.3  # dbar


@pytest.fixture(scope="module")
def constant_pressure():
    """Deck pressure for the first and last 10% of pings, constant between."""
    a = ProcessADCP(ADCPFILE, META_DATA, magdec=0.0)
    time = io.yday0_to_datetime64(a.tsdat.yearbase, a.tsdat.dday)
    n = time.size
    p = np.full(n, AT_DEPTH)
    p[: n // 10] = DECK
    p[-(n // 10) :] = DECK
    return xr.DataArray(p, coords={"time": time}, dims=["time"]), a.dday


def test_explicit_window_with_constant_pressure(constant_pressure):
    p, dday = constant_pressure
    t0, t1 = float(dday[len(dday) // 4]), float(dday[len(dday) // 2])
    a = ProcessADCP(
        ADCPFILE,
        META_DATA,
        magdec=0.0,
        pressure=p,
        tgridparams={"t0": t0, "t1": t1, "dt_hours": 0.5},
    )
    assert a.tgridparams.t0 == pytest.approx(t0)
    assert a.tgridparams.t1 == pytest.approx(t1)


def test_default_window_with_constant_pressure(constant_pressure):
    """Without `t0`/`t1` the default window covers the pings at depth."""
    p, dday = constant_pressure
    a = ProcessADCP(ADCPFILE, META_DATA, magdec=0.0, pressure=p)
    n = len(dday)
    assert a.tgridparams.t0 == pytest.approx(dday[n // 10])
    assert a.tgridparams.t1 == pytest.approx(dday[n - n // 10 - 1])
