"""Tests for external pressure NaN handling.

`_generate_external_pressure_interpolator` patches NaN gaps at the beginning
and the end of an external pressure record with the on-deck (atmospheric)
pressure. These tests drive it with synthetic pressure records built on the
time base of a bundled raw file.
"""

import pathlib

import numpy as np
import pytest
import xarray as xr

from velosearaptor import io
from velosearaptor.madcp import ProcessADCP

DECK_START = 0.3  # dbar, on deck before deployment
DECK_END = 0.4  # dbar, on deck after recovery
IN_WATER = 1500.0  # dbar, deployed


@pytest.fixture(scope="module")
def adcp():
    """A ProcessADCP instance whose time base the synthetic pressure uses."""
    adcpfile = pathlib.Path(__file__).parent / "data/binmap_16670013.000"
    meta_data = {
        "mooring": "Test",
        "project": "Test",
        "lon": 0,
        "lat": 0,
    }
    # magdec is provided explicitly so the test does not shell out to the
    # optional magdec executable.
    return ProcessADCP(adcpfile, meta_data, magdec=0.0)


def _run(adcp, values):
    """Feed `values` in as external pressure and return the patched record.

    The synthetic pressure sits on the instrument time base, so the xarray
    interpolation inside `_generate_external_pressure_interpolator` is the
    identity and the NaN runs put into `values` arrive unchanged. Evaluating
    the resulting interpolator at `dat.dday` therefore returns exactly the
    array the NaN patching produced.
    """
    time = io.yday0_to_datetime64(adcp.tsdat.yearbase, adcp.tsdat.dday)
    assert values.shape == time.shape
    adcp._pressure_provided = xr.DataArray(
        values.copy(), coords={"time": time}, dims=["time"]
    )
    adcp._generate_external_pressure_interpolator(adcp.tsdat)
    return adcp._external_pressure_interpolator(adcp.tsdat.dday)


@pytest.fixture
def npings(adcp):
    return adcp.tsdat.dday.size


def test_leading_and_trailing_nan_runs_on_deck(adcp, npings):
    """Two NaN runs, both on deck, are both patched with the deck medians."""
    p = np.full(npings, IN_WATER)
    p[:25] = DECK_START
    p[-25:] = DECK_END
    p[:5] = np.nan
    p[-5:] = np.nan

    out = _run(adcp, p)

    assert np.all(np.isfinite(out))
    expected = p.copy()
    expected[:5] = DECK_START
    expected[-5:] = DECK_END
    assert np.array_equal(out, expected)


def test_interior_nan_run_raises(adcp, npings):
    """A NaN gap while deployed cannot be filled and must fail loudly."""
    p = np.full(npings, IN_WATER)
    p[:25] = DECK_START
    p[-25:] = DECK_END
    p[:5] = np.nan
    p[-5:] = np.nan
    p[900:950] = np.nan

    with pytest.raises(ValueError, match="interior"):
        _run(adcp, p)


def test_single_leading_run_only(adcp, npings):
    """One leading run: same result as before the multi-run fix."""
    p = np.full(npings, IN_WATER)
    p[:25] = DECK_START
    p[:5] = np.nan

    out = _run(adcp, p)

    expected = p.copy()
    expected[:5] = DECK_START
    assert np.array_equal(out, expected)


def test_trailing_run_does_not_overwrite_good_values(adcp, npings):
    """Good in-water values must survive the trailing patch.

    Fails if the trailing run is located by indexing `p_interpolated` with a
    value that indexes `i_nan` instead: with NaN runs at both ends, that
    slices from the end of the leading run and overwrites the entire deployed
    record with the on-deck median.
    """
    p = np.full(npings, IN_WATER)
    p[:2] = np.nan
    p[2:22] = DECK_START
    p[-22:-2] = DECK_END
    p[-2:] = np.nan

    out = _run(adcp, p)

    # The deployed part of the record is untouched.
    assert np.array_equal(out[22:-22], np.full(npings - 44, IN_WATER))
    assert np.count_nonzero(out == IN_WATER) == npings - 44
    expected = p.copy()
    expected[:2] = DECK_START
    expected[-2:] = DECK_END
    assert np.array_equal(out, expected)


def test_in_water_leading_run_raises(adcp, npings):
    """A leading NaN run recorded in the water is an error, not a deck gap."""
    p = np.full(npings, IN_WATER)
    p[-25:] = DECK_END
    p[:5] = np.nan

    with pytest.raises(ValueError, match="on deck"):
        _run(adcp, p)


def test_in_water_trailing_run_raises(adcp, npings):
    """Same for a trailing NaN run recorded in the water."""
    p = np.full(npings, IN_WATER)
    p[:25] = DECK_START
    p[-5:] = np.nan

    with pytest.raises(ValueError, match="on deck"):
        _run(adcp, p)


def test_no_nan_untouched(adcp, npings):
    """A record without NaN passes through unchanged."""
    p = np.full(npings, IN_WATER)
    p[:25] = DECK_START
    p[-25:] = DECK_END

    out = _run(adcp, p)

    assert np.array_equal(out, p)
