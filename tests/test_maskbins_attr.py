"""`maskbins` must be recorded faithfully whichever form it was given in."""

import numpy as np
import pytest
import xarray as xr

from velosearaptor.madcp import ProcessADCP

META_DATA = {"mooring": "Test", "project": "Test", "lon": 0, "lat": 0}


@pytest.fixture
def proc(rootdir):
    """A processor with just enough state for the attribute writer: `ds` is
    only created by the averaging paths, and none of them is needed here."""
    p = ProcessADCP(rootdir / "data/binmap_16670013.000", META_DATA, magdec=0.0)
    p.ds = xr.Dataset()
    return p


def test_boolean_maskbins_record_the_masked_bins(proc):
    """The documented form: a boolean array indexing into the bins."""
    mask = np.zeros(10, dtype=bool)
    mask[[2, 5]] = True
    proc.parse_editparams({"maskbins": mask})
    proc._add_meta_data_to_ds()
    np.testing.assert_array_equal(proc.ds.attrs["maskbins"], [2, 5])


def test_integer_maskbins_record_the_same_bins(proc):
    """An integer index list masks the same bins as the equivalent boolean
    array, so it must record the same attribute -- `np.flatnonzero` on it
    silently returns positions instead, and drops bin 0 entirely."""
    proc.parse_editparams({"maskbins": [2, 5]})
    proc._add_meta_data_to_ds()
    np.testing.assert_array_equal(proc.ds.attrs["maskbins"], [2, 5])


def test_integer_maskbins_including_bin_zero(proc):
    """Bin 0 is falsy, so the flatnonzero path drops it."""
    proc.parse_editparams({"maskbins": [0, 1]})
    proc._add_meta_data_to_ds()
    np.testing.assert_array_equal(proc.ds.attrs["maskbins"], [0, 1])


def test_empty_maskbins_record_as_none(proc):
    """netCDF attributes cannot hold an empty array."""
    proc.parse_editparams({"maskbins": []})
    proc._add_meta_data_to_ds()
    assert proc.ds.attrs["maskbins"] == "none"


def test_absent_maskbins_record_as_none(proc):
    proc.parse_editparams({})
    proc._add_meta_data_to_ds()
    assert proc.ds.attrs["maskbins"] == "none"
