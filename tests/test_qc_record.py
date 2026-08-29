"""Tests for recording applied QC parameters in the output (issue #83)."""

import numpy as np
import pytest

from velosearaptor.madcp import ProcessADCP

META_DATA = {
    "mooring": "Test",
    "project": "Test",
    "lon": 0,
    "lat": 0,
}


@pytest.fixture
def proc(rootdir):
    adcpfile = rootdir / "data/binmap_16670013.000"
    return ProcessADCP(adcpfile, META_DATA, magdec=0.0)


def test_max_e_applied_recorded(proc):
    """The error-velocity threshold actually applied per averaging interval
    must be stored in the output dataset."""
    proc.average_ensembles()
    ds = proc.ds

    assert "max_e_applied" in ds
    assert ds.max_e_applied.dims == ("time",)
    vals = ds.max_e_applied.values
    assert np.isfinite(vals).any()
    # The applied threshold is min(max_e, sigma-based), so never above max_e.
    assert np.nanmax(vals) <= proc.editparams.max_e + 1e-9


def test_editparams_fully_recorded(proc):
    """maskbins must appear in the output attributes alongside the other
    editing parameters.

    `pg_limit` is recorded too, but only `burst_average_ensembles` applies
    it, so on this path the attribute reads "not applied" rather than the
    number. Recording the number here said the product had been screened at
    30% good when it had not (issue #102); see
    tests/test_provenance_attrs.py.
    """
    maskbins = np.zeros(proc.meta_data.NCells, dtype=bool)
    maskbins[3] = True
    proc.parse_editparams({"maskbins": maskbins, "pg_limit": 30})
    proc.average_ensembles()

    assert proc.ds.attrs["pg_limit"] == "not applied"
    assert list(proc.ds.attrs["maskbins"]) == [3]
