"""Tests for velocity standard-error scaling (issue #80).

The standard error of a time-averaged velocity must be computed from the
number of pings that actually entered the mean in each depth bin, not from
the total number of pings in the averaging interval.
"""

import numpy as np
from pycurrents.system import Bunch

from velosearaptor.madcp import ProcessADCP

META_DATA = {
    "mooring": "Test",
    "project": "Test",
    "lon": 0,
    "lat": 0,
}


def _bare_instance_with_ave():
    """Build a ProcessADCP shell with a synthetic `ave` for testing _ave2nc."""
    p = ProcessADCP.__new__(ProcessADCP)
    p.lat = 0.0

    ntime, ndep = 2, 3
    npings = np.array([100, 50], dtype=np.int16)
    # Good pings per (time, depth) bin; one bin has none.
    ngood = np.array([[100, 50, 25], [50, 10, 0]], dtype=np.int16)
    pg = (100 * ngood // npings[:, np.newaxis]).astype(np.int8)
    u_std = np.full((ntime, ndep), 0.1, dtype=np.float32)

    p.ave = Bunch(
        u=np.ones((ntime, ndep), dtype=np.float32),
        v=np.ones((ntime, ndep), dtype=np.float32),
        w=np.ones((ntime, ndep), dtype=np.float32),
        e=np.zeros((ntime, ndep), dtype=np.float32),
        u_std=u_std.copy(),
        v_std=u_std.copy(),
        w_std=u_std.copy(),
        pg=pg,
        ngood=ngood,
        amp=np.full((ntime, ndep), 100.0, dtype=np.float32),
        temperature=np.full((ntime,), 10.0, dtype=np.float32),
        pressure=np.full((ntime,), 1000.0, dtype=np.float32),
        npings=npings,
        dday=np.array([0.0, 1 / 24]),
        yearbase=2020,
        dep=np.array([100.0, 110.0, 120.0]),
    )
    return p, u_std, ngood


def test_std_error_uses_good_pings():
    """u_error must be u_std / sqrt(ngood), not u_std / sqrt(npings)."""
    p, u_std, ngood = _bare_instance_with_ave()
    p._ave2nc()

    expected = u_std.T / np.sqrt(ngood.T)  # ds is (depth, time)
    got = p.ds.u_error.values

    valid = ngood.T > 0
    assert np.allclose(got[valid], expected[valid])
    # No good pings, no error estimate.
    assert np.all(np.isnan(got[~valid]))


def test_average_ensembles_outputs_ngood(rootdir):
    """average_ensembles must store the per-bin good-ping count that pg is
    derived from."""
    adcpfile = rootdir / "data/binmap_16670013.000"
    proc = ProcessADCP(adcpfile, META_DATA, magdec=0.0)
    proc.average_ensembles()
    ds = proc.ds

    assert "ngood" in ds

    pg = ds.pg.values
    ngood = ds.ngood.values
    npings = ds.npings.values[np.newaxis, :]

    valid = np.isfinite(pg) & (npings > 0)
    assert valid.sum() > 0
    assert np.all(ngood[valid] <= np.broadcast_to(npings, pg.shape)[valid])
    # pg is defined as floor(100 * ngood / npings).
    expected_pg = 100 * ngood[valid] // np.broadcast_to(npings, pg.shape)[valid]
    assert np.array_equal(pg[valid].astype(np.int64), expected_pg.astype(np.int64))
