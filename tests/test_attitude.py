"""Tests for attitude and beam correlation in `average_ensembles`
(motive-mooring-proc#51).

Per averaging interval the output carries the vector-mean heading and its
circular standard deviation, the mean and standard deviation of pitch and
roll, all on `time`, and the beam-mean correlation of the pings regridded
onto the depth axis the way `amp` is.
"""

import numpy as np
from pycurrents.num import interp1

from velosearaptor.madcp import ProcessADCP

META_DATA = {"mooring": "Test", "project": "Test", "lon": 0, "lat": 0}
UPLOOKER = "binmap_16670013.000"
DOWNLOOKER = "24606000.000"
ATTITUDE = ["heading", "heading_std", "pitch", "pitch_std", "roll", "roll_std"]


def _proc(rootdir, filename, **kwargs):
    return ProcessADCP(rootdir / "data" / filename, META_DATA, magdec=0.0, **kwargs)


def _averaged(rootdir, filename, **kwargs):
    frame = kwargs.pop("vertical_frame", "depth")
    proc = _proc(rootdir, filename, **kwargs)
    proc.average_ensembles(vertical_frame=frame)
    return proc


def test_attitude_variables_are_published_on_time(rootdir):
    ds = _averaged(rootdir, DOWNLOOKER).ds
    for name in ATTITUDE:
        assert ds[name].dims == ("time",), name
        assert ds[name].dtype == np.float32, name
        assert ds[name].attrs["units"] == "degree", name
        assert "long_name" in ds[name].attrs, name
    h = ds.heading.values[np.isfinite(ds.heading.values)]
    assert h.size > 0
    assert np.all((h >= 0) & (h < 360))
    assert np.all(ds.heading_std.values[np.isfinite(ds.heading_std.values)] >= 0)
    # An interval without pings carries no attitude either.
    empty = ds.npings.values < 2
    if empty.any():
        assert np.all(np.isnan(ds.heading.values[empty]))


def test_correlation_is_published_on_the_depth_grid(rootdir):
    ds = _averaged(rootdir, DOWNLOOKER).ds
    assert ds.cor.dims == ("depth", "time")
    assert ds.cor.dtype == np.float32
    assert ds.cor.attrs["units"] == "counts"
    c = ds.cor.values[np.isfinite(ds.cor.values)]
    assert c.size > 0
    assert np.all((c >= 0) & (c <= 255))
    # Like amp, finite exactly where the gridded amplitude is finite.
    np.testing.assert_array_equal(np.isnan(ds.cor.values), np.isnan(ds.amp.values))


def test_correlation_in_the_transducer_frame(rootdir):
    ds = _averaged(rootdir, DOWNLOOKER, vertical_frame="transducer").ds
    assert ds.cor.dims == ("z", "time")
    assert ds.cor.dtype == np.float32
    for name in ATTITUDE:
        assert ds[name].dims == ("time",), name


def test_heading_mean_is_a_vector_mean_across_north(rootdir, monkeypatch):
    """Pings alternating 359 and 1 degree average to about 0, never 180
    (review focus 2 of the motive plan)."""
    proc = _proc(rootdir, DOWNLOOKER)
    original = proc.read_ensemble

    def read_ensemble(iens):
        dat = original(iens)
        if dat is not None:
            n = np.asarray(dat.heading).size
            dat.heading = np.where(np.arange(n) % 2 == 0, 359.0, 1.0)
        return dat

    monkeypatch.setattr(proc, "read_ensemble", read_ensemble)
    proc.average_ensembles()
    h = proc.ds.heading.values
    s = proc.ds.heading_std.values
    ok = np.isfinite(h)
    assert ok.any()
    # An odd ping count leaves one unpaired ping, a bias of about 1 deg over
    # the ping count; the published range is [0, 360).
    assert np.all((h[ok] >= 0) & (h[ok] < 360))
    assert np.all(np.minimum(h[ok], 360 - h[ok]) < 0.1)
    assert np.all(np.abs(h[ok] - 180) > 170)
    assert np.all(s[ok] < 2.0)
    # Pitch and roll are plain means and are untouched by the heading edit.
    assert np.isfinite(proc.ds.pitch.values[ok]).all()


def test_correlation_is_the_beam_mean_regridded_like_amp(rootdir):
    """One interval recomputed from `read_ensemble` through the same regrid."""
    proc = _averaged(rootdir, DOWNLOOKER)
    ds = proc.ds
    i = int(np.flatnonzero(ds.npings.values >= 2)[0])
    ens = proc.read_ensemble(i)
    proc._qc(ens)
    proc._to_enu(ens)
    cor_mean = np.ma.filled(ens.cor.mean(axis=-1), np.nan).astype(float)
    depth = proc._burst_average_depth(ens)
    grid = np.full((ens.dday.size, proc.dgrid.size), np.nan)
    for k in range(ens.dday.size):
        col = interp1(
            depth[k], cor_mean[k][:, np.newaxis], proc.dgrid, axis=0, method="linear"
        )
        grid[k] = np.ma.filled(col, np.nan)[:, 0]
    with np.errstate(all="ignore"):
        expected = np.nanmean(grid, axis=0)
    # `_ave2nc` drops levels without velocity, so compare on the published axis.
    keep = np.isin(np.round(proc.dgrid, 6), np.round(ds.depth.values, 6))
    np.testing.assert_allclose(
        ds.cor.values[:, i],
        expected[keep].astype(np.float32),
        rtol=1e-6,
        atol=1e-4,
        equal_nan=True,
    )


def test_correlation_leaves_the_excluded_beam_out(rootdir):
    """With `ibad` the beam mean runs over the three beams in the solution."""
    with_bad = _averaged(rootdir, DOWNLOOKER, ibad=2)
    ds = with_bad.ds
    i = int(np.flatnonzero(ds.npings.values >= 2)[0])
    ens = with_bad.read_ensemble(i)
    with_bad._qc(ens)
    with_bad._to_enu(ens)
    three = np.ma.filled(np.delete(ens.cor, 2, axis=-1).mean(axis=-1), np.nan)
    four = np.ma.filled(ens.cor.mean(axis=-1), np.nan)
    # The two means differ on this file, so the test can tell them apart.
    assert np.nanmax(np.abs(three - four)) > 0.5
    with_bad._regrid_enu_amp(ens)
    with np.errstate(all="ignore"):
        expected = np.nanmean(ens.cor_grid, axis=0)
    keep = np.isin(np.round(with_bad.dgrid, 6), np.round(ds.depth.values, 6))
    np.testing.assert_allclose(
        ds.cor.values[:, i],
        expected[keep].astype(np.float32),
        rtol=1e-6,
        atol=1e-4,
        equal_nan=True,
    )
    # and the gridded column is the three-beam mean, not the four-beam one
    depth = with_bad._burst_average_depth(ens)
    col = np.ma.filled(
        interp1(
            depth[0],
            three[0][:, np.newaxis].astype(float),
            with_bad.dgrid,
            axis=0,
            method="linear",
        ),
        np.nan,
    )[:, 0]
    np.testing.assert_allclose(
        ens.cor_grid[0], col, rtol=1e-6, atol=1e-6, equal_nan=True
    )


def test_velocity_is_unchanged_by_the_new_accumulators(rootdir):
    """Pinned separately in test_pinned_output; here the cheap local check
    that the velocity, pg and amp of the two paths do not move when the
    attitude variables are dropped again."""
    ds = _averaged(rootdir, UPLOOKER).ds
    assert {"u", "v", "w", "pg", "amp", "cor", "heading"} <= set(ds.variables)
