"""Percent good must not report "no data" as "every ping bad" (issue #82).

`burst_average_ensembles` computes `pg` on the instrument-relative bins and
then interpolates the result onto the universal depth grid. Depths outside the
instrument's profile come back from `interp1` masked, and the store

    pg[i] = pgi_grid.astype(np.int8)

turns those into 0 with a `RuntimeWarning: invalid value encountered in cast`.
A grid cell the instrument never sampled is then indistinguishable from one
where every ping was rejected by editing, and both read as "0 % good".

The published `ds.pg` happens to survive today only because `_ave2nc` masks
`pg` wherever `amp` is NaN, and in burst mode `_burst_average_depth` hands the
same tiled depth vector to the amplitude regridding and to this interpolation,
so the two sets of off-grid cells coincide exactly. That coincidence is not a
guarantee; `ave.pg` is a documented output in its own right and is wrong.
"""

import warnings

import numpy as np

from velosearaptor.madcp import ProcessADCP

META_DATA = {
    "mooring": "Test",
    "project": "Test",
    "lon": 0,
    "lat": 0,
}

# The only bundled file that reaches the burst path: 22 pings at 3 s, then an
# 837 s gap, 365 bursts.
BURST_FILE = "data/24606000.000"


def _burst_proc(rootdir, **kwargs):
    return ProcessADCP(
        rootdir / BURST_FILE,
        META_DATA,
        magdec=0.0,
        tgridparams={"burst_average": True},
        **kwargs,
    )


def test_burst_average_pg_is_nan_off_the_instrument_profile(rootdir):
    """Grid cells the burst never reached carry NaN, not 0 % good."""
    proc = _burst_proc(rootdir)
    proc.burst_average_ensembles()

    pg = proc.ave.pg
    # `amp` is regridded from the same tiled depth vector, so its mask is
    # exactly the set of grid cells outside the burst's profile.
    off_grid = np.ma.getmaskarray(proc.ave.amp)
    assert off_grid.sum() > 0

    assert np.issubdtype(np.asarray(pg).dtype, np.floating)
    assert np.all(np.isnan(np.asarray(pg)[off_grid]))


def test_burst_average_does_not_warn_casting_invalid_pg(rootdir):
    """The invalid-value cast fires once per burst on every run today."""
    proc = _burst_proc(rootdir)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        proc.burst_average_ensembles()

    casts = [
        w
        for w in caught
        if issubclass(w.category, RuntimeWarning)
        and "invalid value encountered in cast" in str(w.message)
    ]
    assert casts == []


def test_burst_average_pg_separates_no_data_from_every_ping_bad(rootdir):
    """A bin masked by `maskbins` reads 0 %; a cell with no data reads NaN.

    `maskbins` masks the velocities but leaves the amplitudes alone, so the
    masked bins land on the grid with finite `amp` and zero good pings. That
    is the one case where the two meanings of 0 can be told apart.
    """
    proc = _burst_proc(rootdir)
    proc.editparams.maskbins = proc.generate_binmask([4, 5])
    proc.burst_average_ensembles()

    pg = np.asarray(proc.ave.pg)
    off_grid = np.ma.getmaskarray(proc.ave.amp)

    every_ping_bad = pg == 0
    # Both populations have to be present or the test proves nothing.
    assert every_ping_bad.sum() > 0
    assert off_grid.sum() > 0

    # A cell is either "no data" or "all pings rejected", never both.
    assert not np.any(every_ping_bad & off_grid)
    assert np.all(np.isnan(pg[off_grid]))


def test_burst_average_published_pg_is_unchanged(rootdir):
    """The fix is internal: `ds.pg` keeps its dtype, its NaNs and its values."""
    proc = _burst_proc(rootdir)
    proc.burst_average_ensembles()
    ds = proc.ds

    assert ds.pg.dtype == np.float64

    published = ds.pg.values
    amp_missing = np.isnan(ds.amp.values)
    assert np.array_equal(np.isnan(published), amp_missing)

    finite = ~amp_missing
    assert published[finite].min() >= 0
    assert published[finite].max() <= 100

    # Still an integer-valued percentage, as in the other two averaging
    # methods. Without the `np.floor` that replaces the old `astype(np.int8)`,
    # 4806 of the 6570 finite cells on this file move by up to 1 point.
    assert np.array_equal(published[finite], np.floor(published[finite]))

    # Same numbers as the internal array, on the depth levels that survived
    # the `_ave2nc` drop.
    levels = np.isin(proc.dgrid, ds.depth.values)
    internal = np.asarray(proc.ave.pg)[:, levels].T
    assert np.array_equal(published[finite], internal[finite])


def test_burst_average_ngood_is_nan_off_the_instrument_profile(rootdir):
    """`ngood` conflates the same two cases as `pg`, via `np.nan_to_num`."""
    proc = _burst_proc(rootdir)
    proc.burst_average_ensembles()

    ngood = np.asarray(proc.ave.ngood)
    off_grid = np.ma.getmaskarray(proc.ave.amp)
    assert off_grid.sum() > 0

    assert np.issubdtype(ngood.dtype, np.floating)
    assert np.all(np.isnan(ngood[off_grid]))


def test_published_ngood_separates_no_data_from_zero_good_pings(rootdir):
    """Unlike `pg`, `ngood` reaches the file unmasked, so this one is visible.

    `_ave2nc` masks `pg` on `amp` and leaves `ngood` alone, so every off-grid
    cell was published as "0 good pings". On this file that is 365 of the 366
    zero cells; the remaining one is a real cell where no ping survived.
    """
    proc = _burst_proc(rootdir)
    proc.burst_average_ensembles()
    ds = proc.ds

    ngood = ds.ngood.values
    assert np.issubdtype(ngood.dtype, np.floating)

    no_data = np.isnan(ngood)
    amp_missing = np.isnan(ds.amp.values)
    assert np.array_equal(no_data, amp_missing)

    # The real zero-good cells survive as 0 and stay countable.
    zero_good = ngood == 0
    assert zero_good.sum() > 0
    assert not np.any(zero_good & no_data)

    # A count, so still integer-valued, and never more than the ping count.
    finite = ~no_data
    assert np.array_equal(ngood[finite], np.floor(ngood[finite]))
    npings = np.broadcast_to(ds.npings.values[np.newaxis, :], ngood.shape)
    assert np.all(ngood[finite] <= npings[finite])
