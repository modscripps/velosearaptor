"""Tests for the transducer-relative vertical frame (issue #129).

`average_ensembles` and `burst_average_ensembles` interpolate every ping onto
`self.dgrid` before averaging. `vertical_frame="transducer"` keeps the
instrument's own bin axis instead and publishes it as `z`, the axis
`process_pings` already publishes, with the bin depth of each averaging
interval alongside it as a derived coordinate.
"""

import numpy as np
import pytest

from velosearaptor.madcp import ProcessADCP

META_DATA = {
    "mooring": "Test",
    "project": "Test",
    "lon": 0,
    "lat": 0,
}

UPLOOKER = "binmap_16670013.000"
DOWNLOOKER = "24606000.000"


def _proc(rootdir, filename, **kwargs):
    return ProcessADCP(rootdir / "data" / filename, META_DATA, magdec=0.0, **kwargs)


def _burst(rootdir, **kwargs):
    return _proc(rootdir, DOWNLOOKER, tgridparams={"burst_average": True}, **kwargs)


# --- the frame switch ---------------------------------------------------


@pytest.mark.parametrize("frame", ["depth", "transducer"])
def test_frame_is_recorded(rootdir, frame):
    """A file states which vertical frame it is in."""
    proc = _proc(rootdir, UPLOOKER)
    proc.average_ensembles(vertical_frame=frame)
    assert proc.ds.attrs["vertical_frame"] == frame


def test_process_pings_records_the_transducer_frame(rootdir):
    """`process_pings` has always published this frame; it now says so."""
    proc = _proc(rootdir, UPLOOKER)
    proc.process_pings()
    assert proc.ds.attrs["vertical_frame"] == "transducer"


@pytest.mark.parametrize("method", ["average_ensembles", "burst_average_ensembles"])
def test_unknown_frame_raises(rootdir, method):
    """A typo must not silently fall back to the depth grid.

    Falling back would publish a product in a frame the caller did not ask
    for, with an attribute saying it was the frame they did not want.
    """
    proc = _burst(rootdir)
    with pytest.raises(ValueError, match="vertical_frame"):
        getattr(proc, method)(vertical_frame="instrument")


def test_depth_frame_is_the_default(rootdir):
    """Omitting the keyword grids onto depth, as it always has."""
    proc = _proc(rootdir, UPLOOKER)
    proc.average_ensembles()
    assert proc.ds.attrs["vertical_frame"] == "depth"
    assert "depth" in proc.ds.dims
    assert "z" not in proc.ds.dims


# --- average_ensembles on the bin axis ----------------------------------


def test_average_ensembles_publishes_the_bin_axis(rootdir):
    """`z` is the instrument's range vector, not the universal depth grid."""
    proc = _proc(rootdir, UPLOOKER)
    proc.average_ensembles(vertical_frame="transducer")
    ds = proc.ds

    assert "z" in ds.dims
    assert "depth" not in ds.dims
    # A subset, because `_ave2nc` drops levels carrying no velocity.
    assert np.all(np.isin(ds.z.values, proc.tsdat.dep))
    assert not np.any(np.isin(ds.z.values, proc.dgrid))


def test_average_ensembles_transducer_mean_is_reproducible(rootdir):
    """The published mean is the mean of the ensemble's own bins.

    Reproduced from the pieces rather than compared against the depth-frame
    product, so this asserts that the frame switch changed the axis and
    nothing else.
    """
    proc = _proc(rootdir, UPLOOKER)
    proc.average_ensembles(vertical_frame="transducer")
    ds = proc.ds

    ens = proc.read_ensemble(0)
    proc._qc(ens)
    proc._to_enu(ens)
    expected = np.nanmean(np.ma.filled(ens.enu, np.nan), axis=0)

    idx = np.searchsorted(proc.tsdat.dep, ds.z.values)
    for j, name in enumerate(["u", "v", "w", "e"]):
        np.testing.assert_array_equal(
            ds[name].values[:, 0], expected[idx, j].astype(np.float32)
        )


def test_rejected_cells_stay_out_of_the_transducer_mean(rootdir):
    """The data under the mask of `ens.enu` is finite, not NaN.

    `numpy.ma` only flips mask bits, so a QC-rejected cell still carries the
    raw -32768 fill rotated into earth coordinates. `np.nanmean` reads no
    mask. Averaging `ens.enu` without filling to NaN first would put every
    rejected cell back into the published velocity, with nothing NaN in the
    output to show for it (issue #129).
    """
    proc = _proc(rootdir, UPLOOKER)
    ens = proc.read_ensemble(0)
    proc._qc(ens)
    proc._to_enu(ens)
    # Guard the premise: there is something to reject, and it is finite
    # under the mask.
    rejected = ~ens.valid
    assert rejected.any()
    assert np.isfinite(ens.enu.data[rejected][:, 0]).any()

    proc.average_ensembles(vertical_frame="transducer")
    ds = proc.ds

    masked_mean = np.ma.filled(ens.enu[:, :, 0].mean(axis=0), np.nan)
    idx = np.searchsorted(proc.tsdat.dep, ds.z.values)
    np.testing.assert_allclose(
        ds.u.values[:, 0], masked_mean[idx].astype(np.float32), rtol=1e-6
    )


def test_percent_good_is_a_ping_count_on_the_bin_axis(rootdir):
    """`pg` counts pings per bin, with nothing interpolated in between."""
    proc = _proc(rootdir, UPLOOKER)
    proc.average_ensembles(vertical_frame="transducer")
    ds = proc.ds

    ngood = ds.ngood.values
    npings = np.broadcast_to(ds.npings.values[np.newaxis, :], ngood.shape)
    # An interval with no pings is skipped and divides by zero; `pg` is NaN
    # there and the identity has nothing to say about it.
    good = np.isfinite(ds.pg.values) & (npings > 0)
    assert good.any()
    assert np.all(ngood[good] <= npings[good])
    expected = 100 * ngood[good] // npings[good]
    np.testing.assert_array_equal(ds.pg.values[good], expected)


# --- burst_average_ensembles on the bin axis ----------------------------


def test_burst_transducer_product_is_the_burst_mean(rootdir):
    """The published burst mean is the array the depth frame grids from.

    `burst_average_ensembles` already forms the mean on the instrument bins
    as `uvwe_inst` and then interpolates it onto `self.dgrid`. In the
    transducer frame it stops there, so the two must agree bitwise.
    """
    proc = _burst(rootdir)
    proc.burst_average_ensembles(vertical_frame="transducer")
    ds = proc.ds

    ens = proc.read_ensemble(0)
    proc._qc(ens)
    proc._to_enu(ens)
    expected = ens.enu.mean(axis=0)
    ngood = np.sum(ens.valid, axis=0)
    pgi = 100 * ngood // ens.enu.shape[0]
    expected[pgi < proc.editparams.pg_limit, :] = np.ma.masked

    idx = np.searchsorted(proc.tsdat.dep, ds.z.values)
    np.testing.assert_array_equal(
        ds.u.values[:, 0],
        np.ma.filled(expected[idx, 0], np.nan).astype(np.float32),
    )


def test_burst_pg_limit_and_masking_agree_on_the_bin_axis(rootdir):
    """On the bin axis, screened and masked are the same condition, both ways.

    `pg_limit` masks the bins whose percent good falls below it. In the
    transducer frame `pg` is that same per-bin percent good, so `pg` below
    the limit and a masked velocity are one condition. In the depth frame
    they come apart: `interp1` widens the mask of a screened bin onto both
    grid cells touching it, while the interpolated `pg` at those cells is a
    blend with unscreened neighbours and can still read at or above the
    limit. Asserting only that a screened cell is masked would not
    distinguish the frames, because that direction holds in both.
    """
    proc = _burst(rootdir)
    proc.burst_average_ensembles(vertical_frame="transducer")
    limit = proc.editparams.pg_limit
    pg = proc.ds.pg.values
    u = proc.ds.u.values

    # Transducer frame: the biconditional holds in both directions.
    sampled = np.isfinite(pg)  # exclude intervals with fewer than two pings
    screened = sampled & (pg < limit)
    kept = sampled & (pg >= limit)
    assert screened.any() and kept.any(), "pg_limit should bind somewhere on this file"
    assert np.all(np.isnan(u[screened]))
    assert np.all(np.isfinite(u[kept]))  # the direction that is frame-specific

    # Depth frame: kept cells can still be masked, because a screened
    # neighbour's mask reaches them without dragging pg below the limit.
    proc_depth = _burst(rootdir)
    proc_depth.burst_average_ensembles()
    pg_depth = proc_depth.ds.pg.values
    u_depth = proc_depth.ds.u.values
    sampled_depth = np.isfinite(pg_depth)
    kept_depth = sampled_depth & (pg_depth >= limit)
    assert np.any(np.isnan(u_depth[kept_depth]))


def test_burst_percent_good_identity_is_frame_dependent(rootdir):
    """`pg == 100 * ngood // npings` holds on the bin axis and not on the grid.

    In the depth frame `pg` and `ngood` are interpolated onto the grid and
    floored independently, so the identity breaks. That contrast is the
    clearest statement of what the frame changes about the counts.
    """

    def _identity(ds):
        ngood = ds.ngood.values
        npings = np.broadcast_to(ds.npings.values[np.newaxis, :], ngood.shape)
        good = np.isfinite(ds.pg.values) & np.isfinite(ngood) & (npings > 0)
        assert good.any()
        return ds.pg.values[good], 100 * ngood[good] // npings[good]

    proc = _burst(rootdir)
    proc.burst_average_ensembles(vertical_frame="transducer")
    published, expected = _identity(proc.ds)
    np.testing.assert_array_equal(published, expected)

    proc = _burst(rootdir)
    proc.burst_average_ensembles()
    published, expected = _identity(proc.ds)
    assert not np.array_equal(published, expected)


def test_burst_interpolate_bin_works_on_the_bin_axis(rootdir):
    """`interpolate_bin` is bin-referenced and needs no change.

    It fills the bin from its neighbours before any gridding, so in the
    transducer frame the filled bin reaches the output directly. `pg` is
    deliberately left at 0 there, so the interpolation stays visible.
    """
    proc = _burst(rootdir, editparams={"maskbins": [5]})
    proc.burst_average_ensembles(vertical_frame="transducer", interpolate_bin=5)
    ds = proc.ds

    assert proc.tsdat.dep[5] in ds.z.values
    k = int(np.flatnonzero(ds.z.values == proc.tsdat.dep[5])[0])
    assert np.isfinite(ds.u.values[k]).any()
    assert np.all(ds.pg.values[k][np.isfinite(ds.pg.values[k])] == 0)
    assert ds.attrs["interpolate_bin"] == 5


# --- the derived bin depth ----------------------------------------------


def _uplooker_averaged(rootdir):
    proc = _proc(rootdir, UPLOOKER)
    proc.average_ensembles(vertical_frame="transducer")
    return proc


def _downlooker_burst(rootdir):
    proc = _burst(rootdir)
    proc.burst_average_ensembles(vertical_frame="transducer")
    return proc


@pytest.mark.parametrize(
    ("build", "orientation"),
    [(_uplooker_averaged, "up"), (_downlooker_burst, "down")],
)
def test_bin_depth_is_exactly_the_reconstruction(rootdir, build, orientation):
    """`depth` is `xducer_depth` plus or minus `z`, to the last bit.

    Written so a consumer outside this project does not reimplement the sign
    rule, which is what makes it worth an extra 2-D field. If it were not
    exact, reconstructing would be the better advice and the field would not
    be worth writing. Both orientations, because the sign is the part a
    reader gets wrong.
    """
    proc = build(rootdir)
    ds = proc.ds
    assert proc.orientation == orientation
    sign = -1.0 if orientation == "up" else 1.0

    assert ds.depth.dims == ("z", "time")
    assert "depth" not in ds.dims
    expected = ds.xducer_depth.values[np.newaxis, :] + sign * ds.z.values[:, None]
    np.testing.assert_array_equal(ds.depth.values, expected)


def test_bin_depth_says_it_is_derived(rootdir):
    """The comment has to name the reconstruction and the gsw/seawater gap."""
    proc = _proc(rootdir, UPLOOKER)
    proc.average_ensembles(vertical_frame="transducer")
    comment = proc.ds.depth.attrs["comment"]

    assert "xducer_depth" in comment
    assert "gsw" in comment and "seawater" in comment
    assert proc.ds.depth.attrs["standard_name"] == "depth"
    assert proc.ds.depth.attrs["positive"] == "down"


def test_depth_frame_carries_no_derived_depth(rootdir):
    """A depth-gridded product has a 1-D depth dimension and no comment."""
    proc = _proc(rootdir, UPLOOKER)
    proc.average_ensembles()
    assert proc.ds.depth.dims == ("depth",)
    assert "comment" not in proc.ds.depth.attrs


def test_process_pings_carries_no_derived_depth(rootdir):
    """The single-ping product is the largest this package writes.

    A float64 `depth(z, time)` beside it would more than double it, so the
    reconstruction stays the reader's job there, as the `z` attributes say.
    """
    proc = _proc(rootdir, UPLOOKER)
    proc.process_pings()
    assert "depth" not in proc.ds.variables


def test_depth_offset_moves_the_depths_and_not_the_axis(rootdir):
    """An offset on the instrument's depth says nothing about `z`.

    `z` is a distance from the transducer and is unmoved. `xducer_depth` and
    the derived `depth` both take the offset, and no velocity changes,
    because the offset enters the output only (issue #92).
    """
    offset = 5.0
    plain = _proc(rootdir, UPLOOKER)
    plain.average_ensembles(vertical_frame="transducer")
    shifted = _proc(rootdir, UPLOOKER, depth_offset=offset)
    shifted.average_ensembles(vertical_frame="transducer")

    np.testing.assert_array_equal(plain.ds.z.values, shifted.ds.z.values)
    np.testing.assert_array_equal(plain.ds.u.values, shifted.ds.u.values)
    np.testing.assert_allclose(
        shifted.ds.xducer_depth.values, plain.ds.xducer_depth.values + offset
    )
    np.testing.assert_allclose(
        shifted.ds.depth.values, plain.ds.depth.values + offset
    )
    assert shifted.ds.attrs["depth_offset"] == offset


def test_bin_depth_survives_a_netcdf_round_trip(rootdir, tmp_path):
    """A non-dimension coordinate has to reach the file as a coordinate."""
    import xarray as xr

    proc = _proc(rootdir, UPLOOKER)
    proc.average_ensembles(vertical_frame="transducer")
    f = tmp_path / "z.nc"
    proc.ds.to_netcdf(f)
    back = xr.open_dataset(f)

    assert "depth" in back.coords
    assert back.depth.dims == ("z", "time")
    np.testing.assert_array_equal(back.depth.values, proc.ds.depth.values)


# --- provenance comments ------------------------------------------------


@pytest.mark.parametrize("var", ["pg", "ngood", "nbad_cor"])
def test_counts_say_they_are_ping_counts_on_the_bin_axis(rootdir, var):
    """The count comments are frame-dependent, not only path-dependent."""
    plain = _proc(rootdir, UPLOOKER)
    plain.average_ensembles()
    shifted = _proc(rootdir, UPLOOKER)
    shifted.average_ensembles(vertical_frame="transducer")

    assert plain.ds[var].attrs["comment"] != shifted.ds[var].attrs["comment"]
    assert "average_ensembles" in shifted.ds[var].attrs["comment"]


def test_ngood_does_not_claim_an_off_profile_nan_on_the_bin_axis(rootdir):
    """`io.cf_conventions` gives `ngood` a static off-profile sentence.

    There is no off-profile cell on the bin axis, so that sentence has to be
    replaced rather than appended to.
    """
    proc = _proc(rootdir, UPLOOKER)
    proc.average_ensembles(vertical_frame="transducer")
    comment = proc.ds.ngood.attrs["comment"]
    assert "outside the instrument's profile" not in comment


def test_z_attributes_describe_both_averaging_frames(rootdir):
    """The `z` comment used to say the averaging methods always regrid."""
    proc = _proc(rootdir, UPLOOKER)
    proc.average_ensembles(vertical_frame="transducer")
    comment = proc.ds.z.attrs["comment"]
    assert "average_ensembles" in comment
    assert "burst_average_ensembles" in comment
