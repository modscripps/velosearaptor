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
