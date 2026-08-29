"""Tests for the constant depth offset (issue #92).

An ADCP whose pressure sensor carries a constant bias produces depths that
are offset by a constant amount while the velocities themselves are fine.
`depth_offset` translates the published depth axis by that constant.

The whole point is that it is exactly the post-hoc shift a user would apply
to the finished file, so the velocity field must come out bit-identical to an
unoffset run. That is achieved by gridding entirely in raw pressure-derived
depth and adding the offset only on the way out. Everything below exists to
pin that down, in particular the consequence that a user-supplied
`dgridparams` is read in the RAW frame.
"""

import gsw
import numpy as np
import pytest
import xarray as xr

from velosearaptor import io
from velosearaptor.madcp import ProcessADCP

META_DATA = {
    "mooring": "Test",
    "project": "Test",
    "lon": 0,
    "lat": 0,
}

UPLOOKER = "binmap_16670013.000"
DOWNLOOKER = "24606000.000"

# The offset measured on MOTIVE mooring B in issue #92. Deliberately not a
# round number and not a multiple of the 4 m bin size, so an axis that came
# out of a re-gridded run rather than a translated one would not match.
OFFSET = 5.087


def _proc(rootdir, filename, **kwargs):
    return ProcessADCP(rootdir / "data" / filename, META_DATA, magdec=0.0, **kwargs)


def _averaged(rootdir, filename, **kwargs):
    proc = _proc(rootdir, filename, **kwargs)
    proc.average_ensembles()
    return proc


def _filled(da):
    """Values with NaN made comparable, for bitwise comparison."""
    return np.ma.filled(np.ma.masked_invalid(da.values), -9999.0)


def test_depth_offset_defaults_to_zero(rootdir):
    """Without the argument the output is what it was before issue #92."""
    ds = _averaged(rootdir, UPLOOKER).ds

    assert ds.attrs["depth_offset"] == 0.0
    # No derived pressure variable is invented for a run that corrects nothing.
    assert "pressure_corrected" not in ds.variables


def test_depth_offset_zero_matches_no_argument(rootdir):
    """Passing 0.0 explicitly is the same run, values and axis alike."""
    plain = _averaged(rootdir, UPLOOKER).ds
    zero = _averaged(rootdir, UPLOOKER, depth_offset=0.0).ds

    assert np.array_equal(zero.depth.values, plain.depth.values)
    assert "pressure_corrected" not in zero.variables
    for var in ["u", "v", "w", "pg"]:
        assert np.array_equal(_filled(zero[var]), _filled(plain[var]))


@pytest.mark.parametrize("filename", [UPLOOKER, DOWNLOOKER])
def test_velocity_is_bitwise_identical_under_offset(rootdir, filename):
    """The load-bearing test: an offset run changes the axis and nothing else.

    This is what the post-hoc shift gives and what supplying a corrected
    external pressure series cannot give, because `depth(p)` is nonlinear and
    re-deriving per-ping depths re-samples the interpolation onto the grid.
    """
    plain = _averaged(rootdir, filename).ds
    offset = _averaged(rootdir, filename, depth_offset=OFFSET).ds

    for var in ["u", "v", "w", "pg"]:
        assert np.array_equal(_filled(offset[var]), _filled(plain[var])), (
            f"{var} is not bit-identical under a pure depth translation"
        )


@pytest.mark.parametrize("filename", [UPLOOKER, DOWNLOOKER])
def test_depth_axis_is_translated_exactly(rootdir, filename):
    """`depth` and `xducer_depth` move by exactly the offset, nothing else."""
    plain = _averaged(rootdir, filename).ds
    offset = _averaged(rootdir, filename, depth_offset=OFFSET).ds

    assert offset.depth.size == plain.depth.size
    assert np.array_equal(offset.depth.values, plain.depth.values + OFFSET)
    assert np.array_equal(
        _filled(offset.xducer_depth), _filled(plain.xducer_depth + OFFSET)
    )

    # The raw measurement is untouched: it is what the sensor read.
    assert np.array_equal(_filled(offset.pressure), _filled(plain.pressure))


def test_internal_grid_is_untouched_by_the_offset(rootdir):
    """The offset never reaches the grid the interpolation targets.

    This is the mechanism behind bit-identity, asserted directly rather than
    only through its consequence.
    """
    plain = _proc(rootdir, UPLOOKER)
    offset = _proc(rootdir, UPLOOKER, depth_offset=OFFSET)

    assert np.array_equal(offset.dgrid, plain.dgrid)
    assert offset.default_dgridparams == plain.default_dgridparams
    assert offset.p_median == plain.p_median


def test_user_dgridparams_are_read_in_the_raw_frame(rootdir):
    """A supplied `dtop`/`dbot` means raw pressure-derived depth.

    The alternative reading, that the user means corrected depth, would move
    the internal grid, point `interp1` at different depths and break the
    bit-identity above whenever `dgridparams` is supplied. The frame is
    documented in the class docstring and recorded in the file attributes.
    """
    dgridparams = {"dtop": 100.0, "dbot": 160.0, "d_interval": 4.0}
    plain = _averaged(rootdir, UPLOOKER, dgridparams=dgridparams)
    offset = _averaged(rootdir, UPLOOKER, dgridparams=dgridparams, depth_offset=OFFSET)

    # Same internal grid, so the same velocities, exactly as without
    # dgridparams.
    assert np.array_equal(offset.dgrid, plain.dgrid)
    for var in ["u", "v", "w", "pg"]:
        assert np.array_equal(_filled(offset.ds[var]), _filled(plain.ds[var]))

    # The requested values live on the raw axis, so the published axis is the
    # requested one translated. Asking for 100 and reading 105.087 is the
    # documented behavior, not a bug.
    assert np.all(plain.ds.depth.values >= dgridparams["dtop"])
    assert np.all(plain.ds.depth.values <= dgridparams["dbot"])
    assert np.array_equal(offset.ds.depth.values, plain.ds.depth.values + OFFSET)
    assert offset.ds.depth.values.min() > dgridparams["dtop"]


def test_process_pings_shifts_xducer_depth_but_not_z(rootdir):
    """`z` is a distance from the transducer, so an offset cannot apply to it.

    Only the transducer's own depth moves. Bin depth stays recoverable as
    `xducer_depth` -/+ `z`, now in the corrected frame.
    """
    plain = _proc(rootdir, UPLOOKER)
    plain.process_pings()
    offset = _proc(rootdir, UPLOOKER, depth_offset=OFFSET)
    offset.process_pings()

    assert "z" in offset.ds.coords
    assert "depth" not in offset.ds.coords
    assert np.array_equal(offset.ds.z.values, plain.ds.z.values)
    assert np.array_equal(
        _filled(offset.ds.xducer_depth), _filled(plain.ds.xducer_depth + OFFSET)
    )
    for var in ["u", "v", "w", "pg"]:
        assert np.array_equal(_filled(offset.ds[var]), _filled(plain.ds[var]))


def test_burst_average_shifts_the_depth_axis(rootdir):
    """The third path gets the same treatment as `average_ensembles`."""
    kwargs = {"tgridparams": {"burst_average": True}}
    plain = _proc(rootdir, DOWNLOOKER, **kwargs)
    plain.burst_average_ensembles()
    offset = _proc(rootdir, DOWNLOOKER, depth_offset=OFFSET, **kwargs)
    offset.burst_average_ensembles()

    assert np.array_equal(offset.ds.depth.values, plain.ds.depth.values + OFFSET)
    for var in ["u", "v", "w"]:
        assert np.array_equal(_filled(offset.ds[var]), _filled(plain.ds[var]))


def test_corrected_pressure_is_consistent_with_the_stored_depth(rootdir):
    """`pressure_corrected` is what a sensor at the corrected depth would read.

    `pressure` stays the measurement. Recomputing depth from
    `pressure_corrected` reproduces the published `xducer_depth`, which is the
    inconsistency issue #92 complains about, closed.
    """
    proc = _averaged(rootdir, UPLOOKER, depth_offset=OFFSET)
    ds = proc.ds

    assert "pressure_corrected" in ds.variables
    assert ds.pressure_corrected.dims == ("time",)

    recovered = -gsw.z_from_p(ds.pressure_corrected.values, proc.lat)
    assert np.allclose(recovered, ds.xducer_depth.values, atol=1e-6)

    # It is not the raw pressure plus a constant: depth(p) is nonlinear, so
    # the dbar difference is not constant with depth. That is precisely why
    # the raw measurement is kept alongside it.
    diff = ds.pressure_corrected.values - ds.pressure.values
    assert np.all(diff > 0)
    assert diff.std() > 0


def test_corrected_pressure_has_cf_attributes(rootdir):
    """The derived variable carries its own entry, and says it is derived."""
    ds = _averaged(rootdir, UPLOOKER, depth_offset=OFFSET).ds
    CF = io.cf_conventions()

    assert "pressure_corrected" in CF
    attrs = ds.pressure_corrected.attrs
    assert attrs["units"] == "dbar"
    comment = attrs["comment"].lower()
    assert "derived" in comment
    assert "nonlinear" in comment
    # It must not be mistaken for the measurement.
    assert "`pressure`" in attrs["comment"]


def test_offset_is_recorded_on_depth_and_xducer_depth(rootdir):
    """A reader must be able to see which frame the axis is in."""
    ds = _averaged(rootdir, UPLOOKER, depth_offset=OFFSET).ds

    for name in ["depth", "xducer_depth"]:
        comment = ds[name].attrs.get("comment", "")
        assert "depth_offset" in comment, f"{name} does not name the offset"
        assert str(OFFSET) in comment

    # An unoffset file must not carry the claim.
    plain = _averaged(rootdir, UPLOOKER).ds
    assert "depth_offset" not in plain.depth.attrs.get("comment", "")


def test_offset_round_trips_through_netcdf(rootdir, tmp_path):
    """The offset and its sign convention survive a write and a read."""
    ds = _averaged(rootdir, UPLOOKER, depth_offset=OFFSET).ds
    path = tmp_path / "offset.nc"
    ds.to_netcdf(path)

    with xr.open_dataset(path) as reopened:
        assert reopened.attrs["depth_offset"] == OFFSET
        assert "depth_offset" in reopened.depth.attrs["comment"]
        assert np.array_equal(reopened.depth.values, ds.depth.values)


def test_sign_convention_is_added_not_subtracted(rootdir):
    """Positive `depth_offset` means the instrument was deeper than measured.

    Stated as its own test because a sign error here is silent: the axis
    still looks plausible, it is just wrong by twice the offset.
    """
    plain = _averaged(rootdir, UPLOOKER).ds
    offset = _averaged(rootdir, UPLOOKER, depth_offset=OFFSET).ds

    assert offset.depth.values.min() > plain.depth.values.min()
    assert float(offset.xducer_depth.mean()) > float(plain.xducer_depth.mean())

    negative = _averaged(rootdir, UPLOOKER, depth_offset=-OFFSET).ds
    assert np.array_equal(negative.depth.values, plain.depth.values - OFFSET)
