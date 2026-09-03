"""Tests for the transducer-relative vertical frame (issue #129).

`average_ensembles` and `burst_average_ensembles` interpolate every ping onto
`self.dgrid` before averaging. `vertical_frame="transducer"` keeps the
instrument's own bin axis instead and publishes it as `z`, the axis
`process_pings` already publishes, with the bin depth of each averaging
interval alongside it as a derived coordinate.
"""

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
