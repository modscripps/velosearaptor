"""Tests for the vertical axis of the single-ping product (issue #96).

`process_pings` does not regrid, so its vertical axis is the distance from
the transducer to the center of each bin. That is not water depth, and it is
signed the other way for an uplooker. It is published as `z`; `depth` is
reserved for the averaging methods, which do grid onto water depth.
"""

import re

import numpy as np
import pytest
from pycurrents.data import seawater

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


def _processed(rootdir, filename):
    proc = ProcessADCP(rootdir / "data" / filename, META_DATA, magdec=0.0)
    proc.process_pings()
    return proc


def test_process_pings_publishes_z_not_depth(rootdir):
    """The single-ping vertical axis is `z`, and there is no `depth` axis."""
    ds = _processed(rootdir, UPLOOKER).ds

    assert "z" in ds.coords
    assert "z" in ds.dims
    assert "depth" not in ds.coords
    assert "depth" not in ds.dims
    assert "depth" not in ds.variables


def test_z_carries_cf_z_attributes(rootdir):
    """`z` gets the `z` entry of cf_conventions, not the `depth` entry."""
    ds = _processed(rootdir, UPLOOKER).ds
    CF = io.cf_conventions()

    assert ds.z.attrs == CF["z"]
    assert "standard_name" not in ds.z.attrs
    assert "positive" not in ds.z.attrs

    # `standard_name: depth` belongs to the transducer depth alone on this
    # path. Nothing else in the file is a water depth.
    named_depth = {
        v for v in ds.variables if ds[v].attrs.get("standard_name") == "depth"
    }
    assert named_depth == {"xducer_depth"}


def test_z_attrs_do_not_claim_an_absent_variable(rootdir):
    """Names the `z` attributes quote must exist in the file or be methods.

    The `z` comment describes how to recover bin depth. Every name it puts in
    backticks has to be something the reader can actually find: a variable of
    this dataset, or a `ProcessADCP` method to call. A reference to a `depth`
    variable that the single-ping product does not carry is exactly the claim
    this test exists to catch.
    """
    ds = _processed(rootdir, UPLOOKER).ds

    quoted = set()
    for value in ds.z.attrs.values():
        if isinstance(value, str):
            quoted.update(re.findall(r"`([A-Za-z_][A-Za-z0-9_]*)`", value))

    assert quoted, "the z attributes should name how to recover bin depth"
    for name in quoted:
        assert name in ds.variables or hasattr(ProcessADCP, name), (
            f"`{name}` is neither a variable of this dataset nor a ProcessADCP method"
        )
    assert "xducer_depth" in quoted


@pytest.mark.parametrize(
    ("filename", "orientation"),
    [(UPLOOKER, "up"), (DOWNLOOKER, "down")],
)
def test_xducer_depth_plus_z_reproduces_ping_depth(rootdir, filename, orientation):
    """`xducer_depth` +/- `z` recovers the per-ping bin depth.

    This is what makes the file self-describing after the rename: the depth
    of every bin is recoverable, including under knockdown, and the sign
    depends on the instrument orientation.
    """
    proc = _processed(rootdir, filename)
    ds = proc.ds
    assert proc.orientation == orientation
    sign = -1 if orientation == "up" else 1

    # Read one raw ping and compute its bin depths the way process_pings does.
    idx_start = np.searchsorted(proc.dday, proc.dday_start)
    k = 7
    ens = proc.m.read(start=idx_start + k, stop=idx_start + k + 1)
    ens.pressure = proc._scale_pycurrents_pressure(ens)
    ens_depth = seawater.depth2(ens.pressure, proc.lat)[:, np.newaxis] + sign * ens.dep
    # Guard that ping k of the output really is this raw ping.
    assert np.isclose(float(ds.pressure[k]), ens.pressure[0])

    # Empty levels are dropped from the output, so the published axis is a
    # subset of the instrument's bin distances.
    idx = np.searchsorted(ens.dep, ds.z.values)
    assert np.allclose(ens.dep[idx], ds.z.values)

    recovered = float(ds.xducer_depth[k]) + sign * ds.z.values
    # xducer_depth comes from gsw.z_from_p while the per-ping depth uses
    # seawater.depth2; the two conversions differ by a few cm here.
    assert np.max(np.abs(recovered - ens_depth[0, idx])) < 0.5

    # The other sign is wrong by twice the bin distance, which is metres.
    flipped = float(ds.xducer_depth[k]) - sign * ds.z.values
    assert np.min(np.abs(flipped - ens_depth[0, idx])) > 5


def test_averaging_path_still_publishes_depth(rootdir):
    """Regression guard: `average_ensembles` grids onto water depth."""
    proc = ProcessADCP(rootdir / "data" / UPLOOKER, META_DATA, magdec=0.0)
    proc.average_ensembles()
    ds = proc.ds
    CF = io.cf_conventions()

    assert "depth" in ds.coords
    assert "z" not in ds.coords
    assert "z" not in ds.variables
    assert ds.depth.attrs == CF["depth"]
    # The averaged axis is the universal depth grid, not the instrument's
    # bin distances.
    assert np.all(np.isin(ds.depth.values, proc.dgrid))
