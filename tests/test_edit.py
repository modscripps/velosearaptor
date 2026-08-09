"""Tests for editing/QC behavior in velosearaptor.madcp.ProcessADCP."""

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
def adcpfile(rootdir):
    return rootdir / "data/binmap_16670013.000"


def test_binmap_honors_ibad(adcpfile):
    """A beam excluded via ibad must stay excluded when binmapping
    recomputes the beam-to-xyz transform (issue #79)."""
    proc_4beam = ProcessADCP(adcpfile, META_DATA, magdec=0.0)
    proc_4beam.process_pings(binmap=True)
    u_4beam = proc_4beam.ds.u.values

    proc_3beam = ProcessADCP(adcpfile, META_DATA, magdec=0.0, ibad=0)
    proc_3beam.process_pings(binmap=True)
    u_3beam = proc_3beam.ds.u.values

    # A 3-beam solution must differ from the 4-beam solution. With ibad
    # dropped in the recomputed transform, the two results are identical.
    both_finite = np.isfinite(u_4beam) & np.isfinite(u_3beam)
    assert both_finite.sum() > 0
    assert not np.allclose(u_4beam[both_finite], u_3beam[both_finite])
