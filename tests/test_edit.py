"""Tests for editing/QC behavior in velosearaptor.madcp.ProcessADCP."""

import numpy as np
import pytest
from pycurrents.system import Bunch

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


def _ensemble(dead_beam_cor, min_cor_beam=0, nping=4, nbin=5):
    """A synthetic ensemble whose beam `min_cor_beam` has a low correlation
    and whose other three beams are well above any threshold. Error velocity
    is identically zero so only the correlation test can mask anything."""
    ens = Bunch()
    ens.cor = np.full((nping, nbin, 4), 128.0)
    ens.cor[..., min_cor_beam] = dead_beam_cor
    ens.xyze = np.ma.masked_array(np.zeros((nping, nbin, 4)), mask=False)
    return ens


def test_correlation_test_skips_the_ibad_beam(adcpfile):
    """A beam excluded via `ibad` must not reject cells through the
    correlation test either (issue #90)."""
    proc = ProcessADCP(adcpfile, META_DATA, magdec=0.0, ibad=0)
    proc.parse_editparams({"min_correlation": 64})

    ens = _ensemble(dead_beam_cor=10.0, min_cor_beam=0)
    proc._edit(ens)

    assert not np.ma.getmaskarray(ens.xyze).any()


def test_correlation_test_uses_every_beam_without_ibad(adcpfile):
    """Without `ibad` the any-beam behavior is unchanged: one weak beam
    still rejects the cell."""
    proc = ProcessADCP(adcpfile, META_DATA, magdec=0.0)
    proc.parse_editparams({"min_correlation": 64})

    ens = _ensemble(dead_beam_cor=10.0, min_cor_beam=0)
    proc._edit(ens)

    assert np.ma.getmaskarray(ens.xyze).all()


def test_correlation_test_still_rejects_a_good_beam_dropout(adcpfile):
    """`ibad` narrows the test to the remaining beams; a weak beam among
    those still rejects the cell."""
    proc = ProcessADCP(adcpfile, META_DATA, magdec=0.0, ibad=0)
    proc.parse_editparams({"min_correlation": 64})

    ens = _ensemble(dead_beam_cor=10.0, min_cor_beam=2)
    proc._edit(ens)

    assert np.ma.getmaskarray(ens.xyze).all()


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
