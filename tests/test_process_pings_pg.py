"""Test percent good in single-ping processing (issue #77).

For single-ping data, percent good is binary: 100 where the ping survived
editing, 0 where it was edited out. It must not be zero everywhere.
"""

import numpy as np

from velosearaptor.madcp import ProcessADCP

META_DATA = {
    "mooring": "Test",
    "project": "Test",
    "lon": 0,
    "lat": 0,
}


def test_process_pings_pg(rootdir):
    adcpfile = rootdir / "data/binmap_16670013.000"
    proc = ProcessADCP(adcpfile, META_DATA, magdec=0.0)
    proc.process_pings()
    ds = proc.ds

    good = np.isfinite(ds.u.values)
    pg = ds.pg.values

    assert good.sum() > 0
    assert np.all(pg[good] == 100)

    # Where amplitude is present but the velocity was edited out, pg is 0.
    edited = np.isfinite(ds.amp.values) & ~good
    assert edited.sum() > 0
    assert np.all(pg[edited] == 0)
