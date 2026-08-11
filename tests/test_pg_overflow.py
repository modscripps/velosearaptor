"""Test percent good does not overflow at high ping counts (issue #94).

`average_ensembles` computed percent good as `100 * ngood // nprofs` with
`ngood` stored as `int16`. The product `100 * ngood` stays `int16`, so it wraps
once `ngood` exceeds 327 (100 * 328 = 32800 > 32767). Floor division then
returns a *negative* percent good for the bins with the *most* good pings, and
any `pg` filter discards precisely the best data. `.astype(np.int8)` on the
following line preserves the wrapped value, so nothing raised.

The fixture pings every 2.5 s, so a 0.25 h averaging interval collects exactly
360 pings per ensemble — past the wrap threshold. At 360 pings a fully good bin
reported -83 instead of 100.
"""

import numpy as np

from velosearaptor.madcp import ProcessADCP

META_DATA = {
    "mooring": "Test",
    "project": "Test",
    "lon": 0,
    "lat": 0,
}


def test_pg_in_range_above_int16_wrap_threshold(rootdir):
    adcpfile = rootdir / "data/binmap_16670013.000"
    proc = ProcessADCP(adcpfile, META_DATA, magdec=0.0, tgridparams={"dt_hours": 0.25})
    proc.average_ensembles()
    ds = proc.ds

    pg = ds.pg.values
    ngood = ds.ngood.values
    npings = np.broadcast_to(ds.npings.values[np.newaxis, :], pg.shape)

    # The fixture has to actually cross the threshold, or this tests nothing.
    assert npings.max() > 327

    # pg is a percentage. Nothing outside 0-100 is meaningful.
    assert pg.min() >= 0
    assert pg.max() <= 100

    # Every ping good means 100 % good, whatever the ping count. This is the
    # assertion the overflow broke: at 360 pings it returned -83.
    full = ngood == npings
    assert full.sum() > 0
    assert np.all(pg[full] == 100)

    # pg must equal the percentage computed in arithmetic wide enough not to
    # wrap. Computing the expectation in int16 would reproduce the bug instead
    # of catching it.
    expected = 100 * ngood.astype(np.int64) // np.where(npings > 0, npings, 1)
    assert np.array_equal(pg.astype(np.int64), expected)
