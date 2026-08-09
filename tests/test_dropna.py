"""Tests for empty-level dropping in _ave2nc (issue #84).

A depth level where amplitude is finite but every velocity is NaN must not
survive into the output dataset: amp is not edited, so it can be finite in
levels that carry no velocity data at all.
"""

import numpy as np
from pycurrents.system import Bunch

from velosearaptor.madcp import ProcessADCP


def _bare_instance_with_ave():
    """ProcessADCP shell with a synthetic `ave`: three depth levels, of
    which one is fully valid, one has amp only, one is fully empty."""
    p = ProcessADCP.__new__(ProcessADCP)
    p.lat = 0.0

    ntime, ndep = 2, 3
    u = np.full((ntime, ndep), 0.5, dtype=np.float32)
    amp = np.full((ntime, ndep), 100.0, dtype=np.float32)
    u[:, 1] = np.nan  # level 1: velocity gone, amp still finite
    u[:, 2] = np.nan  # level 2: nothing at all
    amp[:, 2] = np.nan

    p.ave = Bunch(
        u=u.copy(),
        v=u.copy(),
        w=u.copy(),
        e=u.copy(),
        pg=np.full((ntime, ndep), 100, dtype=np.int8),
        amp=amp,
        temperature=np.full((ntime,), 10.0, dtype=np.float32),
        pressure=np.full((ntime,), 1000.0, dtype=np.float32),
        dday=np.array([0.0, 1 / 24]),
        yearbase=2020,
        dep=np.array([100.0, 110.0, 120.0]),
    )
    return p


def test_amp_only_level_is_dropped():
    p = _bare_instance_with_ave()
    p._ave2nc()
    assert list(p.ds.depth.values) == [100.0]
