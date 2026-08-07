"""Tests for ProcessADCP._ensure_monotonic_dday."""

import numpy as np
import pytest
from pycurrents.system import Bunch

from velosearaptor.madcp import ProcessADCP


def _make_stub(dday):
    """Create a minimal ProcessADCP stub with dday and tsdat set."""
    obj = object.__new__(ProcessADCP)
    obj.dday = np.array(dday, dtype=float)
    obj.tsdat = Bunch(
        dday=np.array(dday, dtype=float),
        pressure=np.ones(len(dday)),
        temperature=np.ones(len(dday)) * 15.0,
        ens_num=np.arange(len(dday)),
    )
    return obj


def test_already_monotonic():
    obj = _make_stub([0.0, 0.1, 0.2, 0.3, 0.4])
    obj._ensure_monotonic_dday()
    np.testing.assert_array_equal(obj.dday, [0.0, 0.1, 0.2, 0.3, 0.4])


def test_single_isolated_bad_ping():
    # One ping jumps backward — should be interpolated
    dday = [0.0, 0.1, 0.05, 0.3, 0.4]
    obj = _make_stub(dday)
    obj._ensure_monotonic_dday()
    # Index 2 should be interpolated between 0.1 and 0.3
    assert obj.dday[2] == pytest.approx(0.2)
    assert len(obj.dday) == 5  # length preserved


def test_segment_overlap_truncation():
    # Time resets and never recovers — should truncate
    dday = [0.0, 0.1, 0.2, 0.3, 0.05, 0.06, 0.07, 0.08]
    obj = _make_stub(dday)
    obj._ensure_monotonic_dday()
    np.testing.assert_array_equal(obj.dday, [0.0, 0.1, 0.2, 0.3])
    assert len(obj.tsdat.dday) == 4
    assert len(obj.tsdat.pressure) == 4
    assert len(obj.tsdat.temperature) == 4
    assert len(obj.tsdat.ens_num) == 4


def test_ambiguous_raises():
    # Backward jump that eventually recovers past the pre-jump max
    dday = [0.0, 0.1, 0.2, 0.3, 0.05, 0.06, 0.07, 0.08, 0.09, 0.35, 0.4]
    obj = _make_stub(dday)
    with pytest.raises(ValueError, match="Cannot auto-fix"):
        obj._ensure_monotonic_dday()


def test_bad_ping_at_end():
    # Last ping is non-monotonic — should be interpolated using median dt
    dday = [0.0, 0.1, 0.2, 0.3, 0.1]
    obj = _make_stub(dday)
    obj._ensure_monotonic_dday()
    # Should use median of positive diffs (0.1) added to previous value (0.3)
    assert obj.dday[4] == pytest.approx(0.4)
    assert len(obj.dday) == 5
