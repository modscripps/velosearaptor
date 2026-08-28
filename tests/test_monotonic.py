"""Tests for ProcessADCP._ensure_monotonic_dday."""

import numpy as np
import pytest
from pycurrents.system import Bunch

from velosearaptor.madcp import ProcessADCP


def _make_stub(dday, time_drift_rate=1.0):
    """Create a minimal ProcessADCP stub with dday and tsdat set.

    `tsdat.dday` holds the raw time base and `self.dday` the clock-drift
    corrected copy of it, which is the relationship `parse_driftparams`
    establishes: ``self.dday = self._correct_dday(self.tsdat.dday)``. They are
    separate arrays, so a repair of one does not reach the other.
    """
    obj = object.__new__(ProcessADCP)
    raw = np.array(dday, dtype=float)
    obj.t0 = raw[0]
    obj.time_drift_rate = time_drift_rate
    obj.tsdat = Bunch(
        dday=raw.copy(),
        pressure=np.ones(len(dday)),
        temperature=np.ones(len(dday)) * 15.0,
        ens_num=np.arange(len(dday)),
    )
    obj.dday = obj._correct_dday(raw)
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


# --- issue #98: the repair must actually achieve monotonicity ---------------

# Cases that the classification either repairs into a strictly increasing time
# vector or gives up on with a ValueError.  Reporting success while leaving the
# record non-monotonic is the bug.
_NON_MONOTONIC_CASES = {
    "forward_spike": [0.0, 0.1, 0.9, 0.3, 0.4],
    "isolated_jump_then_segment_overlap": [
        0.0,
        0.1,
        0.05,
        0.2,
        0.3,
        0.4,
        0.1,
        0.11,
        0.12,
        0.13,
    ],
    "duplicate_timestamps": [0.0, 0.1, 0.1, 0.1, 0.2, 0.3],
}


@pytest.mark.parametrize("name", sorted(_NON_MONOTONIC_CASES))
def test_repair_is_monotonic_or_raises(name):
    obj = _make_stub(_NON_MONOTONIC_CASES[name])
    try:
        obj._ensure_monotonic_dday()
    except ValueError:
        return
    assert np.all(np.diff(obj.dday) > 0), f"{name} left dday non-monotonic"


def test_forward_spike_repairs_the_spike_not_its_neighbor():
    # 0.9 is a forward clock spike; 0.3 and 0.4 are correct.  The repair must
    # move 0.9, not rewrite the perfectly good timestamps that follow it.
    dday = [0.0, 0.1, 0.9, 0.3, 0.4]
    obj = _make_stub(dday)
    obj._ensure_monotonic_dday()
    assert np.all(np.diff(obj.dday) > 0)
    assert obj.dday[3] == pytest.approx(0.3)
    assert obj.dday[4] == pytest.approx(0.4)
    assert obj.dday[2] == pytest.approx(0.2)
    assert len(obj.dday) == 5


def test_isolated_jump_then_segment_overlap():
    # Every backward jump must be classified on its own: an isolated bad ping
    # at index 2 and a genuine segment overlap starting at index 6.
    dday = [0.0, 0.1, 0.05, 0.2, 0.3, 0.4, 0.1, 0.11, 0.12, 0.13]
    obj = _make_stub(dday)
    obj._ensure_monotonic_dday()
    assert np.all(np.diff(obj.dday) > 0)
    np.testing.assert_allclose(obj.dday, [0.0, 0.1, 0.15, 0.2, 0.3, 0.4])
    assert len(obj.tsdat.pressure) == 6


def test_duplicate_timestamps():
    # Equal timestamps reach the isolated branch through ``diffs <= 0`` and
    # neighbor averaging cannot repair them.
    dday = [0.0, 0.1, 0.1, 0.1, 0.2, 0.3]
    obj = _make_stub(dday)
    obj._ensure_monotonic_dday()
    assert np.all(np.diff(obj.dday) > 0)
    # First occurrence of the repeated stamp is kept, the copies are spread out
    # towards the next larger stamp.
    assert obj.dday[0] == pytest.approx(0.0)
    assert obj.dday[1] == pytest.approx(0.1)
    assert obj.dday[4] == pytest.approx(0.2)
    assert obj.dday[5] == pytest.approx(0.3)
    assert len(obj.dday) == 6


def test_long_duplicate_run_raises():
    # A stalled clock cannot be repaired without fabricating data.
    dday = [0.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.3]
    obj = _make_stub(dday)
    with pytest.raises(ValueError, match="Cannot auto-fix"):
        obj._ensure_monotonic_dday()


def test_repair_propagates_to_tsdat_dday():
    """A repaired ping must not stay bad in `tsdat.dday`.

    `self.dday` is a drift-corrected *copy* of `self.tsdat.dday`, so repairing
    one in place leaves the other non-monotonic. Only the truncation branch
    used to touch both. `_lowpassfilter_pressure` builds its filter time base
    from `tsdat.dday`, so a record repaired here reached the filter with the
    original bad timestamps still in it.
    """
    obj = _make_stub([0.0, 0.1, 0.9, 0.3, 0.4])
    obj._ensure_monotonic_dday()

    assert np.all(np.diff(obj.dday) > 0)
    assert np.all(np.diff(obj.tsdat.dday) > 0), "tsdat.dday left non-monotonic"
    np.testing.assert_allclose(obj.tsdat.dday, obj.dday)


def test_repair_propagates_through_clock_drift():
    """With a drift correction active, `tsdat.dday` stays the raw time base."""
    obj = _make_stub([0.0, 0.1, 0.9, 0.3, 0.4], time_drift_rate=1.02)
    obj._ensure_monotonic_dday()

    assert np.all(np.diff(obj.tsdat.dday) > 0)
    # The defining relationship still holds after the repair: dday is the
    # corrected version of the raw record, not a copy of it.
    np.testing.assert_allclose(obj.dday, obj._correct_dday(obj.tsdat.dday))


def test_duplicate_repair_propagates_to_tsdat_dday():
    """The repeated-timestamp branch propagates too, not just the spike one."""
    obj = _make_stub([0.0, 0.1, 0.1, 0.1, 0.2, 0.3])
    obj._ensure_monotonic_dday()

    assert np.all(np.diff(obj.tsdat.dday) > 0)
    np.testing.assert_allclose(obj.tsdat.dday, obj.dday)


def test_untouched_pings_are_bit_identical_in_tsdat_dday():
    """Propagation must not perturb pings the repair did not touch."""
    raw = [0.0, 0.1, 0.9, 0.3, 0.4]
    obj = _make_stub(raw)
    obj._ensure_monotonic_dday()

    untouched = [0, 1, 3, 4]
    np.testing.assert_array_equal(
        obj.tsdat.dday[untouched], np.array(raw, dtype=float)[untouched]
    )


def test_truncation_keeps_both_arrays_aligned():
    """The truncation branch already sliced both; guard that it still does."""
    obj = _make_stub([0.0, 0.1, 0.2, 0.05, 0.06, 0.07, 0.08])
    obj._ensure_monotonic_dday()

    assert obj.dday.size == obj.tsdat.dday.size
    np.testing.assert_allclose(obj.tsdat.dday, obj.dday)
