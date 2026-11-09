"""Tests for `driftparams["steps"]`, a clock that jumped once during the record
(motive-mooring-proc#94).

A step is `{"adcp": [Y, M, D, h, m, s], "seconds": s}`: the raw clock jumped by
`seconds` at raw time `adcp` (negative when it fell back). The correction adds
`-seconds` to every raw `dday` at or after `adcp`, then the single linear drift
rate is applied. The rate is computed from `end_adcp`/`end_pc` on the
step-corrected raw end time, so `end_adcp` keeps its meaning as the raw reading
at download.
"""

import numpy as np
import pytest
from pycurrents.codas import to_day

from velosearaptor.madcp import ProcessADCP

META_DATA = {"mooring": "Test", "project": "Test", "lon": 0, "lat": 0}

# 8000 pings at 3 s, raw dday 187.719 to 191.510 of 2021.
DOWNLOOKER = "24606000.000"
YEARBASE = 2021
# Download times after the last ping. The PC reads 100 s later than an
# instrument clock that fell back 100 s, so the corrected rate is exactly one
# and every expectation below is a plain shift.
END_ADCP = [2021, 7, 12, 0, 0, 0]
END_PC = [2021, 7, 12, 0, 1, 40]
STEP_S = -100.0
STEP_MID = [2021, 7, 9, 12, 0, 0]      # raw dday 189.5, inside the record
STEP_BEFORE = [2021, 7, 1, 0, 0, 0]    # before the first ping
STEP_AFTER = [2021, 7, 11, 23, 0, 0]   # after the last ping, before download

DAY = 86400.0


def _proc(rootdir, driftparams):
    return ProcessADCP(
        rootdir / "data" / DOWNLOOKER, META_DATA, magdec=0.0,
        driftparams=driftparams,
    )


def _drift(steps):
    return {"end_adcp": END_ADCP, "end_pc": END_PC,
            "steps": [{"adcp": s, "seconds": STEP_S} for s in steps]}


def test_step_in_the_middle_shifts_only_later_pings(rootdir):
    raw = _proc(rootdir, None).tsdat.dday.copy()
    p = _proc(rootdir, _drift([STEP_MID]))
    t_step = to_day(YEARBASE, *STEP_MID)
    after = raw >= t_step
    assert after.any() and (~after).any()
    assert p.clock_steps == [(pytest.approx(t_step), STEP_S)]
    np.testing.assert_allclose(p.time_drift_rate, 1.0, atol=1e-9)
    np.testing.assert_allclose(p.dday[~after], raw[~after], atol=1e-9)
    np.testing.assert_allclose(p.dday[after], raw[after] - STEP_S / DAY, atol=1e-9)
    # Still strictly increasing: a fallen-back clock is restored, not folded.
    assert np.all(np.diff(p.dday) > 0)


def test_step_before_the_first_ping_shifts_every_ping(rootdir):
    raw = _proc(rootdir, None).tsdat.dday.copy()
    p = _proc(rootdir, _drift([STEP_BEFORE]))
    np.testing.assert_allclose(p.dday, raw - STEP_S / DAY, atol=1e-9)


def test_step_after_the_last_ping_shifts_no_ping(rootdir):
    raw = _proc(rootdir, None).tsdat.dday.copy()
    p = _proc(rootdir, _drift([STEP_AFTER]))
    # The step still precedes the download, so the rate is computed on the
    # step-corrected end time and stays one.
    np.testing.assert_allclose(p.time_drift_rate, 1.0, atol=1e-9)
    np.testing.assert_allclose(p.dday, raw, atol=1e-9)


def test_uncorrect_round_trips_across_the_step(rootdir):
    p = _proc(rootdir, _drift([STEP_MID]))
    raw = p.tsdat.dday
    back = p._uncorrect_dday(p._correct_dday(raw))
    np.testing.assert_allclose(back, raw, atol=1e-9)
    # And the other way round, on both sides of the step.
    fwd = p._correct_dday(p._uncorrect_dday(p.dday))
    np.testing.assert_allclose(fwd, p.dday, atol=1e-9)


def test_rate_is_measured_on_the_step_corrected_end_time(rootdir):
    """Without the step the same download times imply a rate above one."""
    p_nostep = _proc(rootdir, {"end_adcp": END_ADCP, "end_pc": END_PC})
    assert p_nostep.time_drift_rate > 1.0 + 1e-6
    p = _proc(rootdir, _drift([STEP_MID]))
    np.testing.assert_allclose(p.time_drift_rate, 1.0, atol=1e-9)


def test_steps_without_download_times_apply_at_rate_one(rootdir):
    raw = _proc(rootdir, None).tsdat.dday.copy()
    p = _proc(rootdir, {"steps": [{"adcp": STEP_MID, "seconds": STEP_S}]})
    after = raw >= to_day(YEARBASE, *STEP_MID)
    assert p.time_drift_rate == 1
    np.testing.assert_allclose(p.dday[after], raw[after] - STEP_S / DAY, atol=1e-9)
    np.testing.assert_allclose(p.dday[~after], raw[~after], atol=1e-9)


def test_no_steps_leaves_the_correction_unchanged(rootdir):
    """The existing single-rate correction is bit-identical without `steps`."""
    plain = _proc(rootdir, {"end_adcp": END_ADCP, "end_pc": END_PC})
    empty = _proc(rootdir, {"end_adcp": END_ADCP, "end_pc": END_PC, "steps": []})
    assert plain.clock_steps == [] and empty.clock_steps == []
    assert np.array_equal(plain.dday, empty.dday)


def test_clock_steps_attribute_on_the_output(rootdir):
    p = _proc(rootdir, None)
    p.average_ensembles()
    assert p.ds.attrs["clock_steps"] == "none"
    p = _proc(rootdir, _drift([STEP_MID]))
    p.average_ensembles()
    assert p.ds.attrs["clock_steps"] == "2021-07-09T12:00:00 -100.0 s"
    assert p.ds.attrs["clock_correction"] == "applied"
