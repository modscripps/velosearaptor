"""`t0`/`t1` given as strings must be converted against the file's `yearbase`.

`io.datetime64_to_yday0` returns the year of the requested time and drops it,
so a window in a later calendar year than the first ping used to come back as
a day-of-year of that later year and was then read against `yearbase`. See
upstream issue #97.

The bundled file `24606000.000` has `yearbase` 2021 and spans dday
187.72--191.51, i.e. 2021-07-07 17:15 to 2021-07-11 12:15.
"""

import pytest

from velosearaptor.madcp import ProcessADCP

META_DATA = {"mooring": "Test", "project": "Test", "lon": 0, "lat": 0}

# dday of the first and last ping in the test file.
DATA_DDAY0 = 187.71909722222222
DATA_DDAY1 = 191.5104513888889


@pytest.fixture
def proc(rootdir):
    """A processor built with default time gridding, so each test can call
    `parse_tgridparams` itself with the parameters it wants to check."""
    return ProcessADCP(rootdir / "data/24606000.000", META_DATA, magdec=0.0)


def test_cross_year_window_keeps_counting_past_365(proc):
    """A record straddling New Year: the ddays of a window in the second
    calendar year must continue past 365 rather than restart at the day of
    year of that second year."""
    # Shift the record forward so that it spans 2021-12-31 to 2022-01-04
    # while `yearbase` stays 2021.
    proc.dday = proc.dday + 178.0
    assert proc.yearbase == 2021

    proc.parse_tgridparams({"t0": "2022-01-02", "t1": "2022-01-03"})

    # 2021 is not a leap year, so 2022-01-01 is dday 365.
    assert proc.tgridparams.t0 == pytest.approx(366.0)
    assert proc.tgridparams.t1 == pytest.approx(367.0)
    assert proc.tgridparams.t0 > 365


def test_strings_within_the_file_year_are_unchanged(proc):
    """Regression guard: a window inside `yearbase` must convert exactly as
    it did before the fix."""
    proc.parse_tgridparams({"t0": "2021-07-09", "t1": "2021-07-10"})

    assert proc.tgridparams.t0 == pytest.approx(189.0)
    assert proc.tgridparams.t1 == pytest.approx(190.0)


def test_string_with_time_of_day(proc):
    """Sub-daily resolution survives the conversion."""
    proc.parse_tgridparams({"t0": "2021-07-09T06:00:00"})

    assert proc.tgridparams.t0 == pytest.approx(189.25)


def test_floats_pass_through_untouched(proc):
    """`t0`/`t1` given as ddays are not touched by the string conversion."""
    proc.parse_tgridparams({"t0": 188.5, "t1": 189.25, "dt_hours": 1.0})

    assert proc.tgridparams.t0 == 188.5
    assert proc.tgridparams.t1 == 189.25
    assert proc.tgridparams.dt_hours == 1.0


def test_window_without_overlap_raises(proc):
    """A window that misses the data entirely is an error, and the message
    names both ranges in calendar time."""
    with pytest.raises(ValueError) as excinfo:
        proc.parse_tgridparams({"t0": "2022-07-09", "t1": "2022-07-10"})

    msg = str(excinfo.value)
    # requested window
    assert "2022-07-09" in msg
    assert "2022-07-10" in msg
    # actual range of the file
    assert "2021-07-07" in msg
    assert "2021-07-11" in msg


def test_window_without_overlap_as_floats_raises(proc):
    """The check does not care how the window was spelled."""
    with pytest.raises(ValueError):
        proc.parse_tgridparams({"t0": DATA_DDAY1 + 1, "t1": DATA_DDAY1 + 2})


def test_partial_overlap_passes_quietly(proc):
    """A window reaching beyond the data on one side is a normal request."""
    proc.parse_tgridparams({"t0": "2021-07-05", "t1": "2021-07-10"})

    assert proc.tgridparams.t0 == pytest.approx(185.0)
    assert proc.tgridparams.t1 == pytest.approx(190.0)
