"""Test argument paths that no other test exercises (issue #101).

Four documented options crash instead of working. None of them corrupts data;
each one makes the option unusable, and each survived because nothing called it:

- ``process_pings(start=..., stop=...)`` raised ``UnboundLocalError`` because
  ``idx_start``/``idx_stop`` were bound only in the ``is None`` branches, and
  an empty range (``start == stop``) reached the same ``UnboundLocalError`` by
  running the chunk loop zero times;
- ``burst_average=True`` on continuous data produced
  ``pings_per_burst = -2147483648`` from a median over an empty selection and
  then died in an unrelated reduction;
- more than 32767 pings in one averaging interval raised ``OverflowError`` on
  the ``int16`` ``npings``, and would have wrapped ``ngood`` silently;
- ``interpolate_bin`` within two bins of either end of the profile indexed
  ``-1`` or past the last bin.
"""

import numpy as np
import pytest
from pycurrents.num import interp1

from velosearaptor.madcp import ProcessADCP

META_DATA = {
    "mooring": "Test",
    "project": "Test",
    "lon": 0,
    "lat": 0,
}

# tests/data/binmap_16670013.000 is continuous (2.5 s pings throughout).
# tests/data/24606000.000 is burst sampled (22 pings at 3 s, then an 837 s
# gap, 365 bursts) and is the only bundled file that reaches the burst path.
CONTINUOUS_FILE = "data/binmap_16670013.000"
BURST_FILE = "data/24606000.000"


@pytest.fixture
def continuous_proc(rootdir):
    return ProcessADCP(rootdir / CONTINUOUS_FILE, META_DATA, magdec=0.0)


@pytest.fixture
def burst_proc(rootdir):
    return ProcessADCP(
        rootdir / BURST_FILE,
        META_DATA,
        magdec=0.0,
        tgridparams={"burst_average": True},
    )


# --- 0.1 explicit start/stop on process_pings -------------------------------


def test_process_pings_start_and_stop(rootdir):
    """Ping indices given explicitly must select exactly that ping range."""
    proc = ProcessADCP(rootdir / CONTINUOUS_FILE, META_DATA, magdec=0.0)
    start, stop = 10, 60
    proc.process_pings(start=start, stop=stop)
    assert proc.ds.time.size == stop - start
    assert proc.ave.pg.shape[0] == stop - start


def test_process_pings_stop_only(rootdir):
    """`stop` alone must work; `start` falls back to the first ping."""
    proc = ProcessADCP(rootdir / CONTINUOUS_FILE, META_DATA, magdec=0.0)
    proc.process_pings(stop=20)
    assert proc.ds.time.size == 20


def test_process_pings_start_only(rootdir):
    """`start` alone must work; `stop` falls back to the last ping."""
    proc = ProcessADCP(rootdir / CONTINUOUS_FILE, META_DATA, magdec=0.0)
    proc.process_pings(start=10)
    idx_stop = np.searchsorted(proc.dday, proc.dday_end)
    assert proc.ds.time.size == idx_stop - 10


def test_process_pings_start_stop_are_ping_indices(rootdir):
    """The selected pings are the ones at those indices in `dday`."""
    proc = ProcessADCP(rootdir / CONTINUOUS_FILE, META_DATA, magdec=0.0)
    proc.process_pings(start=10, stop=30)
    expected = proc.dday[10:30]
    got = proc.ave.dday
    assert np.allclose(np.asarray(got), np.asarray(expected))


def test_process_pings_rejects_out_of_range_indices(rootdir):
    """A ping index outside the record is named, not turned into a crash.

    Leaving these unchecked is the same class of problem as the
    `UnboundLocalError` this item fixes: a documented option fails somewhere
    far from the argument that caused it.
    """
    proc = ProcessADCP(rootdir / CONTINUOUS_FILE, META_DATA, magdec=0.0)
    npings = proc.dday.size

    with pytest.raises(ValueError, match="start"):
        proc.process_pings(start=npings + 1)
    with pytest.raises(ValueError, match="stop"):
        proc.process_pings(stop=npings + 1)


def test_process_pings_rejects_negative_indices(rootdir):
    """Negative indices are not wrapped Python-style; they are refused.

    `np.searchsorted` never returns a negative index, so the `None` path
    cannot produce one; accepting one here would silently select a different
    ping range from the one the caller asked for.
    """
    proc = ProcessADCP(rootdir / CONTINUOUS_FILE, META_DATA, magdec=0.0)

    with pytest.raises(ValueError, match="negative"):
        proc.process_pings(start=-5)
    with pytest.raises(ValueError, match="negative"):
        proc.process_pings(stop=-5)


def test_process_pings_rejects_start_after_stop(rootdir):
    """`start > stop` gives a negative ping count; say so at the argument."""
    proc = ProcessADCP(rootdir / CONTINUOUS_FILE, META_DATA, magdec=0.0)

    with pytest.raises(ValueError, match="after"):
        proc.process_pings(start=60, stop=10)


def test_process_pings_rejects_an_empty_range(rootdir):
    """`start == stop` selects no pings and must raise at the argument.

    It reached the chunk loop, which never ran, and then died on the
    `UnboundLocalError` for `ens` that this issue is about in the first place.
    """
    proc = ProcessADCP(rootdir / CONTINUOUS_FILE, META_DATA, magdec=0.0)

    for index in (0, 5, proc.dday.size - 1):
        with pytest.raises(ValueError, match="select no pings"):
            proc.process_pings(start=index, stop=index)


def test_process_pings_accepts_the_full_record(rootdir):
    """The bounds are inclusive of the record end, as Python slicing is."""
    proc = ProcessADCP(rootdir / CONTINUOUS_FILE, META_DATA, magdec=0.0)
    npings = proc.dday.size
    proc.process_pings(start=0, stop=npings)
    assert proc.ds.time.size == npings


# --- 0.2 burst_average on continuous data -----------------------------------


def test_burst_average_on_continuous_data_raises(rootdir):
    """Continuous data must be rejected with a message a user can act on."""
    with pytest.raises(ValueError, match="not burst sampled"):
        ProcessADCP(
            rootdir / CONTINUOUS_FILE,
            META_DATA,
            magdec=0.0,
            tgridparams={"burst_average": True},
        )


def test_burst_average_on_burst_data_still_works(burst_proc):
    """Regression guard: the genuinely burst-sampled file is not rejected."""
    assert burst_proc.is_burst_average
    assert burst_proc.start_ddays.size == 365
    assert burst_proc.dt > 0


# --- 0.3 ping counts beyond the int16 range ---------------------------------

# 12 h or 24 h averages at a 1 s ping rate put more than 32767 pings in one
# averaging interval. Synthesizing such a file would take an hour of test time,
# so assert on the arrays the averaging methods actually allocate.
BIG_COUNT = 40000


@pytest.fixture
def averaged(rootdir):
    proc = ProcessADCP(
        rootdir / CONTINUOUS_FILE, META_DATA, magdec=0.0, tgridparams={"dt_hours": 0.25}
    )
    proc.average_ensembles()
    return proc


@pytest.fixture
def burst_averaged(burst_proc):
    burst_proc.burst_average_ensembles(stop=5)
    return burst_proc


@pytest.mark.parametrize("path", ["averaged", "burst_averaged"])
def test_npings_holds_more_than_int16(path, request):
    proc = request.getfixturevalue(path)
    npings = proc.ave.npings
    npings[0] = BIG_COUNT
    assert npings[0] == BIG_COUNT


@pytest.mark.parametrize("path", ["averaged", "burst_averaged"])
def test_ngood_does_not_wrap(path, request):
    proc = request.getfixturevalue(path)
    ngood = proc.ave.ngood
    counts = np.array([40000, 33000, 20000])
    ngood[0, :3] = counts
    assert np.array_equal(np.asarray(ngood[0, :3]), counts)


@pytest.mark.parametrize("path", ["averaged", "burst_averaged"])
def test_pg_from_large_counts_stays_in_range(path, request):
    """The percent good computed from a >int16 count must stay in 0-100."""
    proc = request.getfixturevalue(path)
    ngood = proc.ave.ngood
    nprofs = BIG_COUNT
    ngood[0, :3] = np.array([40000, 33000, 20000])
    pgi = 100 * ngood[0, :3].astype(np.int64) // nprofs
    assert np.array_equal(np.asarray(pgi), np.array([100, 82, 50]))


# --- 0.4 interpolate_bin near the ends of the profile -----------------------


def _nbins(proc):
    return proc.tsdat.dep.size


@pytest.mark.parametrize("offset", [1, 2])
def test_interpolate_bin_near_ends(rootdir, offset):
    """Bins two from either end are the plausible targets and must work."""
    proc = ProcessADCP(
        rootdir / BURST_FILE,
        META_DATA,
        magdec=0.0,
        tgridparams={"burst_average": True},
    )
    nbins = _nbins(proc)
    for zi in (offset, nbins - 1 - offset):
        proc.burst_average_ensembles(stop=3, interpolate_bin=zi)
        assert np.isfinite(proc.ave.u).any()


def test_interpolate_bin_at_the_very_edge_raises(rootdir):
    """The first and last bin have no neighbour on one side."""
    proc = ProcessADCP(
        rootdir / BURST_FILE,
        META_DATA,
        magdec=0.0,
        tgridparams={"burst_average": True},
    )
    nbins = _nbins(proc)
    for zi in (0, nbins - 1):
        with pytest.raises(ValueError, match="interpolate_bin"):
            proc.burst_average_ensembles(stop=3, interpolate_bin=zi)


def test_interpolate_bin_out_of_range_raises(rootdir):
    proc = ProcessADCP(
        rootdir / BURST_FILE,
        META_DATA,
        magdec=0.0,
        tgridparams={"burst_average": True},
    )
    with pytest.raises(ValueError, match="interpolate_bin"):
        proc.burst_average_ensembles(stop=3, interpolate_bin=_nbins(proc))


def test_interpolated_bin_keeps_pg(rootdir):
    """Decision: the interpolated bin keeps its own pg, so the interpolation
    stays visible in the product and is not mistaken for measured data."""
    kw = {"magdec": 0.0, "tgridparams": {"burst_average": True}}
    plain = ProcessADCP(rootdir / BURST_FILE, META_DATA, **kw)
    plain.burst_average_ensembles(stop=3)
    interp = ProcessADCP(rootdir / BURST_FILE, META_DATA, **kw)
    interp.burst_average_ensembles(stop=3, interpolate_bin=2)
    # `pg` and `ngood` carry NaN off the instrument's profile (issue #82), so
    # the comparisons have to treat NaN as equal to NaN.
    assert np.array_equal(plain.ave.pg, interp.ave.pg, equal_nan=True)
    assert np.array_equal(
        np.asarray(plain.ave.ngood), np.asarray(interp.ave.ngood), equal_nan=True
    )


@pytest.mark.parametrize("orientation", ["down", "up"])
def test_neighbor_window_both_orientations(orientation):
    """The clipped window must bracket the target bin and stay interpolable
    for an ascending (downlooker) and a descending (uplooker) depth vector.

    Only the downlooking `24606000.000` is burst sampled among the bundled
    files, so the uplooking depth ordering is covered here rather than end to
    end.
    """
    nbins = 19
    dep = np.arange(nbins) * 8.0 + 12.0
    depth = 500.0 + dep if orientation == "down" else 500.0 - dep
    for zi in (1, nbins - 2):
        neighbors = ProcessADCP._interpolation_neighbors(zi, nbins)
        assert min(neighbors) >= 0
        assert max(neighbors) < nbins
        assert zi not in neighbors
        assert min(neighbors) < zi < max(neighbors)
        values = np.arange(nbins, dtype=float)[:, np.newaxis]
        out = interp1(
            depth[neighbors], values[neighbors], depth[zi], axis=0, method="linear"
        )
        assert np.isfinite(np.ma.filled(out, np.nan)).all()
