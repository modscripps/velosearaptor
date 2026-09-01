"""The adaptive error velocity threshold must see the right samples (issue #100).

`_edit` sets the applied threshold to `min(max_e, max_e_deviation * std(e))`.
Two things put the wrong samples into that standard deviation:

1. `maskbins` was applied *after* `e.std()`, so bins the user declared bad set
   the threshold that decides which of the kept bins survive. When the masked
   bins are noisy their contribution pushes `max_e_deviation * std` past
   `max_e` and the adaptive criterion switches off entirely.

2. In `process_pings`, `_edit` ran once per `ens_size` chunk, so the chunk size
   was the window the standard deviation was estimated over. `ens_size` is
   documented as a memory knob; two people processing the same file with
   different values published different velocities.

The threshold is now computed once over the whole record. `process_pings`
applies the correlation and `maskbins` masks in the chunk loop and defers the
error velocity test until every ping has been read.
"""

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

FILE = "data/binmap_16670013.000"


@pytest.fixture
def adcpfile(rootdir):
    return rootdir / FILE


def _two_population_ensemble(kept_sigma, masked_sigma, nping=200, nbin=20):
    """An ensemble whose upper half is far noisier than its lower half.

    Correlation is well above any threshold, so only the error velocity test
    and `maskbins` can mask anything.
    """
    rng = np.random.default_rng(0)
    ens = Bunch()
    ens.cor = np.full((nping, nbin, 4), 128.0)
    xyze = np.zeros((nping, nbin, 4))
    half = nbin // 2
    xyze[:, :half, 3] = rng.normal(0.0, kept_sigma, (nping, half))
    xyze[:, half:, 3] = rng.normal(0.0, masked_sigma, (nping, nbin - half))
    # Every beam carries data, so `flag_no_data` is empty and the error
    # velocity is the only thing under test.
    ens.vel = np.ma.masked_array(np.zeros((nping, nbin, 4)), mask=False)
    ens.xyze = np.ma.masked_array(xyze, mask=False)
    return ens


def test_masked_bins_do_not_set_the_error_velocity_threshold(adcpfile):
    """The threshold follows the bins that are kept, not the ones masked out."""
    proc = ProcessADCP(adcpfile, META_DATA, magdec=0.0)
    nbin = 20
    half = nbin // 2
    maskbins = np.zeros(nbin, dtype=bool)
    maskbins[half:] = True
    proc.parse_editparams(
        {
            "max_e": 0.2,
            "max_e_deviation": 2,
            "min_correlation": 64,
            "maskbins": maskbins,
        }
    )

    ens = _two_population_ensemble(kept_sigma=0.01, masked_sigma=1.0, nbin=nbin)
    kept_e = ens.xyze[:, :half, 3].copy()
    proc._edit(ens)

    # The noisy half alone would drive 2 * std past max_e and switch the
    # adaptive branch off, leaving the fixed 0.2.
    assert ens.max_e_applied == pytest.approx(2 * np.std(kept_e), rel=1e-6)
    assert ens.max_e_applied < proc.editparams.max_e

    # The masked bins are still masked, and the threshold actually bit.
    assert np.ma.getmaskarray(ens.xyze)[:, half:, :].all()
    assert np.ma.getmaskarray(ens.xyze)[:, :half, :].any()


def test_fixed_threshold_still_wins_when_the_data_are_noisy(adcpfile):
    """Regression guard: the sigma branch must not bind on noisy data."""
    proc = ProcessADCP(adcpfile, META_DATA, magdec=0.0)
    proc.parse_editparams({"max_e": 0.2, "max_e_deviation": 2, "min_correlation": 64})

    ens = _two_population_ensemble(kept_sigma=1.0, masked_sigma=1.0)
    proc._edit(ens)

    assert ens.max_e_applied == proc.editparams.max_e


@pytest.mark.parametrize("binmap", [False, True])
def test_ens_size_does_not_change_process_pings_results(adcpfile, binmap):
    """`ens_size` is a memory knob and must not touch the numbers."""
    # The fixture is 2000 pings, so 250 splits it into eight chunks.
    ens_size = 250

    one_chunk = ProcessADCP(adcpfile, META_DATA, magdec=0.0)
    one_chunk.process_pings(binmap=binmap, ens_size=50000)

    many_chunks = ProcessADCP(adcpfile, META_DATA, magdec=0.0)
    many_chunks.process_pings(binmap=binmap, ens_size=ens_size)

    a, b = one_chunk.ds, many_chunks.ds
    # The split has to be real or this proves nothing.
    assert one_chunk.dday.size // ens_size > 1

    assert np.array_equal(
        a.max_e_applied.values, b.max_e_applied.values, equal_nan=True
    )
    for var in ["u", "v", "w", "e", "pg"]:
        assert np.array_equal(a[var].values, b[var].values, equal_nan=True), var


def test_max_e_applied_is_one_value_for_the_whole_record(adcpfile):
    """A record-wide threshold, not a per-chunk one."""
    proc = ProcessADCP(adcpfile, META_DATA, magdec=0.0)
    proc.process_pings(ens_size=250)

    applied = proc.ds.max_e_applied.values
    assert np.unique(applied[np.isfinite(applied)]).size == 1


def test_masking_before_rotation_equals_masking_after(adcpfile):
    """Deferring the error velocity test past `_to_enu` has to be a no-op.

    The threshold is applied to `xyze` today and to the rotated `enu` after
    this change, so the two have to agree: the error velocity must survive
    `rdi_xyz_enu` untouched, and masking one cell of `xyze` must equal masking
    all four components of the same cell of `enu`.
    """
    proc = ProcessADCP(adcpfile, META_DATA, magdec=0.0)
    ens = proc.m.read(start=0, stop=500)
    ens.dday = proc._correct_dday(ens.dday)
    proc._edit_masks(ens)

    xyze = ens.xyze.copy()
    proc._to_enu(ens)
    enu = ens.enu

    # The error velocity is not rotated.
    assert np.array_equal(
        np.ma.filled(xyze[..., 3], np.nan),
        np.ma.filled(enu[..., 3], np.nan),
        equal_nan=True,
    )

    over = np.abs(np.ma.filled(xyze[..., 3], 0.0)) > 0.13
    assert over.sum() > 0

    before = xyze.copy()
    before[over] = np.ma.masked
    ens_before = Bunch(ens)
    ens_before.xyze = before
    proc._to_enu(ens_before)

    after = enu.copy()
    after[over] = np.ma.masked

    assert np.array_equal(np.ma.getmaskarray(ens_before.enu), np.ma.getmaskarray(after))
    assert np.array_equal(
        np.ma.filled(ens_before.enu, np.nan),
        np.ma.filled(after, np.nan),
        equal_nan=True,
    )
