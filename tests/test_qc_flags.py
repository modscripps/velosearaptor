"""Per-criterion QC flags alongside the existing masking (issue #30).

Editing rejects a cell for one of three reasons, and the output records only
that it was rejected. Issue #30 wants the reason kept. This is the first step
of that: every criterion now produces a boolean flag, and the flags are
carried on `ens` next to the masks they already produce.

Nothing about the output changes, which is the property these tests exist to
hold. The load-bearing assertion is that the mask each editing step adds is
exactly the OR of the flags that step raised, on every configuration and every
call, so the flags are a decomposition of the current behavior and not a
second implementation of it. `tests/test_pinned_output.py` covers the other
half, that the published numbers are untouched.

The three flags, in the order they are applied:

`flag_cor`        correlation below `min_correlation` on any beam that enters
                  the velocity solution, reduced from the per-beam
                  `flag_cor_beam`
`flag_maskbins`   bins the user declared bad through `maskbins`
`flag_max_e`      error velocity above the applied threshold

The lifecycle differs by path and that is the easiest thing to get wrong here.
`average_ensembles` and `burst_average_ensembles` call `_edit`, so all three
flags are raised together, once per ensemble. `process_pings` never calls
`_edit`: it calls `_edit_masks` once per `ens_size` chunk and applies the error
velocity threshold once over the whole record after the loop (issue #100), so
`flag_cor` and `flag_maskbins` are per chunk while `flag_max_e` is per record.
Both call patterns are exercised below.
"""

import numpy as np
import pytest

from velosearaptor import madcp
from velosearaptor.madcp import ProcessADCP

META_DATA = {
    "mooring": "Test",
    "project": "Test",
    "lon": 0,
    "lat": 0,
}

UPLOOKER = "data/binmap_16670013.000"
BURST_FILE = "data/24606000.000"


def _proc(rootdir, filename=UPLOOKER, **kwargs):
    return ProcessADCP(rootdir / filename, META_DATA, magdec=0.0, **kwargs)


def _burst_proc(rootdir, **kwargs):
    return _proc(rootdir, BURST_FILE, tgridparams={"burst_average": True}, **kwargs)


def _cell_mask(xyze):
    """The mask of `xyze`, per cell.

    Every editing step masks all four components of a cell together, which
    `test_editing_masks_whole_cells` pins, so any component gives the same
    answer.
    """
    return np.ma.getmaskarray(xyze)[:, :, 0]


def _record_edit(monkeypatch, method):
    """Capture the cell mask before and after each call to `_edit`/`_edit_masks`.

    Returns the list the records land in. Each is the mask before, the mask
    after, and the ensemble, so a test can compare them against the flags the
    call left behind.
    """
    calls = []
    original = getattr(ProcessADCP, method)

    def recorder(self, ens):
        before = _cell_mask(ens.xyze).copy()
        original(self, ens)
        calls.append((before, _cell_mask(ens.xyze).copy(), ens))

    monkeypatch.setattr(ProcessADCP, method, recorder)
    return calls


def _flags(ens, names):
    return [np.asarray(ens[name]) for name in names]


def _or(flags):
    combined = np.zeros_like(flags[0], dtype=bool)
    for flag in flags:
        combined |= flag
    return combined


# ------------------------------------------------- the load-bearing assertion


@pytest.mark.parametrize("burst", [False, True])
def test_edit_flags_reproduce_the_mask_on_the_averaging_paths(
    rootdir, monkeypatch, burst
):
    """Per ensemble, `_edit` masks exactly the cells its three flags raise."""
    calls = _record_edit(monkeypatch, "_edit")
    if burst:
        proc = _burst_proc(rootdir)
        proc.burst_average_ensembles()
    else:
        proc = _proc(rootdir)
        proc.average_ensembles()

    assert calls, "no ensemble reached _edit"
    for before, after, ens in calls:
        flags = _flags(ens, ("flag_cor", "flag_maskbins", "flag_max_e"))
        for flag in flags:
            assert flag.dtype == bool
            assert flag.shape == before.shape
        assert np.array_equal(after, before | _or(flags))

    # A flag that never fires proves nothing about the decomposition.
    assert any(np.asarray(ens.flag_cor).any() for _, _, ens in calls)
    assert any(np.asarray(ens.flag_max_e).any() for _, _, ens in calls)


@pytest.mark.parametrize("binmap", [False, True])
def test_edit_masks_flags_reproduce_the_chunk_mask_on_the_per_ping_path(
    rootdir, monkeypatch, binmap
):
    """`process_pings` runs the two beam-space criteria once per chunk."""
    calls = _record_edit(monkeypatch, "_edit_masks")
    proc = _proc(rootdir)
    proc.process_pings(binmap=binmap)

    assert calls, "no chunk reached _edit_masks"
    for before, after, ens in calls:
        flags = _flags(ens, ("flag_cor", "flag_maskbins"))
        assert np.array_equal(after, before | _or(flags))
        # The threshold has not run yet on this path, so no third flag exists.
        assert "flag_max_e" not in ens
    assert any(np.asarray(ens.flag_cor).any() for _, _, ens in calls)


def test_the_record_threshold_flag_reproduces_the_mask_on_the_per_ping_path(
    rootdir, monkeypatch
):
    """`process_pings` applies the error velocity threshold after the loop.

    The flag is raised once, over the whole record, and the mask it adds to
    the error velocity is exactly that flag.
    """
    seen = {}
    original = madcp._max_e_flag

    def recorder(e, max_e_applied):
        flag = original(e, max_e_applied)
        seen["before"] = np.ma.getmaskarray(e).copy()
        seen["flag"] = flag.copy()
        return flag

    monkeypatch.setattr(madcp, "_max_e_flag", recorder)
    proc = _proc(rootdir)
    proc.process_pings()

    assert seen, "the record threshold never ran"
    assert seen["flag"].any()
    after = np.ma.getmaskarray(proc.ave.e)
    assert np.array_equal(after, seen["before"] | seen["flag"])


def test_editing_masks_whole_cells(rootdir, monkeypatch):
    """Every criterion masks all four components together.

    `_cell_mask` above reads component 0 alone and the flags are per cell, so
    the tests are only meaningful while this holds.
    """
    calls = _record_edit(monkeypatch, "_edit")
    proc = _proc(rootdir)
    proc.average_ensembles()

    for _, _, ens in calls:
        mask = np.ma.getmaskarray(ens.xyze)
        assert np.array_equal(mask, mask[:, :, :1].repeat(mask.shape[-1], axis=-1))


# ------------------------------------------------------ the individual flags


def test_the_maskbins_flag_marks_the_declared_bins(rootdir, monkeypatch):
    """`maskbins` is None by default, so nothing else here exercises it.

    Follows `test_adaptive_max_e.py`, which masks the same noisy bins of this
    file.
    """
    calls = _record_edit(monkeypatch, "_edit")
    proc = _proc(rootdir)
    masked_bins = [4, 5]
    proc.editparams.maskbins = proc.generate_binmask(masked_bins)
    proc.average_ensembles()

    assert calls
    for before, after, ens in calls:
        flag = np.asarray(ens.flag_maskbins)
        assert flag.any()
        # Exactly the declared bins, in every ping.
        assert np.array_equal(np.flatnonzero(flag.any(axis=0)), masked_bins)
        assert flag[:, masked_bins].all()
        flags = _flags(ens, ("flag_cor", "flag_maskbins", "flag_max_e"))
        assert np.array_equal(after, before | _or(flags))


def test_the_maskbins_flag_takes_bin_numbers_as_well_as_a_mask(rootdir, monkeypatch):
    """`_masked_bins` accepts both forms, and the flag follows it."""
    masked_bins = [0, 7]
    flags = {}
    for form in ("mask", "numbers"):
        calls = _record_edit(monkeypatch, "_edit")
        proc = _proc(rootdir)
        proc.editparams.maskbins = (
            proc.generate_binmask(masked_bins) if form == "mask" else masked_bins
        )
        proc.average_ensembles()
        flags[form] = np.asarray(calls[0][2].flag_maskbins)

    assert np.array_equal(np.flatnonzero(flags["mask"].any(axis=0)), masked_bins)
    assert np.array_equal(flags["mask"], flags["numbers"])


def test_the_per_beam_correlation_flag_reduces_to_the_cell_flag(rootdir, monkeypatch):
    """`flag_cor` is `flag_cor_beam` reduced over the beams that count.

    The per-beam array is what issue #18 (three-beam solutions) and the
    whole-cell rejection of issue #30 both need, and it is free here because
    the correlation test produces it before reducing.
    """
    calls = _record_edit(monkeypatch, "_edit")
    proc = _proc(rootdir)
    proc.average_ensembles()

    for _, _, ens in calls:
        beam = np.asarray(ens.flag_cor_beam)
        assert beam.dtype == bool
        assert beam.shape == ens.cor.shape
        assert np.array_equal(
            beam, np.asarray(ens.cor) < proc.editparams.min_correlation
        )
        assert np.array_equal(np.asarray(ens.flag_cor), beam.any(axis=-1))


def test_the_cell_correlation_flag_ignores_ibad_beams(rootdir, monkeypatch):
    """A beam dropped by `ibad` must not reject cells (issue #90).

    The per-beam flag still records what that beam did, which is the point of
    keeping the full beam axis.
    """
    ibad = 2
    calls = _record_edit(monkeypatch, "_edit")
    proc = _proc(rootdir, ibad=ibad)
    proc.average_ensembles()

    saw_a_difference = False
    for _, _, ens in calls:
        beam = np.asarray(ens.flag_cor_beam)
        assert beam.shape[-1] == np.asarray(ens.cor).shape[-1]
        kept = np.delete(beam, ibad, axis=-1)
        assert np.array_equal(np.asarray(ens.flag_cor), kept.any(axis=-1))
        saw_a_difference |= bool((beam.any(axis=-1) != kept.any(axis=-1)).any())

    # Otherwise the assertion above would hold for the wrong reason.
    assert saw_a_difference, "the excluded beam never differed from the rest"


def test_the_threshold_flag_is_the_error_velocity_test(rootdir, monkeypatch):
    """`flag_max_e` marks the surviving cells above the applied threshold."""
    calls = _record_edit(monkeypatch, "_edit")
    proc = _proc(rootdir)
    proc.average_ensembles()

    for _, _, ens in calls:
        flag = np.asarray(ens.flag_max_e)
        # Raised only where the beam-space criteria let the cell through.
        assert not (flag & np.asarray(ens.flag_cor)).any()
        assert not (flag & np.asarray(ens.flag_maskbins)).any()
        if np.isfinite(ens.max_e_applied):
            assert flag.any()


def test_the_flag_helpers_take_their_parameters_explicitly():
    """No `self` lookups, so the criteria are usable without a ProcessADCP.

    Cheap now and it is what a composable pipeline would need (issue #15).
    """
    cor = np.array([[[70, 70, 70, 10]]], dtype=np.uint8)
    beam = madcp._correlation_flag(cor, 64)
    assert np.array_equal(beam, [[[False, False, False, True]]])
    assert np.array_equal(madcp._correlation_cell_flag(beam), [[True]])
    assert np.array_equal(madcp._correlation_cell_flag(beam, ibad=3), [[False]])

    assert np.array_equal(madcp._maskbins_flag((2, 3), [1]), [[0, 1, 0], [0, 1, 0]])
    assert not madcp._maskbins_flag((2, 3), None).any()

    e = np.ma.masked_array([[0.1, 0.3, 0.5]], mask=[[False, False, True]])
    assert np.array_equal(madcp._max_e_flag(e, 0.2), [[False, True, False]])
    # A threshold that could not be computed rejects nothing.
    assert not madcp._max_e_flag(e, np.nan).any()
