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

The flags, in the order they are applied:

`flag_no_data`    no beam that enters the velocity solution carries a velocity
                  at this cell, so no solution exists before editing looks at
                  it; reduced from the per-beam `flag_no_data_beam`
`flag_cor`        correlation below `min_correlation` on any beam that enters
                  the velocity solution, reduced from the per-beam
                  `flag_cor_beam`
`flag_maskbins`   bins the user declared bad through `maskbins`
`flag_max_e`      error velocity above the applied threshold

`flag_no_data` is the one editing does not cause. Binmapping cannot fill every
cell of every beam and the instrument rejects beams on its own, so cells reach
`_edit` already masked, and the three editing flags alone do not add up to the
mask on `xyze`. Raising it makes `ens.valid`, the complement of the four,
assertable against that mask, which is what lets percent good be counted from
the flags instead of read back out of the damaged array (issue #30).

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
from pycurrents.num import interp1

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


CONFIGURATIONS = (
    "average_ensembles",
    "burst_average_ensembles",
    "process_pings",
    "process_pings_binmap",
)


def _run(rootdir, monkeypatch, config, method="_edit_masks", **kwargs):
    """Run one of the four pinned configurations, recording calls to `method`.

    `_edit` calls `_edit_masks`, so recording the latter reaches every
    configuration; recording the former reaches only the two averaging paths.
    """
    calls = _record_edit(monkeypatch, method)
    if config == "burst_average_ensembles":
        proc = _burst_proc(rootdir, **kwargs)
        proc.burst_average_ensembles()
    elif config == "average_ensembles":
        proc = _proc(rootdir, **kwargs)
        proc.average_ensembles()
    else:
        proc = _proc(rootdir, **kwargs)
        proc.process_pings(binmap=config.endswith("binmap"))
    assert calls, f"no ensemble reached {method}"
    return proc, calls


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

        # The closed form: `before` is `flag_no_data`, so the four flags
        # account for the whole mask and `ens.valid` is their complement.
        all_flags = _flags(
            ens, ("flag_no_data", "flag_cor", "flag_maskbins", "flag_max_e")
        )
        assert np.array_equal(after, _or(all_flags))
        assert np.asarray(ens.valid).dtype == bool
        assert np.array_equal(np.asarray(ens.valid), ~after)

    # A flag that never fires proves nothing about the decomposition.
    assert any(np.asarray(ens.flag_cor).any() for _, _, ens in calls)
    assert any(np.asarray(ens.flag_max_e).any() for _, _, ens in calls)
    assert any(np.asarray(ens.flag_no_data).any() for _, _, ens in calls)


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
        # The threshold has not run yet on this path, so no fourth flag
        # exists and `ens.valid` carries only the three raised so far.
        assert "flag_max_e" not in ens
        all_flags = _flags(ens, ("flag_no_data", "flag_cor", "flag_maskbins"))
        assert np.array_equal(after, _or(all_flags))
        assert np.array_equal(np.asarray(ens.valid), ~after)
    assert any(np.asarray(ens.flag_cor).any() for _, _, ens in calls)
    assert any(np.asarray(ens.flag_no_data).any() for _, _, ens in calls)


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


@pytest.mark.parametrize("config", CONFIGURATIONS)
def test_to_enu_masks_every_component_from_valid(rootdir, monkeypatch, config):
    """`_to_enu` leaves `enu` masked exactly where `ens.valid` is False.

    All four components together, so `_cell_mask` above can read component
    0 alone. On the per-ping path `valid` carries the three beam-space
    criteria at this point and the error velocity flag is applied to the
    whole record after the loop, so the equality holds there too.

    `rdi_xyz_enu` also masks `u` and `v` of any ping whose heading, pitch or
    roll is masked, which `valid` does not record. Neither bundled file has
    such a ping, so the equality is exact here.
    """
    seen = []
    original = ProcessADCP._to_enu

    def recorder(self, ens):
        original(self, ens)
        mask = np.ma.getmaskarray(ens.enu)
        valid = np.asarray(ens.valid)
        seen.append(all(np.array_equal(mask[..., k], ~valid) for k in range(4)))
        seen.append((~valid).any())

    monkeypatch.setattr(ProcessADCP, "_to_enu", recorder)
    if config == "burst_average_ensembles":
        _burst_proc(rootdir).burst_average_ensembles()
    elif config == "average_ensembles":
        _proc(rootdir).average_ensembles()
    else:
        _proc(rootdir).process_pings(binmap=config.endswith("binmap"))

    assert seen, "no ensemble reached _to_enu"
    assert all(seen[0::2])
    # Otherwise the equality holds because nothing was rejected.
    assert any(seen[1::2])


def test_to_enu_under_ibad_masks_e_entirely_and_the_rest_from_valid(
    rootdir, monkeypatch
):
    """With `ibad`, `beam_to_xyz` masks the error velocity everywhere.

    A three-beam solution has no error velocity, so `e` arrives fully masked
    and `_to_enu` can add nothing there. `u`, `v` and `w` are masked from
    `valid` as without `ibad`. `ibad` is outside the four pinned
    configurations, so nothing else here would see this.
    """
    seen = []
    original = ProcessADCP._to_enu

    def recorder(self, ens):
        original(self, ens)
        mask = np.ma.getmaskarray(ens.enu)
        valid = np.asarray(ens.valid)
        seen.append(
            all(np.array_equal(mask[..., k], ~valid) for k in range(3))
            and bool(mask[..., 3].all())
        )

    monkeypatch.setattr(ProcessADCP, "_to_enu", recorder)
    _proc(rootdir, ibad=0).average_ensembles()

    assert seen, "no ensemble reached _to_enu"
    assert all(seen)


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


@pytest.mark.parametrize("config", CONFIGURATIONS)
def test_the_no_data_flag_is_the_invalidity_editing_inherits(
    rootdir, monkeypatch, config
):
    """`flag_no_data` is the mask `xyze` already carries when editing starts.

    Editing does not cause it. Binmapping cannot fill every cell of every
    beam and the instrument rejects beams on its own, so the flag is a
    snapshot of what arrives rather than a criterion applied here. The
    equality against the entry mask is what makes `ens.valid` a decomposition
    of the mask instead of a second opinion about it.
    """
    _, calls = _run(rootdir, monkeypatch, config)

    saw_one = False
    for before, _, ens in calls:
        beam = np.asarray(ens.flag_no_data_beam)
        assert beam.dtype == bool
        assert beam.shape == ens.vel.shape
        assert np.array_equal(beam, np.ma.getmaskarray(ens.vel))
        assert np.array_equal(np.asarray(ens.flag_no_data), beam.any(axis=-1))
        assert np.array_equal(np.asarray(ens.flag_no_data), before)
        saw_one |= bool(beam.any())

    # Otherwise the equality would hold because both sides are empty.
    assert saw_one, "no cell reached editing already masked"


def test_the_no_data_cell_flag_ignores_ibad_beams(rootdir, monkeypatch):
    """A beam dropped by `ibad` must not declare a cell unsolvable (issue #90).

    `Transform.beam_to_xyz` fills the excluded beam from the other three
    before it transforms, so that beam's mask never reaches `xyze`. Reducing
    over all four beams instead disagrees with the mask, and `ibad` is outside
    the four pinned configurations, so nothing else here would catch it.
    """
    ibad = 0
    _, calls = _run(rootdir, monkeypatch, "average_ensembles", ibad=ibad)

    saw_a_difference = False
    for before, _, ens in calls:
        beam = np.asarray(ens.flag_no_data_beam)
        kept = np.delete(beam, ibad, axis=-1)
        assert np.array_equal(np.asarray(ens.flag_no_data), kept.any(axis=-1))
        assert np.array_equal(np.asarray(ens.flag_no_data), before)
        saw_a_difference |= bool((beam.any(axis=-1) != kept.any(axis=-1)).any())

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
    # One reduction, shared by every beam-space criterion.
    assert np.array_equal(madcp._cell_flag(beam), [[True]])
    assert np.array_equal(madcp._cell_flag(beam, ibad=3), [[False]])

    vel = np.ma.masked_array(np.zeros((1, 1, 4)), mask=[[[False, False, False, True]]])
    assert np.array_equal(madcp._no_data_flag(vel), beam)
    assert not madcp._no_data_flag(np.zeros((1, 1, 4))).any()

    assert np.array_equal(madcp._maskbins_flag((2, 3), [1]), [[0, 1, 0], [0, 1, 0]])
    assert not madcp._maskbins_flag((2, 3), None).any()

    e = np.ma.masked_array([[0.1, 0.3, 0.5]], mask=[[False, False, True]])
    assert np.array_equal(madcp._max_e_flag(e, 0.2), [[False, True, False]])
    # A threshold that could not be computed rejects nothing.
    assert not madcp._max_e_flag(e, np.nan).any()


def test_apply_qc_masks_every_component_of_an_invalid_cell():
    """The one write the QC criteria make, in place and idempotent."""
    enu = np.ma.asarray(np.ones((2, 3, 4)))
    valid = np.array([[True, False, True], [True, True, False]])
    madcp._apply_qc(enu, valid)
    mask = np.ma.getmaskarray(enu)
    for k in range(4):
        assert np.array_equal(mask[..., k], ~valid)

    # Narrowing `valid` and applying again adds the new cells and keeps the
    # old ones, which is how `process_pings` applies the record-wide flag.
    narrower = valid & np.array([[True, True, False], [True, True, True]])
    madcp._apply_qc(enu, narrower)
    mask = np.ma.getmaskarray(enu)
    for k in range(4):
        assert np.array_equal(mask[..., k], ~narrower)


# ------------------------------------------------- percent good from the flags


def test_process_pings_percent_good_is_the_assembled_validity(rootdir, monkeypatch):
    """On the per-ping path `pg` is validity, assembled in two stages.

    The beam-space criteria are raised once per `ens_size` chunk and the error
    velocity threshold once over the whole record (issue #100), so the chunks
    have to be stitched back together before the record-wide flag is applied.
    `ens_size` is small here to make the chunking real.
    """
    seen = {}
    original = madcp._max_e_flag

    def recorder(e, max_e_applied):
        flag = original(e, max_e_applied)
        seen["flag"] = flag.copy()
        return flag

    monkeypatch.setattr(madcp, "_max_e_flag", recorder)
    calls = _record_edit(monkeypatch, "_edit_masks")
    proc = _proc(rootdir)
    proc.process_pings(ens_size=500)

    assert len(calls) > 1, "the record was not chunked"
    valid = np.concatenate([np.asarray(ens.valid) for _, _, ens in calls], axis=0)
    valid &= ~seen["flag"]
    assert np.array_equal(proc.ave.pg, 100 * valid.astype(np.int8))
    assert valid.any() and not valid.all()


def test_burst_percent_good_is_the_validity_count(rootdir, monkeypatch):
    """The burst path counts on the instrument bins, so validity counts there.

    Every ping of a burst sits on one tiled depth vector, so the count is
    exact in bin space and the gridding of the count is unchanged (issue #30,
    decision of 2026-08-31).
    """
    _, calls = _run(rootdir, monkeypatch, "burst_average_ensembles", method="_edit")

    counted = False
    for _, _, ens in calls:
        valid = np.asarray(ens.valid)
        assert np.array_equal(valid.sum(axis=0), ens.enu[..., 0].count(axis=0))
        counted |= bool(valid.any())
    assert counted


def test_average_ensembles_percent_good_is_the_gridded_validity(rootdir, monkeypatch):
    """The spaced path counts on the depth grid, so validity is gridded too.

    The transducer moves within an averaging interval there, so the pings of
    one interval share no bin grid and the count has to be taken after
    gridding. Validity rides as one more column of the single `interp1` call
    in `_regrid_enu_amp`, so it costs no extra interpolation.
    """
    _, calls = _run(rootdir, monkeypatch, "average_ensembles", method="_edit")

    for _, _, ens in calls:
        valid_grid = np.asarray(ens.valid_grid)
        assert valid_grid.dtype == bool
        assert np.array_equal(valid_grid, ~np.isnan(ens.enu_grid[..., 0]))


def test_the_extra_columns_leave_the_regridding_alone(rootdir):
    """Ten columns must give exactly what four plus amplitude gave.

    `_regrid_enu_amp` was deliberately reduced to one `interp1` call per
    ping, and the published velocities and amplitudes come out of it. The
    validity column and the four reason columns ride along in that same
    call, so this rebuilds the call without any of them and compares. It is
    the only guard on the one place counting from flags touches a published
    array.
    """
    proc = _proc(rootdir)
    ens = proc.read_ensemble(0)
    proc._edit(ens)
    proc._to_enu(ens)
    proc._regrid_enu_amp(ens)

    depth = proc._burst_average_depth(ens)
    ncols = ens.enu.shape[-1]
    for i in range(ens.dday.size):
        amp_col = ens.amp[i].mean(axis=-1, keepdims=True)
        combined = np.ma.concatenate([ens.enu[i], amp_col], axis=-1)
        without = np.ma.filled(
            interp1(depth[i], combined, proc.dgrid, axis=0, method="linear"), np.nan
        )
        assert np.array_equal(ens.enu_grid[i], without[:, :ncols], equal_nan=True)
        assert np.array_equal(ens.amp_grid[i], without[:, ncols], equal_nan=True)

    assert np.array_equal(ens.valid_grid, ~np.isnan(ens.enu_grid[..., 0]))


def test_an_unreadable_chunk_is_zero_percent_good(rootdir, monkeypatch):
    """Pings the reader cannot return never reach `_edit_masks`.

    Validity assembled from flags has no entry for them, so the array it
    assembles into starts False and they come out 0 % good. Reading `pg` back
    out of `uvwe`, which is what this replaces, got them right for free. None
    of the pinned configurations contains an unreadable chunk, so the reader
    is made to fail one here.
    """
    proc = _proc(rootdir)
    reads = []
    original = proc.m.read

    def read(start=None, stop=None, **kwargs):
        reads.append((start, stop))
        if len(reads) == 2:
            return None
        return original(start=start, stop=stop, **kwargs)

    monkeypatch.setattr(proc.m, "read", read)
    proc.process_pings(ens_size=500)

    assert len(reads) > 2, "the dropped chunk was the last one"
    offset = reads[0][0]
    idx0, idx1 = reads[1][0] - offset, reads[1][1] - offset
    pg = proc.ave.pg
    assert (pg[idx0:idx1] == 0).all()
    assert (pg[:idx0] == 100).any()
    assert (pg[idx1:] == 100).any()


# --------------------------------------------- one reason per rejected cell


def _reasons(ens):
    return [np.asarray(ens[f"reason_{name}"]) for name in madcp.QC_CRITERIA]


@pytest.mark.parametrize("config", ["average_ensembles", "burst_average_ensembles"])
def test_the_reasons_partition_the_rejected_cells(rootdir, monkeypatch, config):
    """Every rejected cell is attributed to exactly one criterion.

    This is what lets the published counts sum to `nprofs - ngood`. The raw
    flags overlap, so publishing them directly would double-count.
    """
    _, calls = _run(rootdir, monkeypatch, config, method="_edit")

    for _, after, ens in calls:
        reasons = _reasons(ens)
        stacked = np.stack(reasons, axis=-1)
        # Disjoint: no cell carries two reasons.
        assert stacked.sum(axis=-1).max() <= 1
        # Complete: their union is exactly the mask editing produced.
        assert np.array_equal(_or(reasons), after)
        assert np.array_equal(_or(reasons), ~np.asarray(ens.valid))


@pytest.mark.parametrize("config", CONFIGURATIONS)
def test_the_correlation_test_fires_on_cells_that_carry_no_data(
    rootdir, monkeypatch, config
):
    """The exclusion in the precedence rule is not defensive.

    `_correlation_flag` compares the fill data underneath already-masked
    cells, and `_mask_binmapped` writes 0 there, which is below any
    threshold. Without this, the partition test above would hold for the
    wrong reason. The burst file is the one configuration where it does not
    happen, so it is excluded from the assertion.
    """
    _, calls = _run(rootdir, monkeypatch, config)

    overlap = sum(
        int((np.asarray(ens.flag_cor) & np.asarray(ens.flag_no_data)).sum())
        for _, _, ens in calls
    )
    if config == "burst_average_ensembles":
        assert overlap == 0
    else:
        assert overlap > 0
        # And the precedence rule removes it from the correlation reason.
        for _, _, ens in calls:
            assert not (np.asarray(ens.reason_cor) & np.asarray(ens.flag_no_data)).any()


def test_the_attribution_helper_takes_its_parameters_explicitly():
    """No `self` and no ensemble, so the rule is testable on its own."""
    nodata = np.array([[True, False, False, False]])
    cor = np.array([[True, True, False, False]])
    maskbins = np.array([[True, True, True, False]])

    a, c, m = madcp._attributed_flags(nodata, cor, maskbins)
    assert np.array_equal(a, [[True, False, False, False]])
    assert np.array_equal(c, [[False, True, False, False]])
    assert np.array_equal(m, [[False, False, True, False]])
    # Disjoint and complete over the union of the inputs.
    assert np.array_equal(a | c | m, nodata | cor | maskbins)
    assert (a.astype(int) + c.astype(int) + m.astype(int)).max() == 1


# ------------------------------------------------ the published nbad counts


NBAD = tuple(f"nbad_{name}" for name in madcp.QC_CRITERIA)


def test_process_pings_publishes_the_reason_for_every_rejected_ping(
    rootdir, monkeypatch
):
    """`pg` says how many pings survived; these say why the rest did not."""
    seen = {}
    original = madcp._max_e_flag

    def recorder(e, max_e_applied):
        flag = original(e, max_e_applied)
        seen["flag"] = flag.copy()
        return flag

    monkeypatch.setattr(madcp, "_max_e_flag", recorder)
    calls = _record_edit(monkeypatch, "_edit_masks")
    proc = _proc(rootdir)
    proc.process_pings(ens_size=500)

    assert len(calls) > 1, "the record was not chunked"
    for name in ("nodata", "cor", "maskbins"):
        expected = np.concatenate(
            [np.asarray(ens[f"reason_{name}"]) for _, _, ens in calls], axis=0
        )
        assert np.array_equal(proc.ave[f"nbad_{name}"], expected.astype(np.int8))
    assert np.array_equal(proc.ave.nbad_max_e, seen["flag"].astype(np.int8))


def test_process_pings_counts_partition_the_rejected_pings(rootdir):
    """Exactly: `pg` is 100 where good, and one count is 1 where not."""
    proc = _proc(rootdir)
    proc.process_pings()

    total = sum(proc.ave[name].astype(np.int32) for name in NBAD)
    good = proc.ave.pg == 100
    assert np.array_equal(total, (~good).astype(np.int32))
    assert total.max() == 1
    for name in NBAD:
        assert proc.ave[name].dtype == np.int8


def test_process_pings_counts_stay_int8_and_unmasked(rootdir):
    """`z` is the instrument bin axis, so there is no off-profile cell here.

    Masking these the way `pg` is masked would promote them to float64 and
    more than double the largest product this package writes.
    """
    proc = _proc(rootdir)
    proc.process_pings()

    for name in NBAD:
        assert proc.ds[name].dtype == np.int8
        assert proc.ds[name].dims == ("z", "time")
        # `int8` cannot carry NaN, so the interesting failure mode is the
        # masking `_ave2nc` applies to the two averaging paths reaching this
        # one by mistake, which would wrap the values in a masked array.
        assert not np.ma.isMaskedArray(proc.ds[name].values)


def test_an_unreadable_chunk_has_no_rejection_reason(rootdir, monkeypatch):
    """Nothing was read, so no criterion rejected anything.

    `pg` is already 0 there and `u` is NaN, so the cell is not ambiguous.
    """
    proc = _proc(rootdir)
    # Declare bad bins so `nbad_maskbins` fires on every ping that was read.
    # Without this it is zero record-wide and the guard below would be
    # vacuous for it, which is the one criterion this test most needs to
    # bite on: `maskbins` rejects every ping it touches, so a dropped chunk
    # reading 0 is only meaningful when the rest of the record reads 1.
    proc.editparams.maskbins = proc.generate_binmask([4, 5])
    reads = []
    original = proc.m.read

    def read(start=None, stop=None, **kwargs):
        reads.append((start, stop))
        if len(reads) == 2:
            return None
        return original(start=start, stop=stop, **kwargs)

    monkeypatch.setattr(proc.m, "read", read)
    proc.process_pings(ens_size=500)

    assert len(reads) > 2, "the dropped chunk was the last one"
    offset = reads[0][0]
    idx0, idx1 = reads[1][0] - offset, reads[1][1] - offset
    for name in NBAD:
        assert (proc.ave[name][idx0:idx1] == 0).all()
        assert proc.ave[name][:idx0].any() or proc.ave[name][idx1:].any()


def test_burst_counts_partition_the_rejected_pings(rootdir, monkeypatch):
    """Counted on the instrument bins, where the partition is exact.

    Checked before the interpolation onto the depth grid, because `interp1`
    on a count is what smears it.
    """
    _, calls = _run(rootdir, monkeypatch, "burst_average_ensembles", method="_edit")

    for _, _, ens in calls:
        counts = [
            np.sum(np.asarray(ens[f"reason_{name}"]), axis=0)
            for name in madcp.QC_CRITERIA
        ]
        nprofs = ens.enu.shape[0]
        ngood = np.sum(np.asarray(ens.valid), axis=0)
        assert np.array_equal(sum(counts), nprofs - ngood)


def test_burst_counts_are_nan_off_the_profile(rootdir):
    """A grid depth the instrument never reached is not "0 pings rejected".

    Same treatment `pg` and `ngood` already get on this path (issue #82).
    """
    proc = _burst_proc(rootdir)
    proc.burst_average_ensembles()

    off = np.isnan(proc.ave.ngood)
    assert off.any(), "no off-profile cell in this file"
    for name in NBAD:
        assert proc.ave[name].dtype == np.float64
        assert np.isnan(proc.ave[name][off]).all()
        assert np.isfinite(proc.ave[name][~off]).all()


def test_burst_counts_ignore_the_interpolated_bin(rootdir):
    """An interpolated bin is not measured data, so its counts stand.

    `pg` and `ngood` are deliberately left alone there and these follow.
    """
    plain = _burst_proc(rootdir)
    plain.burst_average_ensembles()
    interp = _burst_proc(rootdir)
    interp.burst_average_ensembles(interpolate_bin=3)

    for name in NBAD:
        assert np.array_equal(plain.ave[name], interp.ave[name], equal_nan=True), name


def test_every_invalid_grid_cell_on_the_profile_has_a_reason(rootdir, monkeypatch):
    """Gridding loses no attribution.

    `interp1` reads exactly the two bracketing bins, so a grid cell is
    invalid for a ping precisely when at least one bracketing bin was, and
    its reason set is the union of those bins' reasons.
    """
    _, calls = _run(rootdir, monkeypatch, "average_ensembles", method="_edit")

    n_invalid = n_attributed = 0
    for _, _, ens in calls:
        on_profile = ~np.isnan(ens.amp_grid)
        invalid = (~ens.valid_grid) & on_profile
        fired = np.zeros_like(invalid)
        for name in madcp.QC_CRITERIA:
            fired |= ens.reason_grid[name]
        n_invalid += int(invalid.sum())
        n_attributed += int((invalid & fired).sum())

    assert n_invalid > 0
    assert n_attributed == n_invalid


def test_average_ensembles_counts_are_the_gridded_reasons(rootdir, monkeypatch):
    """Each published count is the sum of the gridded reason it names.

    Run with `maskbins` set so `nbad_maskbins` is nonzero and takes part in
    the comparison below; without it a column swap involving that variable
    would pass unnoticed.
    """
    calls = _record_edit(monkeypatch, "_edit")
    proc = _proc(rootdir)
    proc.editparams.maskbins = proc.generate_binmask([4, 5])
    proc.average_ensembles()

    assert calls, "no ensemble reached _edit"
    for _, _, ens in calls:
        for name in madcp.QC_CRITERIA:
            grid = ens.reason_grid[name]
            assert grid.dtype == bool
            assert grid.shape == ens.enu_grid.shape[:2]

    for name in madcp.QC_CRITERIA:
        expected = np.array(
            [np.sum(ens.reason_grid[name], axis=0) for _, _, ens in calls]
        )
        assert np.array_equal(proc.ave[f"nbad_{name}"], expected)
    assert proc.ave["nbad_maskbins"].any(), "maskbins never fired"


def test_burst_counts_match_their_named_reason_in_bin_space(rootdir, monkeypatch):
    """The interpolated count for each criterion is built from its own sum.

    The published values are these bin-space sums interpolated onto the
    depth grid and floored independently, and `floor(a) + floor(b)` can fall
    short of `floor(a + b)` (issue #30), so comparing the published grid
    values against a hand-rolled interpolation would mean reimplementing
    that arithmetic here. Instead this captures the array
    `burst_average_ensembles` hands to `interp1` before any interpolation or
    flooring runs, and checks each of its four reason columns against the
    same sum computed independently by name. That is what would catch the
    counts landing under the wrong `nbad_*` name, for instance a slip
    between the order `counts_inst` is assembled in and the order
    `QC_CRITERIA` is unpacked in afterwards.
    """
    captured = []
    original_interp1 = madcp.interp1
    ncols = 2 + len(madcp.QC_CRITERIA)

    def recorder(*args, **kwargs):
        y_old = args[1]
        if getattr(y_old, "ndim", 0) == 2 and y_old.shape[-1] == ncols:
            captured.append(np.asarray(y_old))
        return original_interp1(*args, **kwargs)

    monkeypatch.setattr(madcp, "interp1", recorder)
    _, calls = _run(rootdir, monkeypatch, "burst_average_ensembles", method="_edit")

    assert len(captured) == len(calls)
    for (_, _, ens), counts_inst in zip(calls, captured):
        for j, name in enumerate(madcp.QC_CRITERIA):
            expected = np.sum(np.asarray(ens[f"reason_{name}"]), axis=0)
            assert np.array_equal(counts_inst[:, 2 + j], expected)


def test_average_ensembles_counts_are_nan_off_the_profile(rootdir):
    """A grid depth the instrument never reached is not "0 pings rejected".

    The bundled file does not reach this on its default grid, so the grid is
    extended past the profile to construct it. This is the `average_ensembles`
    half of the rule `tests/test_qc_flags.py` already checks on the burst path.
    """
    proc = _proc(rootdir, dgridparams={"dtop": 20, "dbot": 400, "d_interval": 4})
    proc.average_ensembles()

    off = np.isnan(proc.ds.amp.values)
    assert off.any(), "the extended grid still lies inside the profile"
    for name in NBAD:
        assert proc.ds[name].dtype == np.float64
        assert np.isnan(proc.ds[name].values[off]).all()
        assert np.isfinite(proc.ds[name].values[~off]).all()
