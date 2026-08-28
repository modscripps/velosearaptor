"""Tests for binmapping and masked arrays (issue #78).

`_binmap_one_beam` flags cells that fall outside the mapped beam range with
NaN. `ens.vel`, `ens.amp` and `ens.cor` are masked arrays, so writing NaN
into them as plain data leaves the mask unset. The NaN then travels into
`xyze` and poisons every statistic computed over it, most importantly the
error-velocity standard deviation that sets the adaptive threshold in
`_edit`. `amp` and `cor` are integer typed on top of that, so their NaN does
not even survive the assignment: it is cast to whatever the platform makes
of `uint8(nan)`.
"""

import numpy as np
import pytest
from pycurrents.adcp.transform import Transform
from pycurrents.system import Bunch

from velosearaptor.madcp import ProcessADCP

META_DATA = {
    "mooring": "Test",
    "project": "Test",
    "lon": 0,
    "lat": 0,
}


@pytest.fixture
def adcpfile(rootdir):
    return rootdir / "data/binmap_16670013.000"


@pytest.fixture
def proc(adcpfile):
    return ProcessADCP(adcpfile, META_DATA, magdec=0.0)


@pytest.fixture
def ensemble(proc):
    """A real ensemble from the bundled uplooker."""
    return proc.m.read(start=0, stop=200)


def _invalid(a):
    """Cells that carry no data: masked, or NaN written as data."""
    return ~np.isfinite(np.ma.filled(a, np.nan))


def test_binmapped_velocity_is_masked_not_nan(proc, ensemble):
    """Out-of-range cells must come back masked, and no NaN may be left in
    the unmasked data."""
    proc._binmap_all_beams(ensemble)
    vel = ensemble.vel

    assert np.ma.getmaskarray(vel).sum() > 0
    assert np.isnan(vel.compressed()).sum() == 0


def test_binmapped_amplitude_and_correlation_are_masked(proc, ensemble):
    """`amp` and `cor` are `uint8`, so NaN written into them is destroyed by
    the cast and silently becomes a valid-looking number (0 on x86). Every
    cell `_binmap_one_beam` returned as NaN must end up masked."""
    # `_binmap_one_beam` does not modify the ensemble, so the invalid cells
    # can be collected before binmapping is applied.
    bad = {
        beam: [_invalid(x) for x in proc._binmap_one_beam(ensemble, beam)]
        for beam in [1, 2, 3, 4]
    }

    proc._binmap_all_beams(ensemble)

    for beam, (bad_vel, bad_amp, bad_cor) in bad.items():
        assert bad_vel.sum() > 0
        for field, invalid in [
            (ensemble.vel, bad_vel),
            (ensemble.amp, bad_amp),
            (ensemble.cor, bad_cor),
        ]:
            mask = np.ma.getmaskarray(field)[..., beam - 1]
            assert mask[invalid].all()
            assert not mask.all()


def test_error_velocity_carries_no_unmasked_nan(proc, ensemble):
    """The mask must survive the beam-to-xyz transform, so that the error
    velocity `_edit` computes its threshold from is a statistic over real
    data. Without it `e.std()` is `masked` and the adaptive threshold in
    `_edit` collapses to `editparams["max_e"]`."""
    proc._binmap_all_beams(ensemble)
    proc._calculate_xyze(ensemble, ibad=proc.ibad)

    assert np.ma.getmaskarray(ensemble.xyze).sum() > 0
    assert np.isnan(ensemble.xyze.compressed()).sum() == 0

    e = ensemble.xyze[:, :, 3]
    assert np.isfinite(e.std())


def _synthetic_ensemble(nping=100, nbin=20, angle=20.0, outlier_e=0.15):
    """An ensemble whose error velocity is small and well behaved except in
    the last two pings.

    The beam velocities are constant with depth, so binmapping leaves the
    in-range cells untouched and NaNs only the cells that map outside the
    beam range. A 10 degree roll puts three of the twenty bins out of range.

    Returns the ensemble and the error velocity written into it, one value
    per ping.
    """
    trans = Transform(angle=angle, geometry="convex")
    # Error velocity produced by a unit velocity on beam 1 alone.
    probe = np.zeros((1, 1, 4))
    probe[..., 0] = 1.0
    scale = trans.beam_to_xyz(probe)[0, 0, 3]

    e_target = np.empty(nping)
    e_target[:-2] = np.linspace(-0.03, 0.03, nping - 2)
    e_target[-2:] = outlier_e

    vel = np.zeros((nping, nbin, 4))
    vel[..., 0] = (e_target / scale)[:, None]

    ens = Bunch(
        vel=np.ma.masked_array(vel, mask=False),
        # Integer typed, as they come off the instrument.
        amp=np.ma.masked_array(np.full((nping, nbin, 4), 100, dtype=np.uint8)),
        cor=np.ma.masked_array(np.full((nping, nbin, 4), 128, dtype=np.uint8)),
        pitch=np.zeros(nping),
        roll=np.full(nping, 10.0),
        dep=np.arange(1, nbin + 1) * 10.0,
        sysconfig=Bunch(angle=angle, convex=True),
    )
    return ens, e_target


def test_sigma_criterion_binds_on_binmapped_data(proc):
    """The test that would have caught the bug.

    The synthetic ensemble is built so the sigma criterion is the binding
    one: std(e) is 0.027, so `max_e_deviation * std(e)` is 0.054, well below
    `max_e` = 0.2. The two outlier pings at 0.15 must be edited out.

    `min_correlation` is set to 0 to disable the correlation test. That test
    is the only reason the NaN cells do not reach the error velocity on x86:
    `uint8(nan)` happens to come out as 0, below any correlation threshold,
    so the cells are rejected by undefined behavior standing in for QC.
    """
    proc.parse_editparams({"min_correlation": 0, "max_e": 0.2, "max_e_deviation": 2})
    ep = proc.editparams
    ens, e_target = _synthetic_ensemble()

    # The ensemble is the intended one: the sigma threshold sits below the
    # absolute one. Checked against the values put in, so the check itself
    # does not depend on the masking under test.
    expected_max_e = ep.max_e_deviation * e_target.std()
    assert expected_max_e < ep.max_e

    proc._binmap_all_beams(ens)
    proc._calculate_xyze(ens, ibad=None)

    in_range = ~_invalid(ens.xyze[:, :, 3])
    assert in_range.sum() > 0
    assert (~in_range).sum() > 0

    proc._edit(ens)

    # The adaptive threshold actually binds, at the level set by the data
    # that is really there.
    assert ens.max_e_applied < ep.max_e
    assert np.isclose(ens.max_e_applied, expected_max_e, rtol=1e-3)

    masked = np.ma.getmaskarray(ens.xyze)[:, :, 3]
    # Every in-range cell of the two outlier pings is edited out ...
    assert masked[-2:][in_range[-2:]].all()
    # ... and no in-range cell of the quiet pings is.
    assert not masked[:-2][in_range[:-2]].any()


def test_process_pings_still_runs_with_binmapping(proc):
    """End-to-end guard: the masking must not break the per-ping path."""
    proc.process_pings(binmap=True)
    ds = proc.ds

    assert np.isfinite(ds.u.values).sum() > 0
    assert np.isfinite(ds.max_e_applied.values).all()
