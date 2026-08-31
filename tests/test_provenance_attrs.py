"""Tests for what the output file says about how it was made (issue #102).

Three separate claims are checked here. That the file does not advertise
editing it did not do (`pg_limit`). That it does not point at variables it
does not carry (`ancillary_variables`). And that the processing parameters
which change the numbers are recoverable from the file itself, rather than
only from a log that does not travel with it.

Issue #82's `pg` comment rides along, because `pg` means a different thing on
each of the three paths and the file has never said which one.
"""

import pathlib

import numpy as np
import pytest
import xarray as xr

import velosearaptor
from velosearaptor import io
from velosearaptor.madcp import ProcessADCP

META_DATA = {
    "mooring": "Test",
    "project": "Test",
    "lon": 0,
    "lat": 0,
}

CONTINUOUS = "binmap_16670013.000"
BURST = "24606000.000"

# Attributes every file carries, whichever method produced it.
UNIVERSAL = [
    "velosearaptor_version",
    "processing_method",
    "orientation",
    "magdec",
    "max_e",
    "max_e_deviation",
    "min_correlation",
    "maskbins",
    "pg_limit",
    "depth_offset",
    "pressure_source",
    "pressure_scale_factor",
    "time_drift_rate",
    "clock_correction",
    "ibad",
    "t0",
    "t1",
]


def _proc(rootdir, filename, **kwargs):
    return ProcessADCP(rootdir / "data" / filename, META_DATA, magdec=0.0, **kwargs)


@pytest.fixture(scope="module")
def datasets():
    """One dataset per processing path, built once.

    Module scoped, so it cannot use the function-scoped `rootdir` fixture.
    """
    rootdir = pathlib.Path(__file__).parent.resolve()
    out = {}

    p = ProcessADCP(rootdir / "data" / CONTINUOUS, META_DATA, magdec=0.0)
    p.process_pings()
    out["process_pings"] = p.ds

    p = ProcessADCP(rootdir / "data" / CONTINUOUS, META_DATA, magdec=0.0)
    p.average_ensembles()
    out["average_ensembles"] = p.ds

    p = ProcessADCP(
        rootdir / "data" / BURST,
        META_DATA,
        magdec=0.0,
        tgridparams={"burst_average": True},
    )
    p.burst_average_ensembles()
    out["burst_average_ensembles"] = p.ds

    return out


# --- 11.1 pg_limit ------------------------------------------------------


def test_pg_limit_is_numeric_only_where_it_is_applied(datasets):
    """`pg_limit` screens percent good in `burst_average_ensembles` alone.

    The class docstring already says so. Writing the number into every file
    tells a reader the product was screened at that threshold when it was
    not, and the file carries the evidence: see the test below.
    """
    assert datasets["burst_average_ensembles"].attrs["pg_limit"] == 50
    for method in ["average_ensembles", "process_pings"]:
        assert datasets[method].attrs["pg_limit"] == "not applied"


def test_unscreened_cells_survive_on_the_averaging_path(datasets):
    """The measurement behind the previous test.

    Velocities below the recorded `pg_limit` are present in the file, which
    is correct behavior and exactly why the attribute must not claim
    otherwise.
    """
    ds = datasets["average_ensembles"]
    pg = ds.pg.values
    unscreened = (pg > 0) & (pg < 50) & np.isfinite(ds.u.values)
    assert unscreened.sum() > 0


def test_pg_limit_none_is_still_reported_as_none(rootdir):
    """A user who switches the screening off gets `none`, not `not applied`.

    The two are different statements and the burst path can make either.
    """
    proc = _proc(
        rootdir,
        BURST,
        tgridparams={"burst_average": True},
        editparams={"pg_limit": None},
    )
    proc.burst_average_ensembles()
    assert proc.ds.attrs["pg_limit"] == "none"


# --- 11.2 dangling ancillary_variables ----------------------------------


@pytest.mark.parametrize(
    "method", ["process_pings", "average_ensembles", "burst_average_ensembles"]
)
def test_no_variable_names_an_absent_ancillary_variable(datasets, method):
    """Asserted generically, so a new variable cannot reintroduce this.

    `io.cf_conventions` names `npings` as an ancillary variable on twelve
    entries, and `process_pings` does not write `npings` at all. A CF checker
    fails on that, and a reader looking for the sample count behind `pg` or
    the error estimates finds nothing.
    """
    ds = datasets[method]
    dangling = {}
    for v in ds.variables:
        names = ds[v].attrs.get("ancillary_variables", "").split()
        missing = [n for n in names if n not in ds.variables]
        if missing:
            dangling[v] = missing
    assert dangling == {}


def test_ancillary_reference_is_kept_where_the_variable_exists(datasets):
    """Stripping is targeted: the averaging paths do write `npings`."""
    ds = datasets["average_ensembles"]
    assert "npings" in ds.variables
    assert ds.u.attrs["ancillary_variables"] == "npings"


def test_cf_conventions_still_declares_the_reference(datasets):
    """The static dict is unchanged; the stripping happens per file.

    Recorded so that a later reader does not conclude the CF entries were
    edited and reintroduce the dangling reference by restoring them.
    """
    assert io.cf_conventions()["u"]["ancillary_variables"] == "npings"


# --- 11.3 provenance ----------------------------------------------------


@pytest.mark.parametrize(
    "method", ["process_pings", "average_ensembles", "burst_average_ensembles"]
)
def test_universal_provenance_attributes_present(datasets, method):
    for att in UNIVERSAL:
        assert att in datasets[method].attrs, f"{method} is missing {att}"


@pytest.mark.parametrize(
    "method", ["process_pings", "average_ensembles", "burst_average_ensembles"]
)
def test_processing_method_names_the_method(datasets, method):
    assert datasets[method].attrs["processing_method"] == method


def test_version_is_recorded(datasets):
    """Without it none of the other attributes is interpretable over time."""
    assert datasets["average_ensembles"].attrs["velosearaptor_version"] == (
        velosearaptor.__version__
    )


def test_depth_grid_recorded_only_where_there_is_one(datasets):
    """`process_pings` does not regrid, so it has no depth grid to report."""
    for att in ["dtop", "dbot", "d_interval", "averaging_interval_hours"]:
        assert att not in datasets["process_pings"].attrs
        for method in ["average_ensembles", "burst_average_ensembles"]:
            assert att in datasets[method].attrs

    ds = datasets["average_ensembles"]
    assert ds.attrs["dtop"] == 10
    assert ds.attrs["dbot"] == 165.0
    assert ds.attrs["d_interval"] == 4.0


def test_dt_hours_recorded_only_where_it_governs(datasets):
    """A burst gets its averaging interval from the ping pattern.

    `make_start_ddays` reads `dt_hours` only when `burst_average` is False,
    so writing it into a burst-averaged file names a parameter that had no
    effect, which is the same fault as the old `pg_limit`. The interval that
    did govern is recorded on both paths instead.
    """
    assert "dt_hours" not in datasets["burst_average_ensembles"].attrs
    ds = datasets["average_ensembles"]
    assert ds.attrs["dt_hours"] == 0.5
    assert ds.attrs["averaging_interval_hours"] == 0.5

    # The burst interval is derived, and it is not the unused dt_hours.
    burst = datasets["burst_average_ensembles"].attrs["averaging_interval_hours"]
    assert burst != 0.5
    assert 0 < burst < 1


def test_binmap_recorded_only_on_the_path_that_offers_it(datasets):
    """Binmapping is an argument of `process_pings` and nothing else."""
    assert datasets["process_pings"].attrs["binmap"] == "False"
    for method in ["average_ensembles", "burst_average_ensembles"]:
        assert "binmap" not in datasets[method].attrs


def test_binmap_true_is_recorded(rootdir):
    proc = _proc(rootdir, CONTINUOUS)
    proc.process_pings(binmap=True)
    assert proc.ds.attrs["binmap"] == "True"


def test_interpolate_bin_recorded_only_on_the_burst_path(datasets):
    assert datasets["burst_average_ensembles"].attrs["interpolate_bin"] == "none"
    for method in ["average_ensembles", "process_pings"]:
        assert "interpolate_bin" not in datasets[method].attrs


def test_interpolate_bin_value_is_recorded(rootdir):
    """An interpolated bin is not measured data, so the file must say so."""
    proc = _proc(rootdir, BURST, tgridparams={"burst_average": True})
    proc.burst_average_ensembles(interpolate_bin=5)
    assert proc.ds.attrs["interpolate_bin"] == 5


@pytest.mark.parametrize(
    ("method", "expected"),
    [
        ("process_pings", "raw"),
        ("average_ensembles", "lowpass"),
        ("burst_average_ensembles", "raw"),
    ],
)
def test_pressure_source_recorded(datasets, method, expected):
    """Which pressure set the depths, since the four sources differ.

    `average_ensembles` low-passes by default, `burst_average_ensembles`
    averages raw pressure over the burst, and `process_pings` uses raw
    pressure per ping. The choice moves `xducer_depth` and the depth grid,
    so it is not recoverable from the numbers alone.
    """
    assert datasets[method].attrs["pressure_source"] == expected


def test_use_raw_pressure_is_reflected(rootdir):
    proc = _proc(rootdir, CONTINUOUS, use_raw_pressure=True)
    proc.average_ensembles()
    assert proc.ds.attrs["pressure_source"] == "raw"


def test_external_pressure_is_reflected(rootdir):
    """An external record overrides every other source on every path."""
    proc = _proc(rootdir, CONTINUOUS)
    time = io.yday0_to_datetime64(proc.tsdat.yearbase, proc.tsdat.dday)
    external = xr.DataArray(
        proc.tsdat.pressure.copy(), coords={"time": time}, dims=["time"]
    )
    proc = _proc(rootdir, CONTINUOUS, pressure=external)
    proc.average_ensembles()
    assert proc.ds.attrs["pressure_source"] == "external"


def test_clock_correction_recorded(rootdir, datasets):
    """`driftparams` is not written to the log file at all today."""
    ds = datasets["average_ensembles"]
    assert ds.attrs["clock_correction"] == "none"
    assert ds.attrs["time_drift_rate"] == 1

    proc = _proc(
        rootdir,
        CONTINUOUS,
        driftparams={
            "end_pc": (2018, 9, 10, 11, 36, 0),
            "end_adcp": (2018, 9, 10, 11, 35, 0),
        },
    )
    proc.average_ensembles()
    assert proc.ds.attrs["clock_correction"] == "applied"
    assert proc.ds.attrs["time_drift_rate"] != 1


def test_requested_time_window_recorded_in_calendar_time(datasets):
    """`t0`/`t1` are stored as year days internally, which few readers know."""
    ds = datasets["average_ensembles"]
    for att in ["t0", "t1"]:
        parsed = np.datetime64(ds.attrs[att])
        assert np.datetime64("2000-01-01") < parsed < np.datetime64("2100-01-01")
    assert np.datetime64(ds.attrs["t0"]) < np.datetime64(ds.attrs["t1"])
    # The window brackets the data it produced.
    assert np.datetime64(ds.attrs["t0"]) <= ds.time.values[0]
    assert np.datetime64(ds.attrs["t1"]) >= ds.time.values[-1]


def test_ibad_recorded(rootdir, datasets):
    assert datasets["average_ensembles"].attrs["ibad"] == "none"

    proc = _proc(rootdir, CONTINUOUS, ibad=2)
    proc.average_ensembles()
    assert proc.ds.attrs["ibad"] == 2


def test_pressure_scale_factor_recorded(rootdir, datasets):
    assert datasets["average_ensembles"].attrs["pressure_scale_factor"] == 1

    proc = _proc(rootdir, CONTINUOUS, pressure_scale_factor=1.02)
    proc.average_ensembles()
    assert proc.ds.attrs["pressure_scale_factor"] == 1.02


@pytest.mark.parametrize(
    "method", ["process_pings", "average_ensembles", "burst_average_ensembles"]
)
def test_attributes_round_trip_through_netcdf(datasets, method, tmp_path):
    """netCDF attributes cannot hold None, booleans or empty arrays.

    Every value written has to survive a write and a read unchanged, which
    is the constraint behind the `"none"` strings and the string-valued
    `binmap`.
    """
    ds = datasets[method]
    path = tmp_path / f"{method}.nc"
    ds.to_netcdf(path)

    with xr.open_dataset(path) as reopened:
        for att, value in ds.attrs.items():
            back = reopened.attrs[att]
            if isinstance(value, np.ndarray):
                assert np.array_equal(back, value), att
            else:
                assert back == value, att


# --- 4a: what pg means on this path -------------------------------------


@pytest.mark.parametrize(
    "method", ["process_pings", "average_ensembles", "burst_average_ensembles"]
)
def test_pg_carries_a_per_path_comment(datasets, method):
    """`pg` is a different quantity on each path and the file must say which.

    `io.cf_conventions` builds one static dict and `_add_names_and_units`
    assigns it wholesale, so the per-path text has to be injected after that.
    """
    comment = datasets[method].pg.attrs["comment"]
    assert method in comment
    # It is not the instrument's own percent good, which velosearaptor reads
    # and never uses. Same name, different quantity.
    assert "RDI" in comment


def test_pg_comments_differ_between_paths(datasets):
    comments = {m: datasets[m].pg.attrs["comment"] for m in datasets}
    assert len(set(comments.values())) == 3


def test_pg_comment_describes_the_binary_single_ping_case(datasets):
    """On `process_pings` percent good is 0 or 100 and nothing between."""
    ds = datasets["process_pings"]
    values = np.unique(ds.pg.values[np.isfinite(ds.pg.values)])
    assert set(values) <= {0.0, 100.0}
    assert "0 or 100" in ds.pg.attrs["comment"]


def test_pg_comment_names_interpolated_bins_on_the_burst_path(datasets):
    """An interpolated bin stays at `pg = 0` so it stays visible."""
    comment = datasets["burst_average_ensembles"].pg.attrs["comment"]
    assert "interpolate_bin" in comment


# The two averaging paths count percent good on different grids, deliberately.
# A comment that says only what each one counts leaves that reading as an
# inconsistency, which is the half of issue #82 that documentation closes.
OTHER_AVERAGING_PATH = {
    "average_ensembles": "`burst_average_ensembles`",
    "burst_average_ensembles": "`average_ensembles`",
}


@pytest.mark.parametrize("method", sorted(OTHER_AVERAGING_PATH))
def test_pg_comment_names_the_criterion_behind_the_ordering(datasets, method):
    """Why this path counts where it does, not only what it counts.

    The criterion is whether the transducer moves appreciably within the
    averaging window. A burst shares one depth vector, so its pings share a
    bin grid and the count is exact in bin space; over the longer intervals
    of `average_ensembles` the transducer moves several meters against a 4 m
    bin, so no common bin grid exists and gridding has to come first.
    """
    comment = datasets[method].pg.attrs["comment"]
    assert "transducer" in comment
    # Each one names the other, so the difference reads as a choice.
    assert OTHER_AVERAGING_PATH[method] in comment


def test_pg_comment_says_the_ordering_does_not_arise_on_single_pings(datasets):
    """`process_pings` averages nothing, so it has no ordering to justify."""
    comment = datasets["process_pings"].pg.attrs["comment"]
    assert "averaging window" in comment
