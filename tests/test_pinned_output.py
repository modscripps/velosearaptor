"""Pinned output of the four processing configurations, as released in v0.3.0.

Every other test in this suite asserts on the one quantity its issue was
about: an axis name, a `pg` dtype, an output attribute, a masked cell. A
restructuring of the processing loops can pass all of them while changing
velocities. This module exists to catch that. It compares the complete output
Dataset of four configurations, variable by variable, against a reference
recorded from v0.3.0.

The reference is `data/pinned_output_v0.3.0.json`, a manifest holding, for
every coordinate and data variable, its dims, shape, dtype, attributes and a
sha256 over the filled values, plus the dataset attributes and their types.
Regenerate it with

    uv run python tests/test_pinned_output.py --write

and read the resulting diff carefully: a changed line there is a changed
number in somebody's output file.

Why a manifest and not a checked-in netCDF or a bare hash. A netCDF reference
of these four configurations is 4.9 MB, git stores a fresh copy of it on every
intentional output change, and comparing against it pins the netCDF encoding
of whichever xarray version wrote it as well as the numbers. A bare hash is
small but reports only that something moved. The manifest is text, so the diff
names the variable that changed and shows its min, max and mean on both sides,
and it is taken from the in-memory Dataset, so no serialization sits between
the pipeline and the pin.

The checksums, dtypes, dims, shapes and attributes are asserted. The summary
statistics are diagnostics for the failure message and for the diff, and are
not asserted; they are there to tell a float reassociation apart from a bug.

Two dataset attributes are excluded. `proc time` is the wall clock at
processing. `velosearaptor_version` changes at every release and reads
`0.4.0.dev0` on the branch that pins `0.3.0`. Everything else is compared.

v0.3.0 is tagged at 625a72b. It is the reference because it is a release a
user can have installed, which is the form the question "did this change any
output" has to take once the QC-flag work of issue #30 reaches its breaking
steps.
"""

import argparse
import hashlib
import json
import pathlib

import numpy as np
import pytest

from velosearaptor.madcp import ProcessADCP

PINNED_RELEASE = "0.3.0"
REFERENCE = "data/pinned_output_v0.3.0.json"

# Attributes that differ between two runs of identical code, or between two
# releases that compute identical numbers.
EXCLUDED_ATTRS = ("proc time", "velosearaptor_version")

META_DATA = {
    "mooring": "Test",
    "project": "Test",
    "lon": 0,
    "lat": 0,
}

UPLOOKER = "data/binmap_16670013.000"
# The only bundled file that reaches the burst path.
BURST_FILE = "data/24606000.000"


def _proc(rootdir, filename, **kwargs):
    # magdec explicitly, so the pin does not depend on the optional magdec
    # executable being installed.
    return ProcessADCP(rootdir / filename, META_DATA, magdec=0.0, **kwargs)


def _build_average_ensembles(rootdir):
    proc = _proc(rootdir, UPLOOKER)
    proc.average_ensembles()
    return proc.ds


def _build_process_pings_binmap_false(rootdir):
    proc = _proc(rootdir, UPLOOKER)
    proc.process_pings(binmap=False)
    return proc.ds


def _build_process_pings_binmap_true(rootdir):
    proc = _proc(rootdir, UPLOOKER)
    proc.process_pings(binmap=True)
    return proc.ds


def _build_burst_average_ensembles(rootdir):
    proc = _proc(rootdir, BURST_FILE, tgridparams={"burst_average": True})
    proc.burst_average_ensembles()
    return proc.ds


CONFIGURATIONS = {
    "average_ensembles": {
        "description": "binmap_16670013.000 via average_ensembles, defaults",
        "build": _build_average_ensembles,
    },
    "process_pings_binmap_false": {
        "description": "binmap_16670013.000 via process_pings(binmap=False)",
        "build": _build_process_pings_binmap_false,
    },
    "process_pings_binmap_true": {
        "description": "binmap_16670013.000 via process_pings(binmap=True)",
        "build": _build_process_pings_binmap_true,
    },
    "burst_average_ensembles": {
        "description": (
            "24606000.000 via burst_average_ensembles, "
            'tgridparams={"burst_average": True}'
        ),
        "build": _build_burst_average_ensembles,
    },
}


# ---------------------------------------------------------------- manifest


def _checksum(values):
    """sha256 over the filled values and, separately, over the invalid mask.

    Two bit patterns are canonicalized first, because neither is reproducible
    across machines and neither carries any information about the data.

    NaN, whose payload bits are not fixed by the standard. Where a cell is
    invalid is hashed separately, so filling cannot hide a change.

    Negative zero. The sign of a computed zero follows the signs of the
    summands and whether the compiler contracted a multiply-add, so it varies
    between hosts. `w` and `e` are sums and differences of beam pairs of
    quantized velocities and hit exact zero often, and on the bundled
    per-ping configuration they were the only two variables to disagree
    between a local run and CI while their minimum, maximum, mean and invalid
    count all matched to the last digit.

    Little-endian throughout, for the same reason.
    """
    a = np.ascontiguousarray(values)
    if a.dtype.kind == "M":
        invalid = np.isnat(a)
        payload = a.view("int64").copy()
    elif a.dtype.kind in "fc":
        invalid = np.isnan(a)
        payload = a.copy()
        payload[payload == 0] = 0
    elif a.dtype.kind in "iub":
        invalid = np.zeros(a.shape, dtype=bool)
        payload = a
    else:
        return hashlib.sha256(repr(a.tolist()).encode()).hexdigest()
    payload[invalid] = 0
    payload = np.ascontiguousarray(payload.astype(payload.dtype.newbyteorder("<")))
    digest = hashlib.sha256()
    digest.update(payload.tobytes())
    digest.update(np.ascontiguousarray(invalid).tobytes())
    return digest.hexdigest()


def _summary(values):
    """Diagnostics for the failure message and the diff. Never asserted on."""
    kind = values.dtype.kind
    if kind == "M":
        valid = ~np.isnat(values)
        if not valid.any():
            return {"n_invalid": int(values.size)}
        return {
            "n_invalid": int((~valid).sum()),
            "first": str(values[valid].min()),
            "last": str(values[valid].max()),
        }
    if kind in "fc":
        valid = np.isfinite(values)
        record = {"n_invalid": int((~valid).sum())}
        if valid.any():
            good = values[valid]
            record["min"] = float(good.min())
            record["max"] = float(good.max())
            record["mean"] = float(np.mean(good, dtype=np.float64))
        return record
    if kind in "iub" and values.size:
        return {
            "n_invalid": 0,
            "min": int(values.min()),
            "max": int(values.max()),
            "mean": float(np.mean(values, dtype=np.float64)),
        }
    return {"n_invalid": 0}


def _canonical(value):
    """A JSON-representable form of an attribute value."""
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, np.ndarray):
        value = value.tolist()
    if isinstance(value, (list, tuple)):
        return [_canonical(item) for item in value]
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    return repr(value)


def _attrs_record(attrs):
    """Attribute values and their Python types.

    The type is recorded because `np.float64(4.0)` and `4.0` serialize
    differently and a refactor can swap one for the other silently.
    """
    return {
        key: {"type": type(value).__name__, "value": _canonical(value)}
        for key, value in sorted(attrs.items())
        if key not in EXCLUDED_ATTRS
    }


def _variable_record(da, kind):
    values = np.asarray(da.values)
    record = {
        "kind": kind,
        "dims": list(da.dims),
        "shape": list(values.shape),
        "dtype": values.dtype.str,
        "checksum": _checksum(values),
        "attrs": _attrs_record(da.attrs),
    }
    record.update(_summary(values))
    return record


def _manifest(ds):
    variables = {
        name: _variable_record(ds[name], "coord") for name in sorted(ds.coords)
    }
    variables.update(
        {name: _variable_record(ds[name], "data_var") for name in sorted(ds.data_vars)}
    )
    return {
        "dims": {name: int(size) for name, size in sorted(ds.sizes.items())},
        "variables": variables,
        "attrs": _attrs_record(ds.attrs),
    }


# -------------------------------------------------------------- comparison


def _summary_line(record):
    fields = [key for key in ("min", "max", "mean", "first", "last") if key in record]
    parts = [f"{key}={record[key]!r}" for key in fields]
    parts.append(f"n_invalid={record.get('n_invalid')}")
    return ", ".join(parts)


def _compare_dims(reference, current):
    lines = []
    for name in sorted(set(reference) | set(current)):
        if reference.get(name) != current.get(name):
            lines.append(
                f"dimension {name}: {reference.get(name)} -> {current.get(name)}"
            )
    return lines


def _compare_variables(reference, current):
    lines = []
    for name in sorted(set(reference) | set(current)):
        if name not in reference:
            lines.append(f"{name}: present now, not in the v{PINNED_RELEASE} pin")
            continue
        if name not in current:
            lines.append(f"{name}: in the v{PINNED_RELEASE} pin, absent now")
            continue
        ref, got = reference[name], current[name]
        for field in ("kind", "dims", "shape", "dtype"):
            if ref[field] != got[field]:
                lines.append(f"{name}.{field}: {ref[field]!r} -> {got[field]!r}")
        if ref["checksum"] != got["checksum"]:
            lines.append(
                f"{name}: values changed\n"
                f"    pinned:  {_summary_line(ref)}\n"
                f"    current: {_summary_line(got)}"
            )
    return lines


def _compare_attrs(reference, current, label):
    lines = []
    for key in sorted(set(reference) | set(current)):
        if key not in reference:
            lines.append(f"{label}: attribute {key!r} is new")
            continue
        if key not in current:
            lines.append(f"{label}: attribute {key!r} was dropped")
            continue
        ref, got = reference[key], current[key]
        if ref["value"] != got["value"]:
            lines.append(f"{label}.{key}: {ref['value']!r} -> {got['value']!r}")
        elif ref["type"] != got["type"]:
            lines.append(f"{label}.{key}: type {ref['type']} -> {got['type']}")
    return lines


def _compare_variable_attrs(reference, current):
    lines = []
    for name in sorted(set(reference) & set(current)):
        lines.extend(
            _compare_attrs(reference[name]["attrs"], current[name]["attrs"], name)
        )
    return lines


def _message(name, lines):
    return "\n".join(
        [
            (
                f"{name}: {len(lines)} difference(s) from the pinned "
                f"v{PINNED_RELEASE} output"
            ),
            *(f"  {line}" for line in lines),
            "",
            "If the change is intended, regenerate the pin with",
            "    uv run python tests/test_pinned_output.py --write",
            "and describe the numerical change in HISTORY.md.",
        ]
    )


# ------------------------------------------------------------------- tests

_DATASETS = {}


def _dataset(rootdir, name):
    """The output Dataset for one configuration, built once per session."""
    if name not in _DATASETS:
        _DATASETS[name] = CONFIGURATIONS[name]["build"](rootdir)
    return _DATASETS[name]


def _reference(rootdir):
    with open(rootdir / REFERENCE) as fp:
        return json.load(fp)


@pytest.mark.parametrize("name", sorted(CONFIGURATIONS))
def test_pinned_variables(rootdir, name):
    """Dims, dtypes and values of every coordinate and data variable."""
    pinned = _reference(rootdir)["configurations"][name]
    current = _manifest(_dataset(rootdir, name))
    lines = _compare_dims(pinned["dims"], current["dims"])
    lines += _compare_variables(pinned["variables"], current["variables"])
    assert not lines, _message(name, lines)


@pytest.mark.parametrize("name", sorted(CONFIGURATIONS))
def test_pinned_variable_attributes(rootdir, name):
    pinned = _reference(rootdir)["configurations"][name]
    current = _manifest(_dataset(rootdir, name))
    lines = _compare_variable_attrs(pinned["variables"], current["variables"])
    assert not lines, _message(name, lines)


@pytest.mark.parametrize("name", sorted(CONFIGURATIONS))
def test_pinned_dataset_attributes(rootdir, name):
    pinned = _reference(rootdir)["configurations"][name]
    current = _manifest(_dataset(rootdir, name))
    lines = _compare_attrs(pinned["attrs"], current["attrs"], name)
    assert not lines, _message(name, lines)


def test_the_pin_covers_every_variable_in_the_output(rootdir):
    """A variable absent from the manifest would be pinned by nobody."""
    reference = _reference(rootdir)
    assert sorted(reference["configurations"]) == sorted(CONFIGURATIONS)
    for name in CONFIGURATIONS:
        ds = _dataset(rootdir, name)
        pinned = set(reference["configurations"][name]["variables"])
        assert pinned == set(ds.coords) | set(ds.data_vars)


def test_the_pin_excludes_only_the_two_varying_attributes(rootdir):
    """`proc time` and the version string, and nothing else, are left out."""
    reference = _reference(rootdir)
    for name in CONFIGURATIONS:
        ds = _dataset(rootdir, name)
        pinned = set(reference["configurations"][name]["attrs"])
        assert set(ds.attrs) - pinned == set(EXCLUDED_ATTRS)


def test_the_checksum_ignores_the_sign_of_zero():
    """Deliberate, and the reason is in `_checksum`. Do not undo it."""
    positive = np.array([1.0, 0.0, -2.0], dtype=np.float32)
    negative = np.array([1.0, -0.0, -2.0], dtype=np.float32)
    assert positive.tobytes() != negative.tobytes()
    assert _checksum(positive) == _checksum(negative)

    # Everything else about a float array still has to move the checksum,
    # including where the invalid cells are.
    assert _checksum(positive) != _checksum(np.array([1.0, 0.0, 2.0], np.float32))
    one_nan = np.array([1.0, np.nan, -2.0], dtype=np.float32)
    other_nan = np.array([np.nan, 0.0, -2.0], dtype=np.float32)
    assert _checksum(one_nan) != _checksum(positive)
    assert _checksum(one_nan) != _checksum(other_nan)


def test_the_harness_detects_a_perturbation(rootdir):
    """The pin is worth nothing unless it fails on a changed output.

    One ulp on one velocity, a widened dtype, a renamed variable and a moved
    attribute, each of which a restructuring could produce, and each of which
    every other test in this suite would let through.
    """
    name = "process_pings_binmap_false"
    reference = _manifest(_dataset(rootdir, name))

    ds = _dataset(rootdir, name).copy(deep=True)
    values = ds.u.values
    finite = np.flatnonzero(np.isfinite(values))
    assert finite.size > 0
    flat = values.ravel()
    flat[finite[0]] = np.nextafter(flat[finite[0]], np.float32(np.inf))
    lines = _compare_variables(reference["variables"], _manifest(ds)["variables"])
    assert [line for line in lines if line.startswith("u: values changed")], lines

    ds = _dataset(rootdir, name).copy(deep=True)
    ds["w"] = ds.w.astype(np.float64)
    lines = _compare_variables(reference["variables"], _manifest(ds)["variables"])
    assert any(line.startswith("w.dtype:") for line in lines), lines

    ds = _dataset(rootdir, name).copy(deep=True).rename({"amp": "amplitude"})
    lines = _compare_variables(reference["variables"], _manifest(ds)["variables"])
    assert any(line.startswith("amp: ") for line in lines), lines
    assert any(line.startswith("amplitude: ") for line in lines), lines

    ds = _dataset(rootdir, name).copy(deep=True)
    ds.attrs["min_correlation"] = 32
    lines = _compare_attrs(reference["attrs"], _manifest(ds)["attrs"], name)
    assert any("min_correlation" in line for line in lines), lines

    ds = _dataset(rootdir, name).copy(deep=True)
    ds.u.attrs["units"] = "cm s-1"
    lines = _compare_variable_attrs(reference["variables"], _manifest(ds)["variables"])
    assert any("units" in line for line in lines), lines


# -------------------------------------------------------------- regenerate


def _write_reference(rootdir):
    reference = {
        "pinned_release": PINNED_RELEASE,
        "generated_by": "uv run python tests/test_pinned_output.py --write",
        "excluded_attributes": list(EXCLUDED_ATTRS),
        "configurations": {},
    }
    for name, configuration in CONFIGURATIONS.items():
        manifest = _manifest(configuration["build"](rootdir))
        manifest["description"] = configuration["description"]
        reference["configurations"][name] = manifest
    path = rootdir / REFERENCE
    with open(path, "w") as fp:
        json.dump(reference, fp, indent=1, sort_keys=True)
        fp.write("\n")
    print(f"wrote {path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--write",
        action="store_true",
        help="overwrite the reference manifest with the current output",
    )
    if parser.parse_args().write:
        _write_reference(pathlib.Path(__file__).parent.resolve())
    else:
        parser.error("nothing to do without --write")
