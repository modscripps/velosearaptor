"""Pinned output of the four processing configurations.

Every other test in this suite asserts on the one quantity its issue was
about: an axis name, a `pg` dtype, an output attribute, a masked cell. A
restructuring of the processing loops can pass all of them while changing
velocities. This module exists to catch that. It compares the complete output
Dataset of four configurations, variable by variable, against a reference
recorded from the last reviewed state of the code.

The reference is `data/pinned_output.json`, a manifest holding, for
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
processing. `velosearaptor_version` changes at every release. Everything else
is compared.

What the pin holds. It was first recorded from the v0.3.0 tag (625a72b) and
has been regenerated since, so it carries the last reviewed output and not
the output of any one release. Every regeneration is a deliberate output
change, is reviewed in the diff of the manifest, and gets a HISTORY.md entry
saying what moved; `git log tests/data/pinned_output.json` lists them. As a
tripwire against a change nobody intended this is as strong as pinning a
release. What it gives up is the question a user with a release installed
would ask, whether the numbers of that version still hold, and answering it
means checking the tag out and regenerating there.
"""

import argparse
import hashlib
import json
import pathlib

import numpy as np
import pytest

from velosearaptor.madcp import ProcessADCP

REFERENCE = "data/pinned_output.json"

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

# Magnitudes below this are folded onto zero before anything is compared.
#
# `e` and `w` are sums and differences of beam pairs and cancel exactly for
# some cells. Whether the cancellation lands on exactly zero or on a residue
# depends on the order the host's numpy chose to add the terms in, so it is a
# property of the machine. Measured on the four configurations pinned here,
# `e` on the per-ping path carries residues down to 4.3e-19, which is 2^-61,
# and the smallest surviving magnitude in any pinned variable is 3.0e-7. The
# threshold sits in the empty range between them, five orders below the
# second number and seven above the first, so it cannot hide a real change.
# A float32 ulp at a typical velocity is 1e-8, four orders above it. The
# snap folds 250 of the 458 near-zero cells of per-ping `e`, whose smallest
# surviving value is then 1.03e-3, the instrument's 1 mm/s velocity quantum.
#
# This folds values onto zero and quantizes nothing else, so a one-ulp change
# anywhere else in the array still moves the checksum.
#
# Without this, `e` and `w` disagreed between CI runs of identical code on
# identical inputs while their minimum, maximum, mean and invalid count
# matched to the last digit. Nothing else in the output has this problem.
NEAR_ZERO = 1e-12


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


def _snap_near_zero(values):
    """Cancellation residues, and the sign of zero, folded onto zero.

    See `NEAR_ZERO`. NaN survives the comparison and is left alone.
    """
    snapped = values.copy()
    snapped[np.abs(snapped) < NEAR_ZERO] = 0
    return snapped


def _checksum(values):
    """sha256 over the filled values and, separately, over the invalid mask.

    Float values are snapped near zero first. NaN is filled, because its
    payload bits are not fixed by the standard; where a cell is invalid is
    hashed separately, so the filling cannot hide a change. Byte order is
    forced little-endian so the manifest does not depend on the host.
    """
    a = np.ascontiguousarray(values)
    if a.dtype.kind == "M":
        invalid = np.isnat(a)
        payload = a.view("int64").copy()
    elif a.dtype.kind in "fc":
        invalid = np.isnan(a)
        payload = _snap_near_zero(a)
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
            # Snapped, like the checksum, so the recorded diagnostics are as
            # host-independent as the quantity they describe.
            good = _snap_near_zero(values)[valid]
            record["min"] = float(good.min())
            record["max"] = float(good.max())
            record["mean"] = float(np.mean(good, dtype=np.float64))
            # `min_abs_nonzero` is what showed that the residues sit twelve
            # orders of magnitude below the smallest real value, which is
            # what makes `NEAR_ZERO` safe. Keep it recorded.
            magnitude = np.abs(good)
            record["n_zero"] = int((magnitude == 0).sum())
            nonzero = magnitude[magnitude > 0]
            record["min_abs_nonzero"] = float(nonzero.min()) if nonzero.size else None
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
    fields = [
        key
        for key in ("min", "max", "mean", "first", "last", "n_zero", "min_abs_nonzero")
        if key in record
    ]
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
            lines.append(f"{name}: present now, not in the pin")
            continue
        if name not in current:
            lines.append(f"{name}: in the pin, absent now")
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
            f"{name}: {len(lines)} difference(s) from the pinned output",
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


def test_the_checksum_snaps_near_zero():
    """Deliberate, and the reason is at `NEAR_ZERO`. Do not undo it."""
    zero = np.array([1.0, 0.0, -2.0], dtype=np.float64)
    for indistinguishable in (-0.0, 4.336808689942018e-19, -1e-13):
        other = np.array([1.0, indistinguishable, -2.0], dtype=np.float64)
        assert zero.tobytes() != other.tobytes()
        assert _checksum(zero) == _checksum(other)

    # The threshold is five orders of magnitude below the smallest real value
    # in any pinned variable, so a change that matters still moves the
    # checksum.
    assert _checksum(zero) != _checksum(np.array([1.0, 1e-10, -2.0]))
    assert _checksum(zero) != _checksum(np.array([1.0, 0.0, 2.0]))

    # So does moving where the invalid cells are.
    one_nan = np.array([1.0, np.nan, -2.0])
    other_nan = np.array([np.nan, 0.0, -2.0])
    assert _checksum(one_nan) != _checksum(zero)
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
