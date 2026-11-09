"""Test tools."""

import numpy as np

import velosearaptor as vr


def test_raw_z_attrs(rootdir):
    """`z` in the raw output is distance from the transducer, not depth.

    The attributes must say so, since the name alone does not.
    """
    ds = vr.io.read_raw_rdi(rootdir / "data/24606000.000")

    assert ds.z.attrs["units"] == "m"
    assert ds.z.attrs["long_name"] == "distance from transducer"
    assert "not water depth" in ds.z.attrs["comment"]
    # No CF standard name exists for this quantity, so claiming one would be
    # wrong. In particular it must not claim to be `depth`.
    assert "standard_name" not in ds.z.attrs

    # Guard the reason the distinction matters: this instrument sits near
    # 1160 dbar while z spans only tens to hundreds of metres, so z cannot be
    # read as water depth.
    assert ds.z.max() < float(ds.pressure.median())


def test_raw_serial_number(rootdir):
    """The raw output carries the instrument serial number from the header.

    It must match the serial number in the file name of the bundled file, as
    an int under the name the processed output uses, on both read paths.
    """
    for aux in [False, True]:
        ds = vr.io.read_raw_rdi(rootdir / "data/24606000.000", auxillary_only=aux)
        assert ds.attrs["sn"] == 24606
        assert type(ds.attrs["sn"]) is int


class TestTimeConversion:
    def test_single_value(self):
        test_dt64 = np.datetime64("2023-06-01 12:30:22")
        year, test_yd0 = vr.io.datetime64_to_yday0(test_dt64)
        dt64_return = vr.io.yday0_to_datetime64(year, test_yd0)
        assert dt64_return == test_dt64

    def test_array(self):
        test_dt64 = np.arange(
            "2023-06-01 12:00:00", "2023-06-01 12:30:00", dtype="datetime64[s]"
        )
        year, test_yd0 = vr.io.datetime64_to_yday0(test_dt64)
        dt64_return = vr.io.yday0_to_datetime64(year, test_yd0)
        np.testing.assert_equal(dt64_return, test_dt64)
