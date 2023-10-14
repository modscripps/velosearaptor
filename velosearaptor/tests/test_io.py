"""Test tools."""
import numpy as np
import pytest

import velosearaptor as vr


class TestTimeConversion:

    def test_single_value(self):
        test_dt64 = np.datetime64("2023-06-01 12:30:22")
        year, test_yd0 = vr.io.datetime64_to_yday0(test_dt64)
        dt64_return = vr.io.yday0_to_datetime64(year, test_yd0)
        assert dt64_return == test_dt64

    def test_array(self):
        test_dt64 = np.arange('2023-06-01 12:00:00', '2023-06-01 12:30:00', dtype='datetime64[s]')
        year, test_yd0 = vr.io.datetime64_to_yday0(test_dt64)
        dt64_return = vr.io.yday0_to_datetime64(year, test_yd0)
        np.testing.assert_equal(dt64_return, test_dt64)

