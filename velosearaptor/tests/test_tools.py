"""Test tools."""
import numpy as np
import pytest

import velosearaptor as vr


def test_timedelta64_to_s():
    test = np.timedelta64(10, "ms")
    s = vr.tools.timedelta64_to_s(test)
    assert s == 0.01


def test_dominant_period_in_s():
    # Generate time vector with 1 minute period
    dt64 = np.arange(
        "2000-01-01 12:00:00", "2000-01-01 12:10:00", dtype="datetime64[m]"
    )
    dom_period = vr.tools.dominant_period_in_s(dt64)
    assert dom_period == 60
