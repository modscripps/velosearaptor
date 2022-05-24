"""A few examples"""
import numpy as np
import pytest
from pycurrents.system import Bunch

from gadcp.madcp import ProcessADCP


def test_binmap(rootdir):
    adcpfile = rootdir / "data/binmap_16670013.000"

    meta_data = {
        "mooring": "Test",
        "project": "Test",
        "lon": 0,
        "lat": 0,
    }
    ADCP = ProcessADCP(adcpfile, meta_data)

    # Make some fake data
    ntime, ndepth, nbeam = 1, 5, 4  # data points in each axis
    dep = np.array([5, 10, 50, 100, 200])  # bin distances in m
    angle_deg = np.array([20])  # Beam angle in degrees
    angel_rad = np.deg2rad(angle_deg)
    vel = np.ones((ntime, ndepth, nbeam)) * np.array([1, 2, 3, 4, 5]).reshape(
        (ntime, ndepth, 1)
    )

    # Negative pitch means that the instrument is tilted in the direction of beam 3
    pitch_deg = np.array([-5])
    # Positive roll is tilt in the direction of beam 1
    roll_deg = np.array([10])
    pitch_rad = np.deg2rad(pitch_deg)
    roll_rad = np.deg2rad(roll_deg)

    # Make a fake ensemble
    ens = Bunch(
        {
            "vel": vel.copy(),
            "amp": vel.copy(),
            "cor": vel.copy(),
            "pitch": pitch_deg.copy(),
            "roll": roll_deg.copy(),
            "dep": dep.copy(),
            "sysconfig": Bunch({"angle": angle_deg.copy()}),
        }
    )

    beam = 1
    dep3 = dep * np.cos(pitch_rad) * np.cos(angel_rad + roll_rad) / np.cos(angel_rad)
    vel3i = np.interp(dep, dep3, vel[0, :, beam - 1])
    outside_range = (dep > dep3[-1]) | (dep < dep3[0])
    inside_range = ~outside_range

    # Edit ensemble in place
    veli, ampi, cori = ADCP._binmap_one_beam(ens, beam)

    assert np.all(np.isclose(veli[0, inside_range], vel3i[inside_range]))
    assert np.all(np.isclose(ampi[0, inside_range], vel3i[inside_range]))
    assert np.all(np.isclose(cori[0, inside_range], vel3i[inside_range]))
    assert np.all(np.isnan(veli[0, outside_range]))
    assert np.all(np.isnan(ampi[0, outside_range]))
    assert np.all(np.isnan(cori[0, outside_range]))

    beam = 2
    dep3 = dep * np.cos(pitch_rad) * np.cos(angel_rad - roll_rad) / np.cos(angel_rad)
    vel3i = np.interp(dep, dep3, vel[0, :, beam - 1])
    outside_range = (dep > dep3[-1]) | (dep < dep3[0])
    inside_range = ~outside_range

    # Edit ensemble in place
    veli, ampi, cori = ADCP._binmap_one_beam(ens, beam)

    assert np.all(np.isclose(veli[0, inside_range], vel3i[inside_range]))
    assert np.all(np.isclose(ampi[0, inside_range], vel3i[inside_range]))
    assert np.all(np.isclose(cori[0, inside_range], vel3i[inside_range]))
    assert np.all(np.isnan(veli[0, outside_range]))
    assert np.all(np.isnan(ampi[0, outside_range]))
    assert np.all(np.isnan(cori[0, outside_range]))

    beam = 3
    dep3 = dep * np.cos(roll_rad) * np.cos(angel_rad - pitch_rad) / np.cos(angel_rad)
    vel3i = np.interp(dep, dep3, vel[0, :, beam - 1])
    outside_range = (dep > dep3[-1]) | (dep < dep3[0])
    inside_range = ~outside_range

    # Edit ensemble in place
    veli, ampi, cori = ADCP._binmap_one_beam(ens, beam)

    assert np.all(np.isclose(veli[0, inside_range], vel3i[inside_range]))
    assert np.all(np.isclose(ampi[0, inside_range], vel3i[inside_range]))
    assert np.all(np.isclose(cori[0, inside_range], vel3i[inside_range]))
    assert np.all(np.isnan(veli[0, outside_range]))
    assert np.all(np.isnan(ampi[0, outside_range]))
    assert np.all(np.isnan(cori[0, outside_range]))

    beam = 4
    dep3 = dep * np.cos(roll_rad) * np.cos(angel_rad + pitch_rad) / np.cos(angel_rad)
    vel3i = np.interp(dep, dep3, vel[0, :, beam - 1])
    outside_range = (dep > dep3[-1]) | (dep < dep3[0])
    inside_range = ~outside_range

    # Edit ensemble in place
    veli, ampi, cori = ADCP._binmap_one_beam(ens, beam)

    assert np.all(np.isclose(veli[0, inside_range], vel3i[inside_range]))
    assert np.all(np.isclose(ampi[0, inside_range], vel3i[inside_range]))
    assert np.all(np.isclose(cori[0, inside_range], vel3i[inside_range]))
    assert np.all(np.isnan(veli[0, outside_range]))
    assert np.all(np.isnan(ampi[0, outside_range]))
    assert np.all(np.isnan(cori[0, outside_range]))

    # Run the all beam solution in case of errors
    ADCP._binmap_all_beams(ens)
    # Check that things were changed but not that the values are correct (since we do that above)
    assert np.all(~np.isclose(ens.vel, vel))

    # Check the full processing MOSTLY changes things too
    close_frac = 0.01  # Max fraction allowed to be close, 1%
    ADCP.process_pings(binmap=True)
    ds_binmap = ADCP.ds.copy()
    ADCP.process_pings(binmap=False)
    ds_no_binmap = ADCP.ds.copy()
    uclose = np.isclose(ds_binmap.u, ds_no_binmap.u)
    assert np.sum(uclose) / uclose.size < close_frac
    vclose = np.isclose(ds_binmap.v, ds_no_binmap.v)
    assert np.sum(vclose) / vclose.size < close_frac
    wclose = np.isclose(ds_binmap.w, ds_no_binmap.w)
    assert np.sum(wclose) / wclose.size < close_frac
    eclose = np.isclose(ds_binmap.e, ds_no_binmap.e)
    assert np.sum(eclose) / eclose.size < close_frac
    ampclose = np.isclose(ds_binmap.amp, ds_no_binmap.amp)
    assert np.sum(ampclose) / ampclose.size < close_frac
