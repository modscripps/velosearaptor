import pathlib

import xarray as xr
from pycurrents.adcp.rdiraw import Multiread, extract_raw

import gadcp


def test_success():
    assert True


# We defined rootdir as a fixture in conftest.py
# and can use it here as input now
def test_read_data(rootdir, tmpdir):
    adcpfile = rootdir / "data/03160000.000"
    assert type(adcpfile) == pathlib.PosixPath
    print(adcpfile)
    assert adcpfile.exists()
    # Read configuration.
    m = Multiread(adcpfile.as_posix(), "wh")
    print(m.sysconfig)

    # # make sure we can write and read the data as netcdf
    # p = pathlib.Path(tmpdir) / "testfile.nc"
    # cx.to_netcdf(p)
    # cx2 = xr.open_dataset(p)
    # assert type(cx2) == xr.core.dataset.Dataset
