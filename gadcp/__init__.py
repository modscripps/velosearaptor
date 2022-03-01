"""Package for dealing with RDI Teledyne ADCP data in various ways. The code
interfaces the UH package `pycurrents` and its `Multiread` for efficient
reading of raw ADCP data. The code is based on the [pycurrents](https://currents.soest.hawaii.edu/hgstage/pycurrents) package, installation instructions can be found [here](https://currents.soest.hawaii.edu/ocn_data_analysis/installation.html).

:mod:`gadcp.io` provides convenience functions for reading raw data either into
an `xarray.Dataset` or into the output structure format provided by the UH
software package.

:mod:`gadcp.adcp` is a collection of functions that are useful to quickly
analyze raw ADCP data.

:mod:`gadcp.madcp` contains functions for processing moored ADCP data. Many thanks to Eric Fiering for sharing his code on moored ADCP data processing [mcm_avp.py](https://currents.soest.hawaii.edu/hgstage/pycurrents/file/tip/pycurrents/adcp/mcm_avg.py) that much of this is based on.

---

TODO: Move :fun:`gadcp.adcp.LADCPYoYoSplit` to its own ladcp submodule.

"""
__all__ = ["io", "adcp", "madcp"]

__author__ = "Gunnar Voet"
__email__ = "gvoet@ucsd.edu"
__version__ = "0.1"

from . import io, adcp, madcp
