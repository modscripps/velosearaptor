"""Package for dealing with RDI Teledyne ADCP data in various ways. The code
interfaces the UH package `pycurrents` and its `Multiread` for efficient
reading of raw ADCP data. The code is based on the [pycurrents](https://currents.soest.hawaii.edu/hgstage/pycurrents) package, installation instructions can be found [here](https://currents.soest.hawaii.edu/ocn_data_analysis/installation.html).

`gadcp.io` provides convenience functions for reading raw data either into
an `xarray.Dataset` or into the output structure format provided by the UH
software package.

`gadcp.adcp` is a collection of functions that are useful to quickly
analyze raw ADCP data.

`gadcp.madcp` contains functions for processing moored ADCP data. Many thanks to Eric Fiering for sharing his code on moored ADCP data processing [mcm_avp.py](https://currents.soest.hawaii.edu/hgstage/pycurrents/file/tip/pycurrents/adcp/mcm_avg.py) that much of this is based on.

---
# Installation

---
# Changes

This release brings a major refactoring of the `gadcp.madcp` module with lots of breaking changes. The module now allows for improved ensemble averages for burst sampling schemes with better control of the editing parameters.

## v0.2.0 (2022 March)

### New Features
- Ensemble-average before depth gridding. This happens automatically when using `gadcp.madcp.ProcessADCP.burst_average_ensembles`.
- Improved gridding for burst sampling schemes in `gadcp.madcp.ProcessADCP.burst_average_ensembles`.
- Apply pg criterion prior to depth gridding in `gadcp.madcp.ProcessADCP.burst_average_ensembles`.
- Interpolate over missing bin prior to depth gridding in `gadcp.madcp.ProcessADCP.burst_average_ensembles`.
- Write log messages to a file and (if desired) to the screen.

### Breaking Changes
- Changed the `gadcp.madcp` architecture and moved from using the function `gadcp.madcp.proc` to the class `gadcp.madcp.ProcessADCP`.

### Bug Fixes
- Correct time stamp calculation for burst averages.
- For burst sampling schemes average pressure before gridding to depth.

### Documentation
- Added lots of documentation to `gadcp.madcp.ProcessADCP`.
- Use [pdoc](https://pdoc.dev/docs/pdoc.html) to generate the package documentation.
- Automatically build the documentation using [GitHub Actions](https://github.com/gunnarvoet/gadcp/actions).
- Automatically deploy the documentation to [gunnarvoet.github.io/gadcp](https://gunnarvoet.github.io/gadcp/gadcp.html) with GitHub Actions.

### Internal Changes
- Improved the pip/conda requirements to automatically install the [pycurrents](https://currents.soest.hawaii.edu/hgstage/pycurrents) package.


---
# To Do
- Pressure time series handling in `gadcp.madcp.ProcessADCP`. See comments in `gadcp.madcp`.
- Move :class:`gadcp.adcp.LADCPYoYoSplit` to its own ladcp submodule.

"""
__all__ = ["io", "adcp", "madcp"]

__author__ = "Gunnar Voet"
__email__ = "gvoet@ucsd.edu"
__version__ = "0.1"

from . import io, adcp, madcp
