gadcp
=====

Library of python modules for reading and processing raw RDI Teledyne ADCP data. The code interfaces the UH package [pycurrents](https://currents.soest.hawaii.edu/hgstage/pycurrents) and its `Multiread` for efficient reading of raw ADCP data.

`gadcp.io` provides convenience functions for reading raw data either into an `xarray.Dataset` or into the output structure format provided by the UH software package.

`gadcp.adcp` is a collection of functions that are useful to quickly analyze raw ADCP data.

`gadcp.madcp` contains functions for processing moored ADCP data. Many thanks to Eric Fiering for sharing his code on moored ADCP data processing [mcm_avp.py](https://currents.soest.hawaii.edu/hgstage/pycurrents/file/tip/pycurrents/adcp/mcm_avg.py) that much of this is based on.

## Installation

### Installing pycurrents

This package depends on the [pycurrents](https://currents.soest.hawaii.edu/hgstage/pycurrents) package. Installation instructions can either be found [here](https://currents.soest.hawaii.edu/ocn_data_analysis/installation.html), or a conda environment including all necessary dependencies can be installed via `conda env create -f environment.yml`. Take a look at `environment.yml` and `requirements.txt` and copy the relevant parts if you would like to incude the package and its dependencies in other environments.


### Installing magdec

Note that this method only works with the editable pip install of `gadcp` (e.g. `pip install -e .`) since it compiles `magdec` locally and does not do a system install. See note at bottom for system install.

**Requirements**:
-  A C compiler, e.g. in OSX with [homebrew](https://brew.sh/) use `brew install gcc`
- [`conda` package manager](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
- The gadcp environment installed via conda, e.g. `conda env create -f environment.yml`

With the requirements satisfied, run the shell script:
```bash
./install_magdec.sh
```

To remove, delete the `geomag/` directory.

For a system install you need to `cd geomag` and do something like `sudo make install`.
