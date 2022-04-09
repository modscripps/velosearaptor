gadcp
=====

Library of python modules for dealing with RDI ADCP data. Heavily based on pycurrents which take a couple extra steps to install, therefore these libraries were moved from `gvpy` to their own module. Follow [pycurrents installation instructions](https://currents.soest.hawaii.edu/ocn_data_analysis/installation.html).

Import as
```
import gadcp
```

To install in regular mode
```
python setup.py install
```

To install in developer mode
```
python setup.py develop
```

or using pip:
```
pip install -e .
```

## Installing magdec

Note that this method only works with the editable pip install of `gadcp` (e.g. `pip install -e .`) since it compiles `magdec` locally and does not do a system install. See note at bottom for system install.

**Requirements**:
-  A C compiler, e.g. in OSX with [homebrew](https://brew.sh/) use `brew install gcc`
- [`conda` package manager](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
- The gadcp environment installed via conda, e.g. `conda env create -f environment`

With the requirements satisfied, run the shell script:
```bash
./install_magdec.sh
```

To remove, delete the `geomag/` directory.

For a system install you need to `cd geomag` and do something like `sudo make install`.