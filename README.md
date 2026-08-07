velosearaptor
=============

Library of python modules for reading and processing raw RDI Teledyne ADCP data. The code interfaces the UH package [pycurrents](https://currents.soest.hawaii.edu/git/uh-currents-group/shipboard-adcp/pycurrents) and its `Multiread` for efficient reading of raw ADCP data.

`velosearaptor.io` provides convenience functions for reading raw data either into an `xarray.Dataset` or into the output structure format provided by the UH software package.

`velosearaptor.adcp` is a collection of functions that are useful to quickly analyze raw ADCP data.

`velosearaptor.madcp` contains functions for processing moored ADCP data. Many thanks to Eric Fiering for sharing his code on moored ADCP data processing [mcm_avg.py](https://currents.soest.hawaii.edu/git/uh-currents-group/shipboard-adcp/pycurrents/-/blob/main/pycurrents/adcp/mcm_avg.py) that much of this is based on.

## Installation

The project is managed with [uv](https://docs.astral.sh/uv/):

```bash
git clone https://github.com/modscripps/velosearaptor.git
cd velosearaptor
uv sync
```

This creates a virtual environment in `.venv/` with `velosearaptor` installed in
editable mode. Run things inside it with `uv run`, e.g. `uv run pytest`, or see
`make help` for the available development tasks.

### pycurrents

`velosearaptor` depends on [pycurrents](https://github.com/gunnarvoet/pycurrents),
which is not on PyPI. `uv sync` installs it automatically from that GitHub
snapshot of the [UH repository](https://currents.soest.hawaii.edu/git/uh-currents-group/shipboard-adcp/pycurrents).
Since pycurrents builds a number of Cython/C extensions, **a C compiler is
required** (on macOS run `xcode-select --install`).

To develop against a local pycurrents checkout, overlay it into the environment
rather than editing `pyproject.toml`:

```bash
uv pip install -e ../pycurrents --no-deps
```

### Installing magdec

`magdec` computes magnetic declination and is optional; it is only needed when
`lon`/`lat` are provided and no declination is passed in explicitly.

**Requirements**: a C compiler and `make`.

```bash
./install_magdec.sh
```

This clones and builds [geomag](https://currents.soest.hawaii.edu/git/Oceanography_Tools/geomag)
into `geomag/` at the repository root, where `velosearaptor` finds it at
runtime. Because this is a local build rather than a system install, it only
works from a source checkout. To remove, delete the `geomag/` directory.

For a system-wide install, `cd geomag` and do something like `sudo make
install`; `magdec` is then picked up from your `PATH` instead.
