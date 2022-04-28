## History

### v0.?? (unreleased)

#### New Features
- Install `magdec` via shell script.
- Add an example notebook.

<!-- #### Breaking Changes -->

#### Bug Fixes
- Fix conda/pip environment.

#### Documentation
- Consolidate readme and history files.

<!-- #### Internal Changes -->



### v0.2.0 (2022 March)
This release brings a major refactoring of the `gadcp.madcp` module with lots of breaking changes. The module now allows for improved ensemble averages for burst sampling schemes with better control of the editing parameters.

#### New Features
- Ensemble-average before depth gridding. This happens automatically when using `gadcp.madcp.ProcessADCP.burst_average_ensembles`.
- Improved gridding for burst sampling schemes in `gadcp.madcp.ProcessADCP.burst_average_ensembles`.
- Apply pg criterion prior to depth gridding in `gadcp.madcp.ProcessADCP.burst_average_ensembles`.
- Interpolate over missing bin prior to depth gridding in `gadcp.madcp.ProcessADCP.burst_average_ensembles`.
- Write log messages to a file and (if desired) to the screen.

#### Breaking Changes
- Changed the `gadcp.madcp` architecture and moved from using the function `gadcp.madcp.proc` to the class `gadcp.madcp.ProcessADCP`.

#### Bug Fixes
- Correct time stamp calculation for burst averages.
- For burst sampling schemes average pressure before gridding to depth.

#### Documentation
- Added lots of documentation to `gadcp.madcp.ProcessADCP`.
- Use [pdoc](https://pdoc.dev/docs/pdoc.html) to generate the package documentation.
- Automatically build the documentation using [GitHub Actions](https://github.com/gunnarvoet/gadcp/actions).
- Automatically deploy the documentation to [gunnarvoet.github.io/gadcp](https://gunnarvoet.github.io/gadcp/gadcp.html) with GitHub Actions.

#### Internal Changes
- Improved the pip/conda requirements to automatically install the [pycurrents](https://currents.soest.hawaii.edu/hgstage/pycurrents) package.


### v0.0.1 (2020-04-28)
* Moved all ADCP-related functions from `gvpy` to this module.
