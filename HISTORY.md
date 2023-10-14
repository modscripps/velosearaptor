## History

### v0.?? (unreleased)

#### New Features
- Install `magdec` via shell script.
- Add an example notebook.
- Read serial number from binary file and compare with SN in meta data ([PR13]( https://github.com/modscripps/velosearaptor/pull/13)). By [Jesse Cusack](https://github.com/jessecusack/).
- Add bin mapping for cases with large pitch and roll. This adds `velosearaptor.madcp.ProcessADCP.process_pings` ([PR17]( https://github.com/modscripps/velosearaptor/pull/17)). By [Jesse Cusack](https://github.com/jessecusack/).
- Add position to moored ADCP meta data ([PR21]( https://github.com/modscripps/velosearaptor/pull/21)).
- Allow for external input of pressure time series in `velosearaptor.madcp.ProcessADCP` ([PR12]( https://github.com/modscripps/velosearaptor/pull/12)).
- Improve default depth grid to also work well with mooring knockdowns ([PR44]( https://github.com/modscripps/velosearaptor/pull/44)).
- Optionally read processing parameters from .yml-file ([PR26]( https://github.com/modscripps/velosearaptor/pull/26)).
- Add CF-compliant meta data to output dataset ([PR26]( https://github.com/modscripps/velosearaptor/pull/26)).

#### Breaking Changes
- Transfer repository from [gunnarvoet](https://github.com/gunnarvoet/) to [modscripps](https://github.com/modscripps/) and rename from gadcp to velosearaptor. Legacy code still exists at [https://github.com/gunnarvoet/gadcp](https://github.com/gunnarvoet/gadcp).
- Change processed dataset coordinate `z` to `depth` ([PR48]( https://github.com/modscripps/velosearaptor/pull/48)).
- Change `vel_std` variables in output dataset to `vel_error` by dividing the standard deviation of each average by the square root of the number of pings ([PR26]( https://github.com/modscripps/velosearaptor/pull/26)).
- Low-pass filter (inherently noisy) pressure before ensemble-averaging continuous ping data ([PR61]( https://github.com/modscripps/velosearaptor/pull/61)). By [Gunnar Voet](https://github.com/gunnarvoet/).

#### Bug Fixes
- Fix conda/pip environment.
- Read correct instrument orientation when a majority of the time series has been recorded outside the water ([PR44]( https://github.com/modscripps/velosearaptor/pull/44)).
- Read paths provided via pathlib.PosixPath objects ([PR55]( https://github.com/modscripps/velosearaptor/pull/55)). By [Gunnar Voet](https://github.com/gunnarvoet/).

#### Documentation
- Consolidate readme and history files.
- Add button with link to source code on GitHub ([PR43]( https://github.com/modscripps/velosearaptor/pull/43)).

#### Internal Changes
- Remove `gvpy` dependency ([PR27]( https://github.com/modscripps/velosearaptor/pull/27)).



### v0.2.0 (2022 March)
This release brings a major refactoring of the `velosearaptor.madcp` module with lots of breaking changes. The module now allows for improved ensemble averages for burst sampling schemes with better control of the editing parameters.

#### New Features
- Ensemble-average before depth gridding. This happens automatically when using `velosearaptor.madcp.ProcessADCP.burst_average_ensembles`.
- Improved gridding for burst sampling schemes in `velosearaptor.madcp.ProcessADCP.burst_average_ensembles`.
- Apply pg criterion prior to depth gridding in `velosearaptor.madcp.ProcessADCP.burst_average_ensembles`.
- Interpolate over missing bin prior to depth gridding in `velosearaptor.madcp.ProcessADCP.burst_average_ensembles`.
- Write log messages to a file and (if desired) to the screen.

#### Breaking Changes
- Changed the `velosearaptor.madcp` architecture and moved from using the function `velosearaptor.madcp.proc` to the class `velosearaptor.madcp.ProcessADCP`.

#### Bug Fixes
- Correct time stamp calculation for burst averages.
- For burst sampling schemes average pressure before gridding to depth.

#### Documentation
- Added lots of documentation to `velosearaptor.madcp.ProcessADCP`.
- Use [pdoc](https://pdoc.dev/docs/pdoc.html) to generate the package documentation.
- Automatically build the documentation using [GitHub Actions](https://github.com/modscripps/velosearaptor/actions).
- Automatically deploy the documentation to [modscripps.github.io/velosearaptor](https://modscripps.github.io/velosearaptor/velosearaptor.html) with GitHub Actions.

#### Internal Changes
- Improved the pip/conda requirements to automatically install the [pycurrents](https://currents.soest.hawaii.edu/hgstage/pycurrents) package.


### v0.0.1 (2020-04-28)
* Moved all ADCP-related functions from `gvpy` to this module.
