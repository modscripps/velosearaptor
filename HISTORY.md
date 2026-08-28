## History

### v0.3.0 (2026 August)

#### New Features
- Install `magdec` via shell script.
- Add an example notebook.
- Read serial number from binary file and compare with SN in meta data ([PR13]( https://github.com/modscripps/velosearaptor/pull/13)). By [Jesse Cusack](https://github.com/jessecusack/).
- Add bin mapping for cases with large pitch and roll. This adds `velosearaptor.madcp.ProcessADCP.process_pings` ([PR17]( https://github.com/modscripps/velosearaptor/pull/17)). By [Jesse Cusack](https://github.com/jessecusack/).
- Add position to moored ADCP meta data ([PR21]( https://github.com/modscripps/velosearaptor/pull/21)).
- Allow for external input of pressure time series in `velosearaptor.madcp.ProcessADCP` ([PR12]( https://github.com/modscripps/velosearaptor/pull/12)).
- Improve default depth grid to also work well with mooring knockdowns ([PR44]( https://github.com/modscripps/velosearaptor/pull/44)).
- Optionally read processing parameters from .yml-file ([PR26]( https://github.com/modscripps/velosearaptor/pull/26)).
- Allow for string format when supplying start and end time in .yml parameter file ([PR63]( https://github.com/modscripps/velosearaptor/pull/63)).
- Add CF-compliant meta data to output dataset ([PR26]( https://github.com/modscripps/velosearaptor/pull/26)).

#### Breaking Changes
- Transfer repository from [gunnarvoet](https://github.com/gunnarvoet/) to [modscripps](https://github.com/modscripps/) and rename from gadcp to velosearaptor. Legacy code still exists at [https://github.com/gunnarvoet/gadcp](https://github.com/gunnarvoet/gadcp).
- Change processed dataset coordinate `z` to `depth` ([PR48]( https://github.com/modscripps/velosearaptor/pull/48)).
- Change `vel_std` variables in output dataset to `vel_error` by dividing the standard deviation of each average by the square root of the number of pings ([PR26]( https://github.com/modscripps/velosearaptor/pull/26)).
- Low-pass filter (inherently noisy) pressure before ensemble-averaging continuous ping data ([PR61]( https://github.com/modscripps/velosearaptor/pull/61)). By [Gunnar Voet](https://github.com/gunnarvoet/).
- Rename the vertical coordinate of `velosearaptor.madcp.ProcessADCP.process_pings` output from `depth` to `z`. Any script indexing `ds.depth` on that path breaks, and the values it was reading were never water depth. For an uplooker the axis also runs the opposite way, so previously published per-ping files are inverted relative to water depth as well as offset from it: on `tests/data/binmap_16670013.000` (mean `xducer_depth` 154.5 m) `depth = 5.15` was 149 m of water depth and `depth = 161.15` was above the sea surface. Reprocess, or flip and offset by hand ([PR103](https://github.com/modscripps/velosearaptor/pull/103)).

#### Bug Fixes
- Fix conda/pip environment.
- Read correct instrument orientation when a majority of the time series has been recorded outside the water ([PR44]( https://github.com/modscripps/velosearaptor/pull/44)).
- Read paths provided via pathlib.PosixPath objects ([PR55]( https://github.com/modscripps/velosearaptor/pull/55)). By [Gunnar Voet](https://github.com/gunnarvoet/).
- Fix xarray warning in `groupby` ([PR55]( https://github.com/modscripps/velosearaptor/pull/64)). By [Gunnar Voet](https://github.com/gunnarvoet/).
- Replace scipy.stats.mode with np.unique for dominant period calculation ([#66] ( https://github.com/modscripps/velosearaptor/pull/66)).
- Fix non-monotonic ADCP time vectors by interpolating isolated bad pings, truncating segment overlaps, or raising on ambiguous cases ([PR70](https://github.com/modscripps/velosearaptor/pull/70)).
- Raise a clear error in `ProcessADCP._parse_sysconfig` when no pressure record exceeds 15 dbar. The existing guard could never fire, so this case surfaced as a confusing `IndexError`.
- Stamp the processing date in the log header in UTC rather than local time.
- `velosearaptor.madcp.ProcessADCP.process_pings` now publishes its vertical axis as `z`, the distance from the transducer to the center of each bin, carrying the `z` attributes from `velosearaptor.io.cf_conventions`. The axis was never water depth on that path. The two averaging methods do grid onto water depth and keep publishing `depth`. Bin depth is recoverable per ping as `xducer_depth` plus `z` for a downlooker and `xducer_depth` minus `z` for an uplooker ([PR103](https://github.com/modscripps/velosearaptor/pull/103)).
- Fix four crashes on documented argument paths that no test exercised ([PR104](https://github.com/modscripps/velosearaptor/pull/104)). `velosearaptor.madcp.ProcessADCP.process_pings(start=..., stop=...)` raised `UnboundLocalError`. Both arguments are ping indices, which the docstring now states along with the fact that they mean something different from `average_ensembles(start, stop)`, and they are now range checked, an empty range included. `burst_average=True` on continuous data now raises a clear "this data is not burst sampled" error rather than computing `-2147483648` pings per burst and failing later in an unrelated reduction. `npings` and `ngood` are stored as `int32`, so averaging intervals longer than 32767 pings (a 12 h average at a 1 s ping rate) no longer raise `OverflowError` or wrap silently. This widens the dtype of both variables in the output file. `interpolate_bin` within two bins of either profile end no longer indexes negative bins or runs off the end. The neighbour window is clipped, and the first and last bin, which have no neighbour on one side, raise a clear error. An interpolated bin keeps its own `pg` and `ngood`, so the interpolation stays visible in the output.
- Binmapping no longer writes unmasked NaN into `vel`, `amp` and `cor` ([PR105](https://github.com/modscripps/velosearaptor/pull/105)). Cells outside the mapped beam range are masked instead. `amp` and `cor` are integer typed, so the NaN written into them was destroyed by the cast and became 0, and those zeros were averaged over the four beams into the published amplitude. On the bundled binmapped test file that biased 4.3% of cells low by a median of 22 counts, so anyone who processed with `binmap=True` should reprocess for amplitude. Velocity is unchanged at the default `min_correlation` of 64, because the `uint8(nan)` cast happened to reject exactly the cells the mask now rejects. That the correlation test rejected them at all rested on `uint8(nan)` evaluating to 0, which is undefined behavior, so the rejection is now explicit rather than incidental.
- Convert `t0`/`t1` given as strings against the `yearbase` of the raw file rather than against the year of the requested time itself ([PR106](https://github.com/modscripps/velosearaptor/pull/106)). A requested window in a later calendar year than the first ping was silently off by a year, which is the normal case for a mooring spanning New Year. `velosearaptor.madcp.ProcessADCP.parse_tgridparams` now also raises when the requested window does not overlap the data at all, and when it is empty or reversed, naming the requested window and the file's actual range in calendar time. A window that reaches beyond the data on one side still passes quietly.
- Fix external pressure NaN patching, which raised `TypeError: only integer scalar arrays can be converted to a scalar index` whenever the external pressure record had more than one contiguous NaN run ([PR107](https://github.com/modscripps/velosearaptor/pull/107)). That is the leading-plus-trailing case the code was written for. Leading and trailing NaN runs recorded on deck are now both patched with atmospheric pressure. A NaN gap between good data, an end run recorded in the water, or an external record that is NaN across the whole ADCP time range now raises rather than flowing into the interpolator as NaN.
- Compute the pressure low-pass cutoff on wall-clock time rather than on the modal ping interval ([PR108](https://github.com/modscripps/velosearaptor/pull/108)). `scipy.signal.filtfilt` treats every ping as equally spaced, so on burst-sampled data the documented 30 minute cap was expressed in a fictitious time base and never bound. A 20 dbar knockdown came back as 8 dbar in a synthetic test. `fs` now comes from the mean sampling interval over the record, and a warning is raised when the ping pattern is detectably non-uniform. Continuous records are unaffected. A record whose mean ping interval puts the cutoff at or below the Nyquist period now raises with a clear message, where it previously failed inside `scipy.signal.butter`.

#### Documentation
- Consolidate readme and history files.
- Add button with link to source code on GitHub ([PR43]( https://github.com/modscripps/velosearaptor/pull/43)).
- Describe the `z` coordinate of the raw output in its attributes ([#52](https://github.com/modscripps/velosearaptor/issues/52)). `z` is the distance from the transducer to each bin, not water depth, and keeps its name for that reason. The averaging methods regrid onto water depth and publish it as `depth`, while `velosearaptor.madcp.ProcessADCP.process_pings` keeps the transducer-relative axis and publishes it as `z`.

#### Internal Changes
- Remove `gvpy` dependency ([PR27]( https://github.com/modscripps/velosearaptor/pull/27)).
- Auto-detect sonar type instead of hardcoding Workhorse ([PR65] (https://github.com/modscripps/velosearaptor/pull/65)).
- Speed up bin-averaging ([PR67] ( https://github.com/modscripps/velosearaptor/pull/67 ).
- Move to a modern `pyproject.toml` / [uv](https://docs.astral.sh/uv/) setup with a `src/` layout ([#69](https://github.com/modscripps/velosearaptor/issues/69)). Removes `setup.py`, `setup.cfg`, `requirements.txt`, and `environment.yml`; tests move to a top-level `tests/` directory.
- Install `pycurrents` from its git snapshot instead of the retired Mercurial repository, which fixes the documentation build ([#68](https://github.com/modscripps/velosearaptor/issues/68)).
- Point the `magdec` install at the new geomag git repository and drop the conda requirement ([#71](https://github.com/modscripps/velosearaptor/issues/71)).
- Replace isort/black/flake8 with [ruff](https://docs.astral.sh/ruff/), and clean up the findings so `make check` passes. Linting, format checking, and the test suite now all run in CI on pull requests.
- Remove `velosearaptor.yoyo`. Every method referenced the `gv` (gvpy) and `plt` names that were dropped when the gvpy dependency was removed in [PR27](https://github.com/modscripps/velosearaptor/pull/27), so the module raised `NameError` on any call.


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
