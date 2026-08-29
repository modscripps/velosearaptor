"""
.. include:: ../../README.md

---

.. include:: ../../HISTORY.md

"""

__all__ = ["adcp", "io", "madcp", "tools"]

__author__ = "velosearaptor Developers"
__email__ = "gvoet@ucsd.edu"
# A `.devN` suffix marks development between releases. Keep this identical to
# `version` in pyproject.toml. `madcp.ProcessADCP._add_meta_data_to_ds` writes
# this string into every output file as the `velosearaptor_version` attribute,
# so a file produced from an untagged commit has to stay distinguishable from
# one produced from the release. Drop the suffix in both files in the release
# commit, then reopen the next `.dev0` after tagging.
__version__ = "0.4.0.dev0"

from . import adcp, io, madcp, tools
