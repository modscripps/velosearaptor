"""
.. include:: ../../README.md

---

.. include:: ../../HISTORY.md

"""

# NOTE: `yoyo` is intentionally omitted. It still references the `gv` (gvpy)
# and `plt` names that were dropped when the gvpy dependency was removed, so it
# raises NameError as soon as any of its methods are called.
__all__ = ["adcp", "io", "madcp", "tools"]

__author__ = "velosearaptor Developers"
__email__ = "gvoet@ucsd.edu"
__version__ = "0.3.0"

from . import adcp, io, madcp, tools
