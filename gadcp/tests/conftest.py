import pathlib

import pytest


@pytest.fixture
def rootdir():
    return pathlib.Path(__file__).parent.resolve()
