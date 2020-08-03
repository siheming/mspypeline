import pytest

from .mock_data import MockData


def pytest_addoption(parser):
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )
    parser.addoption(
        "--runrelease", action="store_true", default=False, help="run release tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")
    config.addinivalue_line("markers", "release: mark test as only relevant for a release")


def pytest_collection_modifyitems(config, items):
    # if config.getoption("--runslow"):
    #     # --runslow given in cli: do not skip slow tests
    #     return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords and not config.getoption("--runslow"):
            item.add_marker(skip_slow)

    skip_release = pytest.mark.skip(reason="need --release option to run")
    for item in items:
        if "release" in item.keywords and not config.getoption("--runrelease"):
            item.add_marker(skip_release)