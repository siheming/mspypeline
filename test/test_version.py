import pytest
import json
from urllib.request import urlopen
from mspypeline import __version__ as version

@pytest.mark.release
def test_version_unique():
    url = "https://pypi.org/pypi/mspypeline/json"
    data = json.load(urlopen(url))
    versions = list(data["releases"].keys())
    assert version not in versions, "version already exists"


@pytest.mark.release
def test_release_version():
    # github_version.txt is created in the release workflow on github
    with open("github_version.txt", "r") as f:
        github_version = f.readline()
    assert f"v{version}" == github_version.rstrip(), "project version does not match github release version"
