"""End-to-end Workflow tests
"""

import contextlib
import os
from pathlib import Path

import mechdriver
import pytest
import yaml
from mechdriver import test_utils
from mechdriver.test_utils import InvalidSignatureError

ROOT_DIR = Path(__file__).parent.parent
TEST_UTILS = test_utils.TestUtils(ROOT_DIR)
TEST_DIRS = [n for n in yaml.safe_load(TEST_UTILS.config_file.read_text())]


TEST_UTILS.extract_archived_tests()


def test_signature():
    """Unpack archive and check provenance"""
    # 1. Check the hashes of the mechdriver repo (this repo)
    signed_commit = TEST_UTILS.signed_mechdriver_commit()
    current_commit = TEST_UTILS.current_mechdriver_commit()
    if not test_utils.are_equivalent_commits(signed_commit, current_commit):
        raise InvalidSignatureError(f"\n   {signed_commit}\n!~ {current_commit}")

    # 2. Check the hashes of the other remote repos
    signed_prov = dict(TEST_UTILS.signed_provenance())
    current_prov = dict(TEST_UTILS.remote_provenance())
    for repo in ("autochem", "autoio", "autofile", "mechanalyzer"):
        signed_commit = signed_prov[repo]
        current_commit = current_prov[repo]
        if not test_utils.are_equivalent_commits(signed_commit, current_commit):
            raise InvalidSignatureError(f"\n   {signed_commit}\n!~ {current_commit}")


@pytest.mark.parametrize("test_dir", TEST_DIRS)
def test_workflow(test_dir: str):
    """Test the entire workflow."""
    print(f"Running in {test_dir}...")

    with contextlib.chdir(TEST_UTILS.tests_dir / test_dir):
        with contextlib.redirect_stdout(test_utils.Logger("out.log")):
            mechdriver.run()


if __name__ == "__main__":
    # test_workflow("quick")
    test_signature()
