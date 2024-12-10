"""End-to-end AutoMech tests
"""

import contextlib
import os
from pathlib import Path

import automech
import pytest
from automech import test_utils
from automech.test_utils import InvalidSignatureError

ROOT_DIR = Path(__file__).parent.parent
TEST_UTILS = test_utils.TestUtils(ROOT_DIR)

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


@pytest.mark.parametrize("test_dir", TEST_UTILS.test_dirs)
def test_workflow(test_dir: str):
    """Test the entire workflow"""
    print(f"Running in {test_dir}...")
    os.chdir(test_dir)

    with contextlib.redirect_stdout(test_utils.Logger("out.log")):
        automech.run()


if __name__ == "__main__":
    # test_workflow("quick")
    test_signature()
