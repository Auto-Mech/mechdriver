"""End-to-end Workflow tests."""

import contextlib

import pytest

import mechdriver
import test_utils as tu
from test_utils import Directory, Logger, Test

TESTS = Test.names()
PROMPT_TESTS = filter(lambda test: "prompt" in test, TESTS)
OTHER_TESTS = filter(lambda test: "prompt" not in test, TESTS)

# Extract tests
tu.extract_archived_tests()


def test_signature():
    """Check signature."""
    sign = tu.read_signature()
    curr_commit = tu.current_commit_line()
    assert tu.commits_are_equivalent(
        sign.signed_commit, curr_commit
    ), f"\n{sign.signed_commit} !~\n{curr_commit}"


@pytest.mark.parametrize("test", OTHER_TESTS)
def test_other_workflow(test: str):
    """Test the entire workflow."""
    print(f"Running in {test}...")

    with contextlib.chdir(Directory.tests / test):
        with contextlib.redirect_stdout(Logger("out.log")):
            mechdriver.run()


@pytest.mark.parametrize("test", PROMPT_TESTS)
def test_prompt_workflow(test: str):
    """Test the entire workflow."""
    print(f"Running in {test}...")

    with contextlib.chdir(Directory.tests / test):
        with contextlib.redirect_stdout(Logger("out.log")):
            mechdriver.run()


if __name__ == "__main__":
    test_other_workflow("quick")
