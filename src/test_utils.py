"""Tests workflow utilities."""

import contextlib
import functools
import os
import re
import shutil
import subprocess
import sys
import tarfile
from collections.abc import Sequence
from pathlib import Path

import yaml

ROOT_PATH = Path(__file__).parent.parent
ARCHIVE_COMMIT_MESSAGE = "Updates tests/archive.tgz"


class Directory:
    """Directories."""

    tests: Path = ROOT_PATH / "tests"
    examples: Path = ROOT_PATH / "examples"


class File:
    """Files."""

    commit: Path = Directory.tests / ".commit"
    tests_yaml: Path = Directory.tests / "tests.yaml"
    archive: Path = Directory.tests / "archive.tgz"


class Test:
    """Test."""

    @classmethod
    @functools.lru_cache
    def names(cls) -> list[str]:
        """Test names."""
        return yaml.safe_load(File.tests_yaml.read_text())

    @classmethod
    def paths(cls) -> list[Path]:
        """Test paths."""
        return [Directory.tests / name for name in cls.names()]

    @classmethod
    def source_paths(cls) -> list[Path]:
        """Test paths."""
        return [Directory.examples / name for name in cls.names()]


# Main testing workflow functions
def setup_tests() -> list[Path]:
    """Set up tests for testing workflow.

    :return: List of test paths
    """
    # 1. Assert that there are no uncommitted Python changes
    changes = uncommitted_python_changes()
    assert not changes, f"You have uncommitted changes:\n{changes}"

    # 2. Write the test commit hash to file
    File.commit.write_text(current_commit_line())

    # 3. Copy input directories over from examples
    test_paths = Test.paths()
    test_source_paths = Test.source_paths()
    for test_path, test_source_path in zip(test_paths, test_source_paths, strict=True):
        if test_path.exists():
            print(f"Removing {test_path}")
            shutil.rmtree(test_path)

        print(f"Creating {test_path} from {test_source_path}")
        shutil.copytree(test_source_path / "inp", test_path / "inp", dirs_exist_ok=True)

    return test_paths


def archive_tests() -> None:
    """Archive tests from testing workflow."""
    exclude = ("subtasks", "run")

    def _filter(obj: tarfile.TarInfo) -> tarfile.TarInfo | None:
        """Filter function for excluding unneeded directories."""
        name = obj.name
        if any(f"{e}/" in name or name.endswith(e) for e in exclude):
            return None
        return obj

    os.chdir(Directory.tests)
    print(os.getcwd())
    print(f"Creating {File.archive}...")
    File.archive.unlink(missing_ok=True)
    with tarfile.open(File.archive, "w:gz") as tar:
        if File.commit.exists:
            tar.add(File.commit, arcname=File.commit.name)
        for test in Test.names():
            tar.add(test, arcname=test, filter=_filter)


def extract_archived_tests() -> None:
    """Extract archived tests."""
    with contextlib.chdir(Directory.tests):
        if File.archive.exists():
            print(f"Unpacking {File.archive}...")
            with tarfile.open(File.archive, "r") as tar:
                tar.extractall()


def commit_test_archive() -> None:
    """Commit the test archive to the MechDriver git repo."""
    subprocess.run(["git", "reset"], cwd=ROOT_PATH)
    subprocess.run(["git", "add", str(File.archive)], cwd=ROOT_PATH)
    subprocess.run(["git", "commit", "-m", ARCHIVE_COMMIT_MESSAGE], cwd=ROOT_PATH)


# Helper functions
def pixi_activation_hook() -> str:
    """Get pixi activation hook."""
    return subprocess.check_output(["pixi", "shell-hook"], text=True, cwd=ROOT_PATH)


def uncommitted_python_changes() -> str:
    """Get uncommitted python changes."""
    return subprocess.check_output(
        ["git", "status", "-s", "*.py"], text=True, cwd=ROOT_PATH / "src"
    )


def current_commit_line(
    skip: Sequence[str] = (
        re.escape(ARCHIVE_COMMIT_MESSAGE),
        r"Merge pull request \S* from \S*",
        r"Merge \S* into \S*",
    ),
) -> str:
    """Get the first commit line from the log (oneline) output.

    :param log: Log output
    :param skip: Regexes to skip
    :return: The first commit line
    """
    log = subprocess.check_output(["git", "log", "--oneline"], text=True)
    if not log:
        return log

    lines = log.splitlines()
    return next((line for line in lines if not re.search("|".join(skip), line)), "")


class Logger:
    """Logger for capturing standard out."""

    def __init__(self, file_name: str = "out.log"):  # noqa: D107
        self.stdout = sys.stdout
        self.file = open(file_name, "w")

    def write(self, message):  # noqa: D102
        self.stdout.write(message)
        self.file.write(message)
        self.stdout.flush()
        self.file.flush()

    def flush(self):  # noqa: D102
        pass
