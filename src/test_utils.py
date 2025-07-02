"""Tests workflow utilities."""

import contextlib
import functools
import os
import re
import shutil
import subprocess
import sys
import tarfile
import textwrap
from pathlib import Path

import pydantic
import yaml

ROOT_PATH = Path(__file__).parent.parent
ARCHIVE_COMMIT_MESSAGE = "Update tests/archive.tgz"
SKIP_COMMITS = (
    re.escape(ARCHIVE_COMMIT_MESSAGE),
    r"Merge pull request \S* from \S*",
    r"Merge \S* into \S*",
    r"Update README.md",
)


class Directory:
    """Directories."""

    tests: Path = ROOT_PATH / "tests"
    examples: Path = ROOT_PATH / "examples"


class File:
    """Files."""

    commit: Path = Directory.tests / ".commit"
    signature: Path = Directory.tests / "signature.yaml"
    tests: Path = Directory.tests / "tests.yaml"
    archive: Path = Directory.tests / "archive.tgz"


class Test:
    """Test."""

    @classmethod
    @functools.lru_cache
    def names(cls) -> list[str]:
        """Test names."""
        return yaml.safe_load(File.tests.read_text())

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
    """Set up tests to prepare local testing workflow.

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


def wrap_up_tests(from_archive: bool = False, allow_override: bool = False) -> None:
    """Archive and sign tests to wrap up local testing workflow."""
    if from_archive:
        extract_archived_tests()

    sign_tests(allow_override=allow_override)
    archive_tests()
    commit_test_archive()


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
        if File.commit.exists():
            tar.add(File.commit, arcname=File.commit.name)
        if File.signature.exists():
            tar.add(File.signature, arcname=File.signature.name)
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


class Signature(pydantic.BaseModel):
    """Sign off on local tests."""

    signed_commit: str
    tested_commit: str
    untested_commits: list[str]
    username: str


def sign_tests(allow_override: bool = False) -> None:
    """Sign off on local tests."""
    # Assert that there are no uncommitted Python changes
    changes = uncommitted_python_changes()
    assert not changes, f"You have uncommitted changes:\n{changes}"

    # Get commits since test commit
    test_commit = File.commit.read_text()
    curr_commit = current_commit_line()
    new_commits = commit_lines_since(test_commit)

    # Interactively approve
    if new_commits:
        if not allow_override:
            print(f"Not signing tests because new commits are present:\n{new_commits}")
            return
        else:
            print("WARNING: New commits since tested version!!")
            print(textwrap.indent("\n".join(new_commits), "    "))
            answer = input(
                "Do you solemnly swear that these changes will not break tests? (yes/no): "
            )
            print()
            if answer != "yes":
                print("Thank you for your honesty.")
                print("Please re-run the tests using `pixi run test local`.")
                sys.exit()

    # Sign
    sign = Signature(
        signed_commit=curr_commit,
        tested_commit=test_commit,
        untested_commits=new_commits,
        username=github_username(),
    )
    print(f"Writing signed repo information to {File.signature}")
    write_signature(sign)


# Signature file I/O
def write_signature(sign: Signature):
    """Write signature file."""
    File.signature.write_text(yaml.safe_dump(sign.model_dump()))


def read_signature() -> Signature:
    """Read signature file."""
    return Signature.model_validate(yaml.safe_load(File.signature.read_text()))


# Helper functions
def pixi_activation_hook() -> str:
    """Get pixi activation hook."""
    return subprocess.check_output(["pixi", "shell-hook"], text=True, cwd=ROOT_PATH)


def uncommitted_python_changes() -> str:
    """Get uncommitted python changes."""
    return subprocess.check_output(
        ["git", "status", "-s", "*.py"], text=True, cwd=ROOT_PATH / "src"
    )


def github_username() -> str:
    """Return the current user's GitHub username as a string.

    Requires the username to be configured:

        git config --global user.name

    :return: The username
    """
    return subprocess.check_output(
        ["git", "config", "--global", "user.name"], text=True
    ).strip()


def current_commit_line() -> str:
    """Get the first commit line from the log (oneline) output.

    :param log: Log output
    :return: The first commit line
    """
    log = subprocess.check_output(["git", "log", "--oneline"], text=True)
    if not log:
        return log

    lines = log.splitlines()
    return next(
        (line for line in lines if not re.search("|".join(SKIP_COMMITS), line)), ""
    )


def commit_lines_since(commit: str) -> list[str]:
    """Get intervening commits between two comments."""
    hash_ = commit_hash(commit)
    log = subprocess.check_output(
        ["git", "log", "--oneline", f"{hash_}..HEAD"], text=True
    )
    return [
        line for line in log.splitlines() if not re.search("|".join(SKIP_COMMITS), line)
    ]


def commit_hash(commit: str) -> str:
    """Get the commit hash from a one-line log summary.

    :param line: A one-line log summary
    :return: The commit hash
    """
    hash_, *_ = commit.split()
    return hash_


def commits_are_equivalent(commit1: str, commit2: str) -> bool:
    """Check if two git commit hashes are equivalent.

    :param commit1: First commit
    :param commit2: Second commit
    :return: `True` if they are, `False` if they aren't
    """
    hash1 = commit_hash(commit1)
    hash2 = commit_hash(commit2)
    nchars = max(4, min(*map(len, [hash1, hash2])))
    return hash1[:nchars] == hash2[:nchars]


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
