"""Utility functions for testing."""

import os
import shutil
import subprocess
import sys
import tarfile
import textwrap
import warnings
from collections import defaultdict
from collections.abc import Sequence
from pathlib import Path

import yaml
from pydantic import BaseModel

# constants
ARCHIVE_COMMIT_MESSAGE = "Updates tests/archive.tgz"
EXCLUDE_GREP_REGEX = r"\|".join(
    (
        ARCHIVE_COMMIT_MESSAGE,
        r"Merge pull request \S* from \S*",
        r"Merge \S* into \S*",
    )
)
EXCLUDE_GREP_ARGS = ["--invert-grep", "--grep", EXCLUDE_GREP_REGEX]
EXCLUDE_COMMAND_DCT = {"mechdriver": EXCLUDE_GREP_ARGS}


# exceptions
class InvalidSignatureError(Exception):
    def __init__(self, message, *_):
        # Call the base class constructor with the parameters it needs
        lines = [
            message,
            "  You need to sign off on your local tests using `pixi run test sign`",
            "  (This must be done for an up-to-date version of the code)",
        ]
        super().__init__("\n".join(lines))


class UncommittedChangesError(Exception):
    pass


START_YELLOW = "\033[93m"
END_COLOR = "\033[0m"
warnings.formatwarning = lambda m, c, *_: f"{START_YELLOW}{c.__name__}{END_COLOR}: {m}"


# data classes
class Provenance(BaseModel):
    autochem: str
    autoio: str
    autofile: str
    mechanalyzer: str
    mechdriver: str


class Signature(BaseModel):
    provenance: Provenance
    overrides: dict[str, list[str]]
    username: str


# test helper object (just avoids having to pass root dir around)
class TestUtils:
    def __init__(self, root_dir: str | Path):
        self.root_dir = Path(root_dir)

    # directories
    @property
    def src_dir(self) -> Path:
        return self.root_dir / ".."

    @property
    def mechdriver_dir(self) -> Path:
        return self.root_dir

    @property
    def tests_dir(self) -> Path:
        return self.root_dir / "tests"

    @property
    def examples_dir(self) -> Path:
        return self.root_dir / "examples"

    # data / config files
    @property
    def archive_file(self) -> Path:
        return self.tests_dir / "archive.tgz"

    @property
    def signature_file(self) -> Path:
        return self.tests_dir / "signature.yaml"

    @property
    def provenance_file(self) -> Path:
        return self.tests_dir / "provenance.yaml"

    @property
    def config_file(self) -> Path:
        return self.tests_dir / "config.yaml"

    # test examples
    @property
    def test_names(self) -> list[str]:
        return [n for n in yaml.safe_load(self.config_file.read_text())]

    @property
    def test_dirs(self) -> list[Path]:
        return [self.tests_dir / n for n in self.test_names]

    @property
    def example_dirs(self) -> list[Path]:
        return [self.examples_dir / n for n in self.test_names]

    # git commands
    def repos_output(
        self,
        command: list[str],
        command_dct: dict[str, list[str]] | None = None,
    ) -> dict[str, str]:
        """Get output from a command in each repo

        :param command: A command to run for each repo
        :param command_dct: Additional repo-specific commands, by repo name
        :return: One-line summaries of most recent commits
        """
        command_dct = {} if command_dct is None else command_dct
        output_dct = {}
        for repo in Provenance.model_fields:
            repo_dir = self.src_dir / repo
            full_command = command + command_dct.get(repo, [])
            output_dct[repo] = subprocess.check_output(
                full_command, text=True, cwd=repo_dir
            ).strip()
        return output_dct

    def check_for_uncommited_python_changes(self, throw_error: bool = False) -> None:
        """Assert that there are no uncommited Python changes in any repo.

        :param throw_error: Throw an error if there are uncommitted changes?
        """
        status_dct = self.repos_output(["git", "status", "-s", "*.py"])
        for repo, status in status_dct.items():
            if status:
                message = f"{repo} has uncommitted Python changes:\n{status}\n"
                if throw_error:
                    message += "Please commit changes before running tests!\n"
                    raise UncommittedChangesError(message)
                else:
                    warnings.warn(message)

    def current_provenance(self) -> Provenance:
        """Get information about the current version of each repo

        :return: One-line summaries of most recent commits
        """
        prov_dct = self.repos_output(
            ["git", "log", "--oneline", "-1"], command_dct=EXCLUDE_COMMAND_DCT
        )
        return Provenance(**prov_dct)

    def current_mechdriver_commit(self, hash_only: bool = False) -> str:
        line = subprocess.check_output(
            ["git", "log", "--oneline", "-1", *EXCLUDE_GREP_ARGS],
            text=True,
            cwd=self.mechdriver_dir,
        )
        if hash_only:
            return commit_hash_from_line(line)
        return line

    def signed_provenance(self) -> Provenance:
        signature = self.read_signature()
        return signature.provenance

    def signed_mechdriver_commit(self, hash_only: bool = False) -> str:
        prov = self.signed_provenance()
        line = prov.mechdriver
        if hash_only:
            return commit_hash_from_line(line)
        return line

    def commit_test_archive(self) -> None:
        """Commit the test archive to the MechDriver git repo."""
        subprocess.run(["git", "restore", "--staged", "."], cwd=self.mechdriver_dir)
        subprocess.run(["git", "add", str(self.archive_file)], cwd=self.mechdriver_dir)
        subprocess.run(
            ["git", "commit", "-m", ARCHIVE_COMMIT_MESSAGE], cwd=self.mechdriver_dir
        )

    # read/write files
    def write_provenance(self, prov: Provenance) -> None:
        print(f"\nWriting provenance to {self.provenance_file}")
        self.provenance_file.write_text(yaml.safe_dump(prov.model_dump()))

    def read_provenance(self) -> Provenance:
        print(f"\nReading provenance from {self.provenance_file}")
        return Provenance(**yaml.safe_load(self.provenance_file.read_text()))

    def write_signature(self, sign: Signature) -> None:
        print(f"\nWriting signature to {self.signature_file}")
        self.signature_file.write_text(yaml.safe_dump(sign.model_dump()))

    def read_signature(self) -> Signature:
        print(f"\nReading signature from {self.signature_file}")
        return Signature(**yaml.safe_load(self.signature_file.read_text()))

    # local test workflow functions
    def setup_tests(self) -> None:
        """Sets up clean test directories and saves provenance information."""
        self.check_for_uncommited_python_changes(throw_error=True)
        self.write_provenance(self.current_provenance())
        # Copy input directories over from examples
        for test_dir, example_dir in zip(
            self.test_dirs, self.example_dirs, strict=True
        ):
            if test_dir.exists():
                print(f"Removing {test_dir}")
                shutil.rmtree(test_dir)

            print(f"Creating {test_dir} from {example_dir}")
            shutil.copytree(example_dir / "inp", test_dir / "inp", dirs_exist_ok=True)

    def archive_tests(self) -> None:
        """Archives tests."""
        exclude = ("subtasks", "run")

        def _filter(obj: tarfile.TarInfo) -> tarfile.TarInfo | None:
            """Filter function for excluding unneeded directories."""
            name = obj.name
            if any(f"{e}/" in name or name.endswith(e) for e in exclude):
                return None
            return obj

        os.chdir(self.tests_dir)
        print(f"Creating {self.archive_file}...")
        self.archive_file.unlink(missing_ok=True)
        with tarfile.open(self.archive_file, "w:gz") as tar:
            if self.provenance_file.exists():
                tar.add(self.provenance_file, arcname=self.provenance_file.name)
            if self.signature_file.exists():
                tar.add(self.signature_file, arcname=self.signature_file.name)
            for test in self.test_dirs:
                tar.add(test, arcname=test.name, filter=_filter)

    def extract_archived_tests(self) -> None:
        """Extract archived tests."""
        os.chdir(self.tests_dir)
        if self.archive_file.exists():
            print(f"Unpacking {self.archive_file}...")
            with tarfile.open(self.archive_file, "r") as tar:
                tar.extractall()

    # signing workflow functions
    def provenance_diff(self) -> dict:
        """Determine the differences between the current provenance and the one stored
        in the provenance file."""
        prov = self.read_provenance()
        command_dct = defaultdict(list)
        command_dct.update(EXCLUDE_COMMAND_DCT)
        for module, commit_line in dict(prov).items():
            commit_hash = commit_hash_from_line(commit_line)
            command_dct[module].append(f"{commit_hash}..HEAD")
        diff_dct = self.repos_output(
            ["git", "log", "--oneline"], command_dct=command_dct
        )
        diff_dct = {k: v.splitlines() for k, v in diff_dct.items() if v}
        return diff_dct

    def interactively_approve_provenance_diff_overrides(self) -> dict[str, list[str]]:
        """Interactively override provenance diffs, with user input."""
        diff_dct = self.provenance_diff()
        overrides = {}
        print("This assumes that you first ran local tests using `pixi run test local`")
        print(
            "Checking that the tests were run against the current repository versions..."
        )
        for repo, diffs in diff_dct.items():
            print(
                f"WARNING: {repo} does not match the tested version!! "
                f"It has {len(diffs)} additional commits:"
            )
            print(textwrap.indent("\n".join(diffs), "    "))
            answer = input(
                "Do you solemnly swear that these changes will not affect the tests? (yes/no): "
            )
            print()
            if answer == "yes":
                overrides[repo] = diffs
            else:
                print("Thank you for your honesty.")
                print("Please re-run the tests using `pixi run test local`.")
                return
        return overrides

    def sign_tests(self) -> None:
        """Sign off on local tests."""
        # Compare to tested version and prompt for override as needed
        overrides = self.interactively_approve_provenance_diff_overrides()

        # Sign
        signature = Signature(
            provenance=self.current_provenance(),
            overrides=overrides,
            username=github_username(),
        )
        print(f"\nWriting signed repo information to {self.signature_file}")
        self.write_signature(signature)

    def remote_provenance(self) -> Provenance:
        """Get provenance information from remote repos."""
        commit_dct = {}
        for repo in Provenance.model_fields:
            git_url = f"https://github.com/Auto-Mech/{repo}.git"
            commit_dct[repo] = subprocess.check_output(
                ["git", "ls-remote", git_url, "HEAD"], text=True
            ).strip()
        return Provenance(**commit_dct)


# logger
class Logger(object):
    def __init__(self, file_name: str = "out.log"):
        self.stdout = sys.stdout
        self.file = open(file_name, "a")

    def write(self, message):
        self.stdout.write(message)
        self.file.write(message)
        self.stdout.flush()
        self.file.flush()

    def flush(self):
        pass


# paths
# archive


# git commands
def github_username() -> str:
    """Return the current user's GitHub username as a string

    Requires the username to be configured:

        git config --global user.name

    :return: The username
    """
    return subprocess.check_output(
        ["git", "config", "--global", "user.name"], text=True
    ).strip()


# commit hash utilities
def are_equivalent_commits(hash1: str, hash2: str) -> bool:
    """Check if two git commit hashes are equivalent.

    :param hash1: First commit hash or line
    :param hash2: Second commit hash or line
    :return: `True` if they are, `False` if they aren't
    """
    hash1, hash2 = map(commit_hash_from_line, (hash1, hash2))
    nchars = max(4, min(*map(len, [hash1, hash2])))
    return hash1[:nchars] == hash2[:nchars]


def commit_hash_from_line(line: str) -> str:
    """Get the commit hash from a one-line log summary

    :param line: A one-line log summary
    :return: The commit hash
    """
    return line.split()[0]
