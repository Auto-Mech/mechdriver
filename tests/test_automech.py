"""End-to-end AutoMech tests
"""

import contextlib
import os
import subprocess
import sys
import tarfile
from pathlib import Path

import automech
import pytest
import yaml

TEST_ROOT_DIR = Path(__file__).parent
TESTS = [t for t in yaml.safe_load((TEST_ROOT_DIR / "config.yaml").read_text())]
ARCHIVE_FILE = TEST_ROOT_DIR / "archive.tgz"
SIGNATURE_FILE = TEST_ROOT_DIR / "signature.yaml"


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


class InvalidSignatureError(Exception):
    def __init__(self, message, *_):
        # Call the base class constructor with the parameters it needs
        lines = [
            message,
            "  You need to sign off on your local tests using `pixi run test sign`",
            "  (This must be done for an up-to-date version of the code)",
        ]
        super().__init__("\n".join(lines))


def test_signature():
    """Unpack archive and check provenance"""
    os.chdir(TEST_ROOT_DIR)

    # 1. Unpack the test archive
    unpack_archive()

    # 2. Read test provenance from signature file
    sign_dct = yaml.safe_load(SIGNATURE_FILE.read_text())
    prov_dct: dict[str, str] = sign_dct["provenance"]

    # 3. Check the commit hash of this repo
    # (must match most recent commit that isn't a merge PR or test archive update)
    print("Checking whether mechdriver matches the signed version...")
    exclude_messages = ("Merge pull request", "Updates tests/archive.tgz")
    exclude_grep = r"\|".join(exclude_messages)
    self_repo_prov = subprocess.check_output(
        ["git", "log", "--oneline", "-1", "--invert-grep", "--grep", exclude_grep],
        text=True,
    )
    self_sign_prov = prov_dct.pop("mechdriver")
    self_repo_commit = self_repo_prov.split()[0]
    self_sign_commit = self_sign_prov.split()[0]
    if not are_equivalent_commit_hashes(self_repo_commit, self_sign_commit):
        raise InvalidSignatureError(f"{self_repo_commit} !~ {self_sign_commit}")

    # 4. Check the commit hashes of the other repos
    # (must match most recent commit)
    for repo, sign_prov in prov_dct.items():
        print(f"Checking whether {repo} matches the signed version...")

        repo_prov = subprocess.check_output(
            ["git", "ls-remote", f"https://github.com/Auto-Mech/{repo}.git", "HEAD"],
            text=True,
        )
        repo_commit = repo_prov.split()[0]
        sign_commit = sign_prov.split()[0]
        if not are_equivalent_commit_hashes(repo_commit, sign_commit):
            raise InvalidSignatureError(f"{repo_commit} !~ {sign_commit}")


@pytest.mark.parametrize("name", TESTS)
def test_workflow(name: str):
    """Test the entire workflow"""
    print(f"Running in {name}...")

    test_dir = TEST_ROOT_DIR / name
    os.chdir(test_dir)

    with contextlib.redirect_stdout(Logger("out.log")):
        automech.run()


def unpack_archive() -> None:
    """Un-tar a directory in its current location

    :param path: The path where the AutoMech subtasks were set up
    """
    os.chdir(TEST_ROOT_DIR)
    if ARCHIVE_FILE.exists():
        print(f"Unpacking {ARCHIVE_FILE}...")
        with tarfile.open(ARCHIVE_FILE, "r") as tar:
            tar.extractall()


def are_equivalent_commit_hashes(hash1: str, hash2: str) -> bool:
    """Check if two git commit hashes are equivalent.

    :param hash1: First hash
    :param hash2: Second hash
    :return: `True` if they are, `False` if they aren't
    """
    nchars = max(4, min(*map(len, [hash1, hash2])))
    return hash1[:nchars] == hash2[:nchars]


if __name__ == "__main__":
    # test_workflow("quick")
    test_signature()
