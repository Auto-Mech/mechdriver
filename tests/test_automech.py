"""End-to-end AutoMech tests
"""

import contextlib
import os
import sys
import tarfile
from pathlib import Path

import automech
import pytest
import yaml

TEST_ROOT_DIR = Path(__file__).parent
TESTS = [t for t in yaml.safe_load((TEST_ROOT_DIR / "config.yaml").read_text())]
ARCHIVE_FILE = TEST_ROOT_DIR / "archive.tgz"


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


def test_sign():
    """Unpack archive and check provenance"""
    os.chdir(TEST_ROOT_DIR)
    if ARCHIVE_FILE.exists():
        print(f"Unpacking {ARCHIVE_FILE}...")
        with tarfile.open(ARCHIVE_FILE, "r") as tar:
            tar.extractall()


@pytest.mark.parametrize("name", TESTS)
def test_workflow(name: str):
    """Test the entire workflow"""
    print(f"Running in {name}...")

    test_dir = TEST_ROOT_DIR / name
    os.chdir(test_dir)

    with contextlib.redirect_stdout(Logger("out.log")):
        automech.run()


if __name__ == "__main__":
    # test_workflow("quick")
    test_sign()
