"""End-to-end AutoMech tests
"""

import contextlib
import os
import sys
from pathlib import Path

import automech
import pytest
import yaml

DIR = Path(__file__).parent
TESTS = [t["name"] for t in yaml.safe_load((DIR / "config.yaml").read_text())]


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


@pytest.mark.parametrize("name", TESTS)
def test_workflow(name: str):
    """Test the entire workflow"""
    print(f"Running in {name}...")

    test_dir = DIR / name
    os.chdir(test_dir)

    with contextlib.redirect_stdout(Logger("out.log")):
        automech.subtasks.untar_subtask_data()
        automech.run()


if __name__ == "__main__":
    test_workflow("quick")
