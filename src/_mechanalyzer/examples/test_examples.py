import subprocess
from pathlib import Path

import pytest

EXAMPLE_DIR = Path(__file__).parent
EXAMPLE_NAMES = [
    "sort_basic",
]


@pytest.mark.parametrize("name", EXAMPLE_NAMES)
def test__examples(name):
    path = EXAMPLE_DIR / name
    out = subprocess.run(
        "./command.sh", shell=True, cwd=path, capture_output=True, text=True, check=True
    )
    print(out.stdout)
