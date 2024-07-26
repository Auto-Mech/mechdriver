""" Check an AutoMech log file for errors
"""

import enum
import re
from pathlib import Path


class Status(enum.Enum):
    SUCCESS = "SUCCESS"
    WARNING = "WARNING"
    ERROR = "ERROR"


def print_status(status: Status, path: str | Path) -> None:
    path = str(path)
    color_start_code = {
        Status.SUCCESS: "\033[92m", # green
        Status.WARNING: "\033[93m", # yellow
        Status.ERROR: "\033[91m",   # red
    }.get(status)
    color_end_code = "\033[0m"
    print(f"{path:.<80} {color_start_code}{status.value:^7}{color_end_code}")


def main(path: str = "."):
    """Check an AutoMech log file to see if it succeeded

    :param path: The path to the log file or directory. If the path is a directory, the
        log file must be called `out.log`.
    """
    path: Path = Path(path)
    assert path.exists(), f"Path does not exist: {path}"
    if path.is_dir():
        path /= "out.log"
    assert path.is_file(), f"File does not exist: {path}"

    log = path.read_text()
    has_exit_message = re.search("EXITING AUTOMECHANIC", log)
    has_warning = re.search("(?<!Future)Warning", log, flags=re.IGNORECASE)

    status = (
        Status.ERROR
        if not has_exit_message
        else Status.WARNING if has_warning else Status.SUCCESS
    )
    print_status(status, path)
