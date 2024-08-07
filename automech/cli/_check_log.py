""" Check an AutoMech log file for errors
"""

import enum
import re
from pathlib import Path


class Status(enum.Enum):
    # Status codes for one log file that exists
    OK = "OK"
    WARNING = "WARNING"
    ERROR = "ERROR"
    # Status codes for multiple or non-existent log files
    RUNNING = "RUNNING"  # Currently running
    TBD = "TBD"  # Not yet started
    OK_1E = "OK_1E"  # All but 1 log file succeeded
    OK_2E = "OK_2E"  # All but 2 log files succeeded


STATUS_WIDTH = 7


def parse_log_status(log_path: str | Path) -> Status:
    """Parse the status of a log file

    :param log_path: The log file path
    :return: The status
    """
    log_path = Path(log_path)
    if not log_path.exists():
        return Status.TBD

    log = log_path.read_text()
    has_exit_message = re.search("EXITING AUTOMECHANIC", log)
    has_is_running_file = Path(f"{log_path}_IS_RUNNING").exists()
    if not has_exit_message:
        return Status.RUNNING if has_is_running_file else Status.ERROR

    has_warning = re.search("(?<!Future)Warning", log, flags=re.IGNORECASE)
    return Status.WARNING if has_warning else Status.OK


def colored_status_string(status: Status) -> str:
    """Get a colored status string

    :param status: The status
    :return: The colored status string
    """
    color_start_code = {
        Status.OK: "\033[92m",  # bright green
        Status.WARNING: "\033[93m",  # bright yellow
        Status.ERROR: "\033[91m",  # bright red
        Status.RUNNING: "\033[96m",  # bright cyan
        Status.TBD: "\033[90m",  # gray
        Status.OK_1E: "\033[95m",  # bright magenta
        Status.OK_2E: "\033[95m",  # bright magenta
    }.get(status)
    color_end_code = "\033[0m"
    return f"{color_start_code}{status.value:^{STATUS_WIDTH}}{color_end_code}"


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

    status = parse_log_status(path)
    print(f"{str(path) + ' ':.<80} {colored_status_string(status)}")
