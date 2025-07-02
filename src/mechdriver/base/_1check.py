"""Check an AutoMech log file for errors"""

import enum
import re
import time
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


class Extension:
    running = ".running"


def check_log(path: str | Path = ".", log: bool = False) -> tuple[Status, str | None]:
    """Check an AutoMech log file to see if it succeeded

    :param path: The path to the log file or directory. If the path is a directory, the
        log file must be called `out.log`.
    :param log: Whether to print the result to the terminal.
    """
    path = Path(path)
    assert path.exists(), f"Path does not exist: {path}"
    if path.is_dir():
        path /= "out.log"
    assert path.is_file(), f"File does not exist: {path}"

    status, line = _check_log(path)

    if log:
        print(f"{str(path) + ' ':.<80} {colored_status_string(status)}")
        if line is not None:
            print(line)

    return (status, line)


def _check_log(log_path: str | Path) -> tuple[Status, str | None]:
    """Check a log file, returning the status and the line that triggered it.

    :param log_path: The log file path
    :return: The status and the line triggering the status, if applicable
    """
    log_path = Path(log_path)
    lock_path = log_path.with_suffix(Extension.running)
    line = None
    if not log_path.exists():
        status = Status.TBD
        return (status, line)

    log = log_path.read_text().strip()
    has_exit_message = re.search("EXITING AUTOMECHANIC", log)
    has_lock_file = lock_path.is_file()
    has_recent_log = (time.time() - log_path.stat().st_mtime) < 60
    error_match = re.search("ERROR", log)
    has_recent_log_no_error = has_recent_log and not error_match
    if not has_exit_message:
        status = (
            Status.RUNNING if has_lock_file or has_recent_log_no_error else Status.ERROR
        )
        line = log.splitlines()[-1] if log else ""
        return (status, line)

    warning_match = re.search(r".*(?<!Future)Warning.*", log, flags=re.IGNORECASE)
    status = Status.WARNING if warning_match else Status.OK
    status = Status.ERROR if error_match else status
    line = warning_match.group(0) if warning_match else None
    line = error_match.group(0) if error_match else line
    return (status, line)


def colored_status_string(status: Status, width: int | None = None) -> str:
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
    width = max(len(s.value) for s in Status) if width is None else width
    return f"{color_start_code}{status.value:^{width}}{color_end_code}"
