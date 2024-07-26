import click

from automech.cli import _check_log, _run, _subtasks_setup


@click.group()
def main():
    """AutoMech CLI"""
    pass


@main.command()
@click.option(
    "-p", "--path", default=".", show_default=True, help="The job run directory"
)
@click.option("-S", "--safemode-off", is_flag=True, help="Turn off safemode?")
def run(path: str = ".", safemode_off: bool = False):
    """Run central workflow
    Central Execution script to launch a MechDriver process which will
    parse all of the user-supplied input files in a specified directory, then
    launches all of the requested electronic structure, transport,
    thermochemistry and kinetics calculations via their associated
    sub-drivers.

    The AutoMech directory must contain an `inp/` subdirectory with the following
    required files: run.dat, theory.dat, models.dat, species.csv, mechanism.dat
    """
    _run.main(path=path, safemode_off=safemode_off)


@main.command()
@click.option(
    "-p", "--path", default=".", show_default=True, help="The path to the log file"
)
def check_log(path: str = "."):
    """Check an AutoMech log file to see if it succeeded

    The path must point either directly to the log file, or to a directory where the log
    file is named "out.log"
    """
    _check_log.main(path=path)


@main.group()
def subtasks():
    """Run AutoMech subtasks in parallel"""
    pass


@subtasks.command()
@click.option(
    "-p", "--path", default=".", show_default=True, help="The job run directory"
)
def setup(path: str = "."):
    """Set-up subtasks from a user-supplied AutoMech directory

    The AutoMech directory must contain an `inp/` subdirectory with the following
    required files: run.dat, theory.dat, models.dat, species.csv, mechanism.dat
    """
    _subtasks_setup.main(path=path)


@main.command()
def greetme():
    """Hello world function, for CLI testing purposes"""
    print("Hello, world!")
