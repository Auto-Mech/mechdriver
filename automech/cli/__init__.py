import click
from automech.cli import _run
from automech.cli import _setup_subtasks


@click.group()
def main():
    """AutoMech CLI"""
    pass


@main.command()
@click.option('-p', '--path', default='.', show_default=True, help='The job run directory')
@click.option('-S', '--safemode-off', is_flag=True, help='Turn off safemode?')
def run(path: str = ".", safemode_off: bool = False):
    """Central Execution script to launch a MechDriver process which will
    parse all of the user-supplied input files in a specified directory, then
    launches all of the requested electronic structure, transport,
    thermochemistry and kinetics calculations via their associated
    sub-drivers.

    The job run directory must contain an `inp/` subdirectory with the following
    required files: run.dat, theory.dat, models.dat, species.csv
    """
    _run.main(path=path, safemode_off=safemode_off)


@main.command()
@click.option('-p', '--path', default='.', show_default=True, help='The job run directory')
def setup_subtasks(path: str = "."):
    """Take the user-supplied input files in a job run directory and split them into
    subtasks for parallel execution

    The job run directory must contain an `inp/` subdirectory with the following
    required files: run.dat, theory.dat, models.dat, species.csv
    """
    _setup_subtasks.main(path=path)


@main.command()

def greetme():
    """Hello world function, for CLI testing purposes"""
    print("Hello, world!")
