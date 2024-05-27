import click
from automech.cli import automech


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
    automech.main(path=path, safemode_off=safemode_off)


@main.command()

def greetme():
    """Hello world function, for CLI testing purposes"""
    print("Hello, world!")
