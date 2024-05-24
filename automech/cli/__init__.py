import click
from automech.cli import automech


@click.group()
def main():
    """AutoMech CLI"""
    pass


@main.command()
def run():
    """Expand stereochemistry for a mechanism
    """
    automech.main()


@main.command()
def greetme():
    """Hello world function, for CLI testing purposes"""
    print("Hello, world!")
