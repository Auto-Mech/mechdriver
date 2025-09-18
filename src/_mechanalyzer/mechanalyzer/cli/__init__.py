import click
import numpy as np
from mechanalyzer.cli import pssa, prompt, sort, ste_mech
from mechanalyzer.cli import compare_rates as compare_rates_
from mechanalyzer.cli import compare_thermo as compare_thermo_
from mechanalyzer.cli import pes_diagram_from_mess as pes_diagram_from_mess_

@click.group()
def main():
    """MechAnalyzer CLI"""
    pass


# sort
@main.command()
@click.option(
    "-m",
    "--mech",
    default="mechanism.dat",
    show_default=True,
    help="Input mechanism file name",
)
@click.option(
    "-s",
    "--spc",
    default="species.csv",
    show_default=True,
    help="Input species file name",
)
@click.option(
    "-t",
    "--therm",
    default="therm.dat",
    show_default=True,
    help="Input thermo file name",
)
@click.option(
    "-i",
    "--sortopts",
    default="sort.dat",
    show_default=True,
    help="Input sort file name",
)
@click.option(
    "-o",
    "--outmech",
    default="outmech.dat",
    show_default=True,
    help="Output mechanism file name",
)
@click.option(
    "-c",
    "--outspc",
    default="outspc.csv",
    show_default=True,
    help="Output species file name",
)
@click.option(
    "-g",
    "--outgroups",
    default="pes_groups.dat",
    show_default=True,
    help="Output PES groups file name",
)
def sortmech(
    mech: str = "mechanism.dat",
    spc: str = "species.csv",
    therm: str = "therm.dat",
    sortopts: str = "sort.dat",
    outmech: str = "outmech.dat",
    outspc: str = "outspc.csv",
    outgroups: str = "pes_groups.dat",
    ):
    """Sort the reactions in a mechanism"""
    sort.main(
        mech=mech,
        spc=spc,
        therm=therm,
        sortopts=sortopts,
        outmech=outmech,
        outspc=outspc,
        outgroups=outgroups,
    )

# compare-rates
@main.command()
@click.option(
    "-m",
    "--mechs_yaml",
    default="mechs.yaml",
    show_default=True,
    help="YAML filename for mechanism specification. Format:\\mech1:\\\trate_file: 'rate.ckin'\\\ttherm_file: 'therm.ckin'\\\tspecies_csv: 'species.csv'\\mech2: ..."
)
@click.option(
    "-o",
    "--plot_fname",
    default="rate_plots.pdf",
    show_default=True,
    help="PDF filename for plot output."
)
@click.option(
    "-f",
    "--out_txt_fname",
    default="rate_ordering.txt",
    show_default=True,
    help="Text filename for sorting output."
)
@click.option(
    "-d",
    "--job_path",
    default="",
    show_default=True,
    help="Directory for input/output files."
)
@click.option(
    "-t",
    "--temps_lst",
    default=None,
    type=lambda s: [float(temp) for temp in s.split(',')],
    show_default=True,
    help="Comma separated first and last temperatures (K). None defaults to 500,1500."
)
@click.option(
    "-p",
    "--pressures",
    default=None,
    type=lambda s: [float(press) for press in s.split(',')],
    show_default=True,
    help="Comma separated array of pressures (atm). None defaults to '1,10,100'."
)
@click.option(
    "-s",
    "--sort_method",
    default="ratios",
    show_default=True,
    help="Sort the plots by the ratio of the differences with 'ratios', or not at all with None."
)
@click.option(
    "-r",
    "--rev_rates",
    default=True,
    show_default=True,
    help="If True, reverses rates of remaining mechanisms to make them match the direction in the first mechanism."
)
@click.option(
    "-l",
    "--remove_loners",
    default=1,
    show_default=True,
    help="1 only plots rates that are in 2 mechanisms, 2 only plots rates that are in ALL mechanisms, 0 plots ALL rates in all mechanisms."
)
def compare_rates(
    mechs_yaml: str,
    plot_fname: str,
    out_txt_fname: str,
    job_path: str,
    temps_lst: list,
    pressures: list,
    sort_method: str,
    rev_rates: bool,
    remove_loners: int
):
    """Compare the rate constants in a mechanism"""
    compare_rates_.main(
        mechs_yaml,
        plot_fname,
        out_txt_fname,
        job_path,
        temps_lst,
        pressures,
        sort_method,
        rev_rates,
        remove_loners,
    )

# compare-thermo
@main.command()
@click.option(
   "-m",
   "--mechs_yaml",
   default="mechs.yaml",
   show_default=True,
   help="YAML filename for mechanism specification. Format:\\mech1:\\\ttherm_file: 'therm.ckin'\\\tspecies_csv: 'species.csv'\\mech2: ...",
)
@click.option(
   "-o",
   "--plot_fname",
   default="thermo_plots.pdf",
   show_default=True,
   help="PDF filename for plot output."
)
@click.option(
   "-f",
   "--out_txt_fname",
   default="thermo_ordering.txt",
   show_default=True,
   help="Text filename for sorting output."
)
@click.option(
   "-d",
   "--job_path",
   default=".",
   show_default=True,
   help="Directory for input/output files."
)
@click.option(
   "-t",
   "--temps_lst",
   default=None,
   type=lambda s: [float(temp) for temp in s.split(',')],
   show_default=True,
   help="Array of temperatures (K). None defaults to 500--1500."
)
@click.option(
   "-s",
   "--sort_method",
   default="lnq",
   show_default=True,
   help=("Method for sorting. Sorts by max difference in enthalpy ('h'), \\"
         "entropy ('s'), Gibbs ('g'), c_p ('cp'), natural log of the \\"
         "partition function ('lnq'), or not at all (None)."
    )
)
@click.option(
   "-st",
   "--sort_temp",
   default=None,
   show_default=True,
   help="If sorting, specifies the temperature at which to sort (K). None (default) specifies to sort by maximum difference."
)
@click.option(
   "-l",
   "--remove_loners",
   default=True,
   show_default=True,
   help="1 only plots species that are in at least 2 mechanisms, 2 only plots species that are in ALL mechanisms, and 0 plots ALL species in all mechanisms."
)
@click.option(
   "-p",
   "--print_missing",
   default=True,
   show_default=True,
   help="True prints a warning for any species that are not in the species.csv file."
)
def compare_thermo(
   mechs_yaml: str,
   plot_fname: str,
   out_txt_fname: str,
   job_path: str,
   temps_lst: list,
   sort_method: str,
   sort_temp: float,
   remove_loners: int,
   print_missing: bool
):
   """Compare the thermo properties in a mechanism"""
   compare_thermo_.main(
       mechs_yaml,
       plot_fname,
       out_txt_fname,
       job_path,
       temps_lst,
       sort_method,
       sort_temp,
       remove_loners,
       print_missing
   )

# expand
@main.command()
def expand():
    """Expand stereochemistry for a mechanism"""
    ste_mech.main()


@main.command()
@click.option(
    "-f",
    "--flds",
    default= ['',],
    show_default=True,
    required=True,
    multiple=True,
    help="directories with mess files; pass each directory as a single -f argument",
)
@click.option(
    "-i",
    "--messinput",
    default="mess.inp",
    show_default=True,
    help="mess input",
)
@click.option(
    "-o",
    "--messoutput",
    default="mess.out",
    show_default=True,
    help="mess rate output",
)
@click.option(
    "-ke",
    "--outputmicro",
    default="ke.out",
    show_default=True,
    help="mess microcanonical rate output",
)
@click.option(
    "-l",
    "--log",
    default="mess.log",
    show_default=True,
    help="mess log file",
)
@click.option(
    "-m",
    "--model",
    default="rovib_dos",
    show_default=True,
    help="Model to compute prompt branching fractions",
)
@click.option(
    "-b",
    "--bfthresh",
    default=0.1,
    show_default=True,
    help="keep reactions if contributing more than this threshold",
)
@click.option(
    "-fit",
    "--fitmethod",
    default="plog",
    show_default=True,
    help="Fitting method",
)
@click.option(
    "-or",
    "--outputrates",
    default="rates_prompt.txt",
    show_default=True,
    help="Output prompt rates file name",
)
def promptcalc(
    flds: list = ['',],
    messinput: str = 'mess.inp',
    messoutput: str = 'mess.out',
    outputmicro: str = 'ke.out',
    log: str = 'mess.log',
    model: str = 'rovib_dos',
    bfthresh: float = 0.1,
    fitmethod: str = 'plog',
    outputrates: str = 'rates_prompt.txt'
):
    """Compute prompt effects from mess output files"""
    prompt.main(
        flds=flds,
        messinput=messinput,
        messoutput=messoutput,
        outputmicro=outputmicro,
        log=log,
        model=model,
        bfthresh=bfthresh,
        fitmethod=fitmethod,
        outputrates=outputrates
    )

@main.command()
@click.option(
    "-i",
    "--startmech",
    default='kin.CKI',
    show_default=True,
    required=True,
    help="input mechanism in chemkin format",
)
@click.option(
    "-s",
    "--pssa_spcs",
    default=['',],
    show_default=True,
    required=True,
    multiple=True,
    help="list of species to apply pssa to; provide each -s argument separately",
)
@click.option(
    "-T",
    "--temprange",
    default=[500, 2000],
    show_default=True,
    multiple=True,
    help="temperature range for rate constant calculation; provide 2 limits as -T arguments separately",
)
@click.option(
    "-P",
    "--pvect",
    default=[1.0,],
    show_default=True,
    multiple=True,
    type=float,
    help="pressure vector for rate constant calculation; provide pressure values as multiple -P arguments separately",
)
@click.option(
    "-b",
    "--bfthresh",
    default=1e-4,
    show_default=True,
    help="keep reactions if contributing more than this threshold",
)
@click.option(
    "-t",
    "--thermofile",
    default=None,
    show_default=True,
    help="thermochemistry file to compute backward rate constants",
)
@click.option(
    "-or",
    "--outputrates",
    default="rates_pssa.txt",
    show_default=True,
    help="Output pssa rates file name",
)
@click.option(
    "-tol",
    "--fitduptol",
    default=15.,
    show_default=True,
    help="% Tolerance to switch from single to double arrhenius fit",
)
def runpssa(
    startmech: str='kin.CKI',
    pssa_spcs: list=['',],
    temprange: list = [500, 2000],
    pvect: list =[1.0,],
    bfthresh: float = 1e-4,
    thermofile: str = None,
    outputrates: str = 'pssa_rates.txt',
    fitduptol: float = 15.,
):
    """Sort the reactions in a mechanism"""
    pssa.main(
    startmech=startmech,
    pssa_spcs=pssa_spcs,
    temprange=temprange,
    pvect=pvect,
    bfthresh=bfthresh,
    thermofile=thermofile,
    outputrates=outputrates,
    fitduptol=fitduptol,
    )

# PES diagram
@main.command()
@click.option(
    "--input_file",
    "-i",
    type=str,
    help="Path to the MESS input file",
    show_default=True,
    default="mess.inp")
@click.option(
    "--well_threshold",
    "-w",
    type=int,
    help="How many connections a species must have to be centered as well",
    show_default=True,
    default=2)
@click.option(
    "--colors_on",
    "-c",
    type=bool,
    help="True/False colorful PES, automatically makes each well and their connections a unique color",
    show_default=True,
    default=True)
@click.option(
    "--gravity",
    "-g",
    type=int,
    help="How much the species are pulled together in the spring layout",
    show_default=True,
    default=1)
@click.option(
    "--spring_iterations",
    "-s",
    type=int,
    help="How many iterations to run the spring layout algorithm",
    show_default=True,
    default=20000)
@click.option(
    "--nudge_iterations",
    "-n",
    type=int,
    show_default=True,
    help="How many iterations to nudge the species to minimize overlap",
    default=10)
@click.option(
    "--min_distance",
    "-d",
    type=float,
    help="Minimum distance between non-neighboring nodes, aka whats considered overlap",
    show_default=True,
    default=1.0)
@click.option(
    "--labels",
    "-l",
    type=bool,
    help="True/False whether to label the species in the PES",
    show_default=True,
    default=True)
@click.option(
    "--output_file",
    "-o",
    type=str,
    help="Name of the output figure file",
    show_default=True,
    default="pes_diagram")
@click.option(
    "--format",
    "-f",
    type=str,
    help="Format of the output figure file (e.g., svg, png)",
    show_default=True,
    default="svg")
@click.option(
    "--aspect_ratio",
    "-a",
    type=float,
    help="Aspect ratio of the output figure (width / height)",
    show_default=True,
    default=1.0)
@click.option(
    "--remove_fake",
    "-r",
    type=bool,
    help="remove fake vDW",
    show_default=True,
    default=True)
@click.option(
    "--shift_energy",
    "-e",
    type=bool,
    help="Whether to shift energies relative to the lowest well",
    show_default=True,
    default=True)


def pes_diagram(
    input_file: str = "mess.inp",
    well_threshold: int = 2,
    colors_on: bool = True,
    gravity: int = 1,
    spring_iterations: int = 20000,
    nudge_iterations: int = 10,
    min_distance: float = 1.0,
    max_distance: float = 30.0,
    output_file: str = "pes_diagram",
    format: str = "svg",
    aspect_ratio: float = 1.0,
    labels: bool = True,
    remove_fake: bool = True,
    shift_energy: bool = True
):
    """Generate a PES diagram from a MESS input file"""

    pes_diagram_from_mess_.main(
        input_file,
        well_threshold,
        colors_on,
        gravity,
        spring_iterations,
        nudge_iterations,
        min_distance,
        max_distance,
        output_file,
        format,
        aspect_ratio,
        labels,
        remove_fake,
        shift_energy)
