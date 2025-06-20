"""
  Library of Runtime Messages
"""

import random
from mechlib.amech_io.printer._lib import obj
from mechlib.amech_io.printer._print import message


def program_header(driver):
    """ print the header for a program
    """
    header_dct = {
        'amech': AMECH_MSG,
        'inp': INP_MSG,
        'ktp': KTP_MSG,
        'thermo': THM_MSG,
        'trans': TRANS_MSG,
        'es': ES_MSG,
        'proc': PRINT_MSG
    }
    message(header_dct[driver])


def program_exit(driver):
    """ print the header for a program
    """
    header_dct = {
        'amech': AMECH_EXIT_MSG,
        'inp': INP_EXIT_MSG,
        'ktp': KTP_EXIT_MSG,
        'thermo': THM_EXIT_MSG,
        'trans': TRANS_EXIT_MSG,
        'es': ES_EXIT_MSG,
        'proc': PRINT_EXIT_MSG
    }
    message(header_dct[driver]+'\n')


def driver_tasks(
        run_es, write_messpf, run_messpf, run_nasa,
        write_messrate, run_messrate, run_fits, run_trans):
    """ a
    """

    if run_es:
        message('  - ESDriver')
        # Add the tasks for the ESDriver
    if write_messpf or run_messpf or run_nasa:
        message('  - ThermoDriver')
        if write_messpf:
            message('    - write_messpf')
        if run_messpf:
            message('    - run_messpf')
        if run_nasa:
            message('    - run_nasa')
    if write_messrate or run_messrate or run_fits:
        message('  - kTPDriver')
        if write_messrate:
            message('    - write_messrate')
        if run_messrate:
            message('    - run_messrate')
        if run_fits:
            message('    - run_fits')
    if run_trans:
        message('  - TransportDriver')
        # Add the tasks for the ESDriver


# LIBRARY OF MESSAGES
AMECH_MSG = """
          ================================================================
          ==                        AUTOMECHANIC                        ==
          ===         Andreas Copan, Sarah Elliott, Kevin Moore,       ===
          ===            Luna Pratali Maffei, Daniel Moberg,           ===
          ===            Carlo Cavallotti, Yuri Georgievski,           ===
          ==       Ahren Jasper, Murat Keceli, Stephen Klippenstein     ==
          ================================================================"""

AMECH_EXIT_MSG = """
          ================================================================
          ==                    EXITING AUTOMECHANIC                    ==
          ================================================================"""

INP_MSG = """
=========================================================
PARSING INPUT
========================================================="""

INP_EXIT_MSG = """
=========================================================
FINISHED INPUT PARSING
========================================================="""

KTP_MSG = """
=========================================================
kTPDRIVER

Sarah Elliott, Kevin Moore, Andreas Copan,
Daniel Moberg, Carlo Cavallotti, Yuri Georgievski,
Ahren Jasper, Stephen Klippenstein
========================================================="""

KTP_EXIT_MSG = """
=========================================================
EXITING kTPDRIVER
========================================================="""

THM_MSG = """
=========================================================
THERMODRIVER

Sarah Elliott, Kevin Moore, Andreas Copan,
Murat Keceli, Yuri Georgievski, Stephen Klippenstein
========================================================="""

THM_EXIT_MSG = """
=========================================================
EXITING THERMODRIVER
========================================================="""

TRANS_MSG = """
=========================================================
TRANSPORTDRIVER
Kevin Moore, Ahren Jasper, Stephen Klippenstein
========================================================="""

TRANS_EXIT_MSG = """
=========================================================
EXITING TRANSPORTDRIVER
========================================================="""

ES_MSG = """
=========================================================
ESDRIVER
Sarah Elliott, Andreas Copan, Kevin Moore,
Carlo Cavolotti, Stephen Klippenstein
========================================================="""

ES_EXIT_MSG = """
=========================================================
EXITING ESDRIVER
========================================================="""

PRINT_MSG = """
=========================================================
OUTPUT PREPARATION
Sarah Elliott, Andreas Copan, Kevin Moore,
Carlo Cavolotti, Stephen Klippenstein
========================================================="""

PRINT_EXIT_MSG = """
=========================================================
EXITING OUTPUT PREP
========================================================="""


def random_cute_animal():
    """ Print a picture of a fun, cute animal at random
    """
    msg = random.choice([
        r"""
                                        _,--._
                                      ,'      `.
                              |\     /          \     /|
                              )o),/ ( ,--,  ,--, ) \.(o(
                             /o/// /|            |\ \\ \\o\\
                            / / |\ \(   .----,   )/ /| \ \\
                            | | \o`-/    `--'    \-'o/ | |
                            \ \  `,'              `.'  / /
                         \.  \ `-'  ,'|   /\   |`.  `-' /  ,/
                          \`. `.__,' /   /  \   \ `.__,' ,'/
                           \o\     ,'  ,'    `.  `.     /o/
                            \o`---'  ,'        `.  `---'o/
                             `.____,'           `.____,'  """,

        r"""
                                 ,,,         ,,,
                               ;"   ^;     ;'   ",
                              ;    s$$$$$$$s      ;
                               ,  ss$$$$$$$$$$s  ,'
                               ;s$$$$$$$$$$$$$$$
                               $$$$$$$$$$$$$$$$$$
                              $$$$P""Y$$$Y""W$$$$$
                              $$$$  0"$$$"0  $$$$$
                              $$$$  .$$$$$.  $$$$
                               $$$$$$$$$$$$$$$$$
                                "Y$$$"'*'"$$$Y"
                                   "$$b.d$$"        """,
        r"""
                                   _.---~-~-~~-..
               ..       __.    .-~               ~-.
               ((\     /   `}.~                     `.
                \\\\\   {     }              /     \   \\
            (\   \\\\~~       }             |       }   \\
             \`.-~-^~     }  ,-,.         |       )    \\
             (___,    ) _}  (    :        |    / /      `.
              `----._-~.     _\ \ |_       \   / /- _     -.
                     ~~----~~  \ \| ~~--~~~(  + /     ~-.   '--~.
                               /  /         \  \         `~-..__ `~__
                            __/  /          _\  )               `~~---'
                          .<___.'         .<___/  """])
    message(msg)
    obj('vspace')
