"""
  Library of Runtime Messages
"""

import random


def program_header(driver):
    """ print the header for a program
    """
    header_dct = {
        'amech': AMECH_MSG,
        'inp': INP_MSG,
        'ktp': KTP_MSG,
        'thermo': THM_MSG,
        'es': ES_MSG
    }
    print(header_dct[driver])


def program_exit(driver):
    """ print the header for a program
    """
    header_dct = {
        'amech': AMECH_EXIT_MSG,
        'inp': INP_EXIT_MSG,
        'ktp': KTP_EXIT_MSG,
        'thermo': THM_EXIT_MSG,
        'es': ES_EXIT_MSG
    }
    print(header_dct[driver]+'\n')


AMECH_MSG = """
          ================================================================
          ==                        AUTOMECHANIC                        ==
          ===         Andreas Copan, Sarah Elliott, Kevin Moore,       ===
          ===     Daniel Moberg, Carlo Cavallotti, Yuri Georgievski,   ===
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
    print(msg)
    print('\n\n')
