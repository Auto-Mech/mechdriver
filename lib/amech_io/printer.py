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


AMECH_MSG = """
          ================================================================
          ==                        AUTOMECHANIC                        ==
          ===         Andreas Copan, Sarah Elliott, Kevin Moore,       ===
          ===     Daniel Moberg, Carlo Cavallotti, Yuri Georgievski,   ===
          ==       Ahren Jasper, Murat Keceli, Stephen Klippenstein     ==
          ================================================================"""

INP_MSG = """
=========================================================
INPUT PARSER
========================================================="""

KTP_MSG = """
=========================================================
kTPDRIVER

Sarah Elliott, Kevin Moore, Andreas Copan,
Daniel Moberg, Carlo Cavallotti, Yuri Georgievski,
Ahren Jasper, Stephen Klippenstein 
========================================================="""

THM_MSG = """
=========================================================
THERMODRIVER

Sarah Elliott, Kevin Moore, Andreas Copan,
Murat Keceli, Yuri Georgievski, Stephen Klippenstein
========================================================="""

ES_MSG = """
=========================================================
ESDRIVER
Sarah Elliott, Andreas Copan, Kevin Moore,
Carlo Cavolotti, Stephen Klippenstein
========================================================="""


def random_cute_animal():
    msg = random.choice(["""
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

     """
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
    """
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
                          .<___.'         .<___/
    """])
    print(msg)
    print('\n\n')
# _KTP_MSG = """
#           ================================================================
#           ==                          KTPDRIVER                         ==
#           ===         Sarah Elliott, Kevin Moore, Andreas Copan,       ===
#           ===      Daniel Moberg, Carlo Cavallotti, Yuri Georgievski,  ===
#           ==            Ahren Jasper, Stephen Klippenstein              ==
#           ================================================================"""
# 
# _THM_MSG = """
#           ================================================================
#           ==                        THERMODRIVER                        ==
#           ===         Sarah Elliott, Kevin Moore, Andreas Copan,       ===
#           ===    Murat Keceli, Yuri Georgievski, Stephen Klippenstein   ==
#           ================================================================"""
# 
# _ES_MSG = """
#           ================================================================
#           ==                          ESDRIVER                          ==
#           ====        Sarah Elliott, Andreas Copan, Kevin Moore,      ====
#           ==            Carlo Cavolotti, Stephen Klippenstein           ==
#           ================================================================"""
